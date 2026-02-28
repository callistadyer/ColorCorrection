import torch, numpy as np
from abc import ABC, abstractmethod
import warnings

# base class for measurement matrix
class Measurement(ABC):
    @abstractmethod
    def measure(self, x):
        pass

    @abstractmethod
    def recon(self, x):
        pass

# simple measurement matrix with convolution
# use ConvTranspose2d and Conv2d for projection
class ConvMatrix(Measurement):
    def __init__(self, kernel_size, stride, device, channels=3):
        # use to record image size
        self.imsize = None

        # conv_transpose work for stride >= kernel_size
        if stride < kernel_size:
            raise Warning(f'''Stride ({stride}) should NOT be smaller than kernel size ({kernel_size})
                                for ConvTranspose2d to act properly as a linear projection''')

        # sampling kernel
        self.conv = torch.nn.Conv2d(channels, channels, kernel_size, stride,
                                    groups=channels, bias=False, device=device)

        # inverse operation
        self.conv_tr = torch.nn.ConvTranspose2d(channels, channels, kernel_size, stride,
                                                groups=channels, bias=False, device=device)

        # make an averaging conv kernal
        kernel = torch.ones(channels, 1, kernel_size, kernel_size, device=device)
        kernel = torch.nn.Parameter(kernel / torch.norm(kernel[0, 0, :, :]))

        self.conv.weight = kernel
        self.conv_tr.weight = kernel

    # assume calculation requires no gradient
    # linear measurement
    def measure(self, x):
        with torch.no_grad():
            self.imsize = x.shape[1:]
            x = self.conv(x.unsqueeze(0))
            return x.squeeze(0)

    # linear projection
    def recon(self, x):
        with torch.no_grad():
            x = x.unsqueeze(0)
            x = self.conv_tr(x, output_size=self.imsize)
            return x.squeeze(0)

class RenderMatrix(Measurement):
    def __init__(self, R, im_size, device):
        self.im_size = im_size
        #
        self.R = R.to(device)

    def to(self, device):
        self.R = self.R.to(device)
        return self

    def measure(self, x):
        '''
        Given (orthogonalized) render matrix R
        and image x, compute the measurement

        A transpose is required due to different
        in convention between MATLAB and torch
        (Row-Major vs Column-Major)
        '''

        # transpose such that flatten() result is column-major (MATLAB/Fortran)
        return torch.matmul(self.R, x.transpose(1, 2).flatten())

    def recon(self, msmt):
        '''
        From measurement to image space
        '''

        # transpose such that flatten() result is column-major (MATLAB/Fortran)
        return torch.matmul(self.R.T, msmt).reshape(self.im_size).transpose(1, 2)


class DichromatMatrix(Measurement):
    """
    For each pixel i we apply the a 2x3 matrix A to that pixel's RGB (Converts RGB -> LMS -> dichromat)
    That is: y_i(2x1) = A(2x3) @ x_i(3x1).

    If you stacked all pixels into a single long vector x_big of shape (3N,),
    then you would need to build that large matrix R_big: (2N, 3N)
    but this is unweildy so this is an attempt to do it another way:
    """

    def __init__(self, A_2x3, im_size, device):
        """
        Inputs
        ----------
        A_2x3 : torch.Tensor
            Shape (2,3). Per-pixel measurement matrix
            Maps one pixel RGB (3,) -> 2D measurement (2,)
        im_size : tuple
            Expected image size in channel-first format: (3, H, W)
        device : torch.device
            Device where this should live (cpu/cuda/mps)
        """
        self.im_size = im_size                 # stores (3, H, W) for reshaping in recon()
        self.device  = device                  # stores current device the operator is on
        self.A       = A_2x3.to(device).float()# move A to appropriate device, cast to float; shape (2,3)
        self.AT      = self.A.T                # transpose of A; shape (3,2)

    def to(self, device):
        """
        Move A and AT to appropriate device
        Weird stuff is happening where different variables are "living" on different devices so this
        is an attempt to make sure they are all on the same one
        """
        self.device = device                   # remember device
        self.A      = self.A.to(device)        # move A; still shape (2,3)
        self.AT     = self.A.T                 # recompute transpose on same device; shape (3,2)
        return self                           
    
    def measure(self, x):
        """
        Apply the measurement operator

        Inputs
        -----
        x : torch.Tensor
            Image tensor of shape (3, H, W)

        Output
        ------
        msmt : torch.Tensor
            Flattened measurement vector of shape (2N,), where N = H*W 
            This corresponds to concatenating each pixel's 2 measurements (2 because of dichromat)
        """

        # If x is on GPU (mps/cuda) but A is on CPU (or vice versa), matmul will error
        # This makes sure A is on x's device 
        if self.A.device != x.device:
            self.to(x.device)

        # Read spatial dimensions from x
        # x.shape is (3, H, W), so:
        H, W = x.shape[1], x.shape[2]          # H, W are ints

        # Reorder dimensions so pixels are rows:
        # x: (3, H, W)  ->  (H, W, 3)
        # Then flatten spatial dims: (H, W, 3) -> (N, 3), where N = H*W
        xN3 = x.permute(1, 2, 0).reshape(-1, 3)# shape (N,3)

        # Apply the same 2x3 matrix A to every pixel:
        # xN3: (N,3)
        # A.T: (3,2)
        # result yN2: (N,2)
        yN2 = xN3 @ self.A.T                   # shape (N,2)

        # Flatten (N,2) into a single vector (2N,)
        # This matches the convention used by RenderMatrix: return a 1D measurement vector
        return yN2.reshape(-1)                 # shape (2N,)

    def recon(self, msmt):
        """
        Apply the reconstruction mapping

        Input
        -----
        msmt : torch.Tensor
            Measurement vector of shape (2N,)

        Output
        ------
        x : torch.Tensor
            Reconstructed??? image of shape (3, H, W)
        """

        # Same device-safety guard as measure():
        # if msmt is on GPU but A is on CPU, reshape+matmul will error later
        if self.A.device != msmt.device:
            self.to(msmt.device)

        # Use stored image size for reshape:
        # self.im_size is (3, H, W)
        H, W = self.im_size[1], self.im_size[2]

        # Undo the flattening of measurements:
        # msmt: (2N,) -> (N,2)
        yN2 = msmt.reshape(-1, 2)               # shape (N,2)

        # Map measurement back to RGB space using A (or A^T)
        # Here we want (N,3):
        #   (N,2) @ (2,3) = (N,3)
        #
        # self.AT is (3,2), so self.A is (2,3)
        xN3 = yN2 @ self.A                   # (N,2) @ (2,3) -> (N,3)

        # Reshape rows back to image grid
        # (N,3) -> (H, W, 3) -> (3, H, W)
        x = xN3.reshape(H, W, 3).permute(2, 0, 1)  # shape (3,H,W)

        return x

class ArrayMatrix(Measurement):
    '''
    Generalization of the RenderMatrix class to an array
    of matrices that tile through larger images
    '''

    def __init__(self, array, array_size, im_size, device):
        self.device = device
        self.nx = array_size[0]
        self.ny = array_size[0]

        # im_size[1] == im_size[2]
        self.im_size = im_size
        self.edge = im_size[1]

        for idx in range(self.nx):
            for idy in range(self.ny):
                array[idx][idy] = torch.tensor(array[idx][idy].astype('float32')).to(device)

        self.array = array

    def measure(self, x):
        '''
        Measurement array from a set of matrices
        '''

        # init
        msmt = [[0 for y in range(self.ny)]
                      for x in range(self.nx)]

        # loop through matrices
        for idx in range(self.nx):
            for idy in range(self.ny):
                sliced = x[:, idy * self.edge : (idy + 1) * self.edge,
                          idx * self.edge : (idx + 1) * self.edge]

                msmt[idx][idy] = \
                    torch.matmul(self.array[idx][idy],
                                 sliced.transpose(1, 2).flatten())
        return msmt

    def recon(self, msmt):
        '''
        From an array of measurements to image space
        '''

        # init
        recon = torch.empty(size=(3, self.ny * self.edge,
                                  self.nx * self.edge),
                            device=self.device)

        # loop through measurements
        for idx in range(self.nx):
            for idy in range(self.ny):
                sliced = torch.matmul(self.array[idx][idy].T, msmt[idx][idy])

                recon[:, idy * self.edge : (idy + 1) * self.edge,
                      idx * self.edge : (idx + 1) * self.edge] = \
                sliced.reshape(self.im_size).transpose(1, 2)

        return recon

# sample from prior with linear constraint (render matrix)
def linear_inverse(model, render, input, h_init=0.01, beta=0.01, sig_end=0.01,
    t_max=float('inf'), stride=10, seed=None, msmt_flag=False, with_grad=False):

    if not (seed is None):
        torch.manual_seed(seed)

    # Callista changed this - think I had a problem with device 
    # device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    if torch.backends.mps.is_available() and torch.backends.mps.is_built():
        device = torch.device("mps")
    elif torch.cuda.is_available():
        device = torch.device("cuda")
    else:
        device = torch.device("cpu")

    # Put the denoiser / score model in evaluation mode and move to device
    model = model.eval().to(device)
    
    if hasattr(render, "to"):
        render = render.to(device)
    if torch.is_tensor(input):
        input = input.to(device)

    # helper function for pytorch image to numpy image
    numpy_image = lambda y: y.detach().cpu() \
                .squeeze(0).permute(1, 2, 0).numpy()

    # the network calculates the noise residual
    # Define log prior gradient function:
    # model(y) predicts noise residual; minus sign gives gradient of log prior
    if with_grad:
        log_grad = lambda y: - model(y)
    else:
        def log_grad(y):
            with torch.no_grad():
                return - model(y)

    # measurement matrix calculation
    # Measurement operator R and backwards R_T
    # R(x)    : image -> measurements
    # R_T(a)  : measurements -> image-shaped backprojection
    R = render.measure
    R_T = render.recon

    # init variables
    if msmt_flag or input.dim() == 1:
        proj = R_T(input)
    elif input.dim() == 3:
        proj = R_T(R(input))

    # Create an all-ones image with same shape as proj
    e = torch.ones_like(proj)
    n = torch.numel(e)

    # Mean of the initial Gaussian sample:
    mu = 0.5 * (e - R_T(R(e))) + proj
    # Sample initial image y ~ N(mu, I)
    y = torch.normal(mean=mu, std=1.0).unsqueeze(0).to(device)
    sigma = torch.norm(log_grad(y)) / np.sqrt(n)

    t = 1
    all_ys = []
    # Main loop for denoising
    while sigma > sig_end:
        # update step size
        h = (h_init * t) / (1 + h_init * (t - 1))

        # projected log prior gradient
        # (denoising direction)
        d = log_grad(y).squeeze(0)
        # Project the update
        # (I - A^T A) d       : remove measurable component of prior gradient
        # proj - A^T A y      : data-consistency correction
        d = (d - R_T(R(d)) + proj -
            R_T(R(y.squeeze(0)))).unsqueeze(0)

        # noise magnitude
        sigma = torch.norm(d) / np.sqrt(n)

        # protect against divergence
        div_thld  = 1e2
        iter_thld = 1e3
        if sigma > div_thld or t > iter_thld:
            warnings.warn('Divergence detected, resample with \
                larger step size and tolerance.', RuntimeWarning)

            return linear_inverse(model, render, input, h_init, beta * 2, sig_end * 2,
                                        t_max, stride, seed, msmt_flag, with_grad)

        # inject noise
        gamma = np.sqrt((1 - beta * h) ** 2 - (1 - h) ** 2) * sigma
        noise = torch.randn(size=y.size(), device=device)

        # update image
        y = y + h * d + gamma * noise

        if stride > 0 and (t - 1) % stride == 0:
            print('iter %d, sigma %.2f' % (t, sigma.item()))
            all_ys.append(numpy_image(y))

        t += 1

        # safe guard for iteration limit
        # (typically) use in conjection with grad=True
        # for GPU memory limit
        if t > t_max:
            break

    final = y + log_grad(y)
    all_ys.append(numpy_image(final))

    if with_grad:
        return final.squeeze(0), t, all_ys

    return all_ys

# linear inverse with linear noisy measurements
def noise_inverse(model, render, input, weight,
    h_init=0.01, beta=0.01, sig_end=0.01, stride=10):

    # device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    if torch.backends.mps.is_available() and torch.backends.mps.is_built():
        device = torch.device("mps")
    elif torch.cuda.is_available():
        device = torch.device("cuda")
    else:
        device = torch.device("cpu")
    model = model.eval().to(device)

    # helper function for pytorch image to numpy image
    numpy_image = lambda y: y.detach().cpu() \
        .squeeze(0).permute(1, 2, 0).numpy()

    # the network calculates the noise residual
    def log_grad(y):
        with torch.no_grad():
            return - model(y)

    # measurement matrix calculation
    R = render.measure
    R_T = render.recon

    # init variables
    if input.dim() == 1:
        proj = R_T(weight * input)
    elif input.dim() == 3:
        proj = R_T(weight * R(input))

    e = torch.ones_like(proj)
    n = torch.numel(e)

    mu = 0.5 * (e - R_T(weight * R(e))) + proj
    y = torch.normal(mean=mu, std=1.0).unsqueeze(0).to(device)
    sigma = torch.norm(log_grad(y)) / np.sqrt(n)

    t = 1
    all_ys = []
    while sigma > sig_end:
        # update step size
        h = (h_init * t) / (1 + h_init * (t - 1))

        # projected log prior gradient
        d = log_grad(y).squeeze(0)

        # projected log prior gradient
        # weighted by prior / measurement precision
        d_prior = d - R_T(weight * R(d))
        d_msmt = proj - R_T(weight * R(y.squeeze(0)))

        d = (d_prior + d_msmt).unsqueeze(0)

        # noise magnitude
        sigma = torch.norm(d) / np.sqrt(n)

        # protect against divergence
        div_thld  = 1e2
        iter_thld = 1e3
        if sigma > div_thld or t > iter_thld:
            warnings.warn('Divergence detected, resample with \
                larger step size and tolerance.', RuntimeWarning)

            return noise_inverse(model, render, input, weight,
                        h_init, beta * 2, sig_end * 2, stride)

        # inject noise
        gamma = np.sqrt((1 - beta * h) ** 2 - (1 - h) ** 2) * sigma
        noise = torch.randn(size=y.size(), device=device)

        # update image
        y = y + h * d + gamma * noise

        if stride > 0 and (t - 1) % stride == 0:
            print('iter %d, sigma %.2f' % (t, sigma.item()))
            all_ys.append(numpy_image(y))

        # iteration counter
        t += 1

    final = y + log_grad(y)
    all_ys.append(numpy_image(final))

    return all_ys
