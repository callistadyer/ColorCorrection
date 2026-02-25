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
    Implicit (2N x 3N) block-diagonal measurement operator using a per-pixel matrix A (2x3)
    Equivalent to R_big = I_N ⊗ A, without explicitly building R_big
    """

    def __init__(self, A_2x3, im_size, device):
        """
        A_2x3: torch.Tensor of shape (2,3)  (your ortho_mtx_small)
        im_size: (3, H, W)
        """
        self.im_size = im_size
        self.device  = device
        self.A       = A_2x3.to(device).float()     # (2,3)
        self.AT      = self.A.T                     # (3,2)

    def to(self, device):
        self.device = device
        self.A = self.A.to(device)
        self.AT = self.A.T
        return self
    
    def measure(self, x):
        """
        x: (3,H,W)
        returns msmt: (2N,) flattened
        """
        if self.A.device != x.device:
            self.to(x.device)
        H, W = x.shape[1], x.shape[2]
        # (3,H,W) -> (N,3)
        xN3 = x.permute(1, 2, 0).reshape(-1, 3)
        # (N,3) -> (N,2)
        yN2 = xN3 @ self.A.T
        # return (2N,)
        return yN2.reshape(-1)

    def recon(self, msmt):
        """
        msmt: (2N,)
        returns: (3,H,W)
        """
        if self.A.device != msmt.device:
            self.to(msmt.device)
        H, W = self.im_size[1], self.im_size[2]
        # (2N,) -> (N,2)
        yN2 = msmt.reshape(-1, 2)
        # (N,2) -> (N,3)
        xN3 = yN2 @ self.AT.T   # (N,2) @ (2,3) = (N,3)
        # (N,3) -> (3,H,W)
        x = xN3.reshape(H, W, 3).permute(2, 0, 1)
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
    model = model.eval().to(device)
    
    if hasattr(render, "to"):
        render = render.to(device)
    if torch.is_tensor(input):
        input = input.to(device)

    # helper function for pytorch image to numpy image
    numpy_image = lambda y: y.detach().cpu() \
                .squeeze(0).permute(1, 2, 0).numpy()

    # the network calculates the noise residual
    if with_grad:
        log_grad = lambda y: - model(y)
    else:
        def log_grad(y):
            with torch.no_grad():
                return - model(y)

    # measurement matrix calculation
    R = render.measure
    R_T = render.recon

    # init variables
    if msmt_flag or input.dim() == 1:
        proj = R_T(input)
    elif input.dim() == 3:
        proj = R_T(R(input))

    e = torch.ones_like(proj)
    n = torch.numel(e)

    mu = 0.5 * (e - R_T(R(e))) + proj
    y = torch.normal(mean=mu, std=1.0).unsqueeze(0).to(device)
    sigma = torch.norm(log_grad(y)) / np.sqrt(n)

    t = 1
    all_ys = []
    while sigma > sig_end:
        # update step size
        h = (h_init * t) / (1 + h_init * (t - 1))

        # projected log prior gradient
        d = log_grad(y).squeeze(0)
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
