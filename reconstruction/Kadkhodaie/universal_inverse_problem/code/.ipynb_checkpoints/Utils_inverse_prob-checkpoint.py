import numpy as np
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import torch.nn as nn
import os
from skimage.metrics import peak_signal_noise_ratio, structural_similarity
import torch.fft
import gzip
import argparse
from network import BF_CNN

################################################# Helper Functions #################################################
def load_denoiser(architecture,grayscale, training_data, training_noise): 
    if architecture=='BF_CNN': 
        model = load_BF_CNN(grayscale, training_data, training_noise)

        
    return model

def load_BF_CNN(grayscale, training_data, training_noise): 
    '''
    @ grayscale: if True, number of input and output channels are set to 1. Otherwise 3
    @ training_data: models provided in here have been trained on {BSD400, mnist, BSD300}
    @ training_noise: standard deviation of noise during training the denoiser
    '''
    parser = argparse.ArgumentParser(description='BF_CNN_color')
    parser.add_argument('--dir_name', default= '../noise_range_')
    parser.add_argument('--kernel_size', default= 3)
    parser.add_argument('--padding', default= 1)
    parser.add_argument('--num_kernels', default= 64)
    parser.add_argument('--num_layers', default= 20)
    if grayscale is True: 
        parser.add_argument('--num_channels', default= 1)
    else:
        parser.add_argument('--num_channels', default= 3)
    
    args = parser.parse_args('')

    model = BF_CNN(args)
    if torch.cuda.is_available():
        model = model.cuda()
    model_path = os.path.join('denoisers/BF_CNN',training_data,training_noise,'model.pt')
    if torch.cuda.is_available():
        learned_params =torch.load(model_path)

    else:
        learned_params =torch.load(model_path, map_location='cpu' )
    model.load_state_dict(learned_params)
    return model

#################################################
def single_image_loader(data_set_dire_path, image_number):
    
    if 'mnist' in data_set_dire_path.split('/'): 
        f = gzip.open(data_set_dire_path + '/t10k-images-idx3-ubyte.gz','r')
        f.read(16)
        buf = f.read(28 * 28 *10000)
        data = np.frombuffer(buf, dtype=np.uint8).astype(float)/255
        x = torch.tensor(data.reshape( 10000,28, 28).astype('float32'))[image_number:image_number+1]
        
    else: 
        all_names = os.listdir(data_set_dire_path)
        file_name = all_names[image_number]
        im = plt.imread(data_set_dire_path + file_name)
        if len(im.shape) == 3:
            x = torch.tensor(im).permute(2,0,1)
        elif len(im.shape) == 2:
            x = torch.tensor(im.reshape(1, im.shape[0], im.shape[1]))
            
    return x

class test_image: 
    def __init__(self, grayscale,path, image_num):
        super(test_image, self).__init__()
        
        self.grayscale = grayscale
        self.path = path
        self.image_num = image_num
        
        self.im = single_image_loader(self.path,self.image_num)
        if self.im.dtype == torch.uint8: 
            self.im = self.im/255
        if self.im.size()[0] == 3 and grayscale==True: 
            raise Exception('model is trained for grayscale images. Load a grayscale image')
        elif self.im.size()[0] == 1 and grayscale==False: 
            raise Exception('model is trained for color images. Load a color image')
        if torch.cuda.is_available():
            self.im = self.im.cuda()
        
    def show(self):
        if self.grayscale is True: 
            if torch.cuda.is_available():
                plt.imshow(self.im.squeeze(0).cpu(), 'gray', vmin=0, vmax = 1)
            else: 
                plt.imshow(self.im.squeeze(0), 'gray', vmin=0, vmax = 1)                
        else: 
            if torch.cuda.is_available():
                plt.imshow(self.im.permute(1,2,0).cpu(), vmin=0, vmax = 1)
            else: 
                plt.imshow(self.im.permute(1,2,0), vmin=0, vmax = 1)

        plt.title('test image')
        plt.colorbar()
#         plt.axis('off');

    def crop(self, x0,y0,h,w):
        self.cropped_im = self.im[:, x0:x0+h, y0:y0+w]             
        if self.grayscale is True: 
            if torch.cuda.is_available():
                plt.imshow(self.cropped_im.squeeze(0).cpu(), 'gray', vmin=0, vmax = 1)
            else: 
                plt.imshow(self.cropped_im.squeeze(0), 'gray', vmin=0, vmax = 1)
                
        else: 
            if torch.cuda.is_available():            
                plt.imshow(self.cropped_im.permute(1,2,0).cpu(), vmin=0, vmax = 1)
            else: 
                plt.imshow(self.cropped_im.permute(1,2,0), vmin=0, vmax = 1)

        plt.title('cropped test image')
        plt.colorbar()
#         plt.axis('off');        
        return self.cropped_im


#################################################
def rescale_image(im):
    if type(im) == torch.Tensor: 
        im = im.numpy()
    return ((im - im.min()) * (1/(im.max() - im.min()) * 255)).astype('uint8')


def plot_synthesis(intermed_Ys, sample):
    f, axs = plt.subplots(1,len(intermed_Ys), figsize = ( 4*len(intermed_Ys),4))
    axs = axs.ravel()

    #### plot intermediate steps
    for ax in range(len(intermed_Ys)):
        if torch.cuda.is_available():
            intermed_Ys[ax] = intermed_Ys[ax].cpu()
            
        x = intermed_Ys[ax].permute(1,2,0).detach().numpy() 
        if x.shape[2] == 1: # if grayscale
            fig = axs[ax].imshow(x.squeeze(-1), 'gray')
        else: # if color
            fig = axs[ax].imshow(rescale_image(x))
        axs[ax].axis('off')

    #### plot final sample
    if torch.cuda.is_available():
        sample =sample.cpu()
        
    sample = sample.permute(1,2,0).detach().numpy()
    if sample.shape[2] == 1: # if grayscale
        fig = axs[-1].imshow(sample.squeeze(-1),'gray' )
    else: # if color
        fig = axs[-1].imshow(rescale_image(sample))

    axs[-1].axis('off')
    print('value range', np.round(np.min(sample ),2), np.round(np.max(sample),2) )


def plot_sample(x, corrupted, sample):
    if torch.cuda.is_available():
        x = x.cpu()
        corrupted = corrupted.cpu()
        sample = sample.cpu()
        
    x = x.permute(1,2,0)
    corrupted = corrupted.permute(1,2,0)
    sample = sample.detach().permute(1,2,0)
        
    if x.size()!=corrupted.size():    
        h_diff = x.size()[0] - corrupted.size()[0]
        w_diff = x.size()[1] - corrupted.size()[1]
        x = x[0:x.size()[0]-h_diff,0:x.size()[1]-w_diff,: ]
        print('NOTE: psnr and ssim calculated using a cropped original image, because the original image is not divisible by the downsampling scale factor.')
        
    f, axs = plt.subplots(1,3, figsize = (15,5))
    axs = axs.ravel()        
    if x.shape[2] == 1: # if gray scale image
        fig = axs[0].imshow( x.squeeze(-1), 'gray', vmin=0, vmax = 1)
        axs[0].set_title('original')
        
        fig = axs[1].imshow(corrupted.squeeze(-1), 'gray',vmin=0, vmax = 1)
        ssim = np.round(structural_similarity(x.squeeze(-1).numpy(), corrupted.squeeze(-1).numpy()  ) ,3 )
        psnr = np.round(peak_signal_noise_ratio(x.numpy(), corrupted.numpy() ),2)
        axs[1].set_title('corrupted image \n psnr: '+str( psnr) + '\n ssim '+ str(ssim) );  
        
        fig = axs[2].imshow(sample.squeeze(-1),'gray' ,vmin=0, vmax = 1)
        ssim = np.round(structural_similarity(x.squeeze(-1).numpy(), sample.squeeze(-1).numpy()  ) ,3 )
        psnr = np.round(peak_signal_noise_ratio(x.numpy(), sample.numpy() ),2)
        axs[2].set_title('reconstructed \n psnr: '+str( psnr) + '\n ssim '+ str(ssim) );

            
    else: # if color image
        fig = axs[0].imshow( x, vmin=0, vmax = 1)
        axs[0].set_title('original')        
        
        fig = axs[1].imshow( torch.clip(corrupted,0,1), vmin=0, vmax = 1)
        # ssim = np.round(structural_similarity(x.numpy(), corrupted.numpy(), multichannel=True  ) ,3 )
        ssim = np.round(structural_similarity(x.numpy(), corrupted.numpy(), channel_axis=-1, data_range=1.0), 3)

        psnr = np.round(peak_signal_noise_ratio(x.numpy(), corrupted.numpy() ),2)
        axs[1].set_title('corrupted image \n psnr: '+str( psnr) + '\n ssim '+ str(ssim) );  
        
        fig = axs[2].imshow(torch.clip(sample, 0,1),vmin=0, vmax = 1)
        # ssim = np.round(structural_similarity(x.numpy(), sample.numpy() , multichannel=True) ,3)
        ssim = np.round(structural_similarity(x.numpy(), sample.numpy(), channel_axis=-1, data_range=1.0), 3)
        psnr = np.round(peak_signal_noise_ratio(x.numpy(), sample.numpy() ),2)   
        axs[2].set_title('reconstructed \n psnr: '+str( psnr) + '\n ssim '+ str(ssim) );
            
            
    for i in range(3): 
        axs[i].axis('off')
    
    


def plot_all_samples(sample, intermed_Ys):
    n_rows = int(np.ceil(len(intermed_Ys)/4))

    f, axs = plt.subplots(n_rows,4, figsize = ( 4*4, n_rows*4))
    axs = axs.ravel()

    #### plot intermediate steps
    for ax in range(len(intermed_Ys)):
        if torch.cuda.is_available():
            intermed_Ys[ax] = intermed_Ys[ax].cpu()
            
        x = intermed_Ys[ax].detach().permute(1,2,0).numpy()
        if x.shape[2] == 1:
            fig = axs[ax].imshow(x.squeeze(-1), 'gray')
        else:
            fig = axs[ax].imshow(rescale_image(x))
        axs[ax].axis('off')
    
    #### plot final sample
    if torch.cuda.is_available():
        sample =sample.cpu()
        
    sample = sample.detach().permute(1,2,0).numpy()
    if sample.shape[2] == 1:
        fig = axs[-1].imshow(sample.squeeze(-1),'gray' )
    else:
        fig = axs[-1].imshow(rescale_image(sample))
    axs[-1].axis('off')
    plt.colorbar(fig, ax=axs[-1], fraction=.05)


    for ax in range(len(intermed_Ys),n_rows*4 ):
        axs[ax].axis('off')


def plot_corrupted_im(x_c): 
    try:

        if torch.cuda.is_available():
            plt.imshow(x_c.squeeze(0).cpu(), 'gray', vmin=0, vmax = 1)
        else: 
            plt.imshow(x_c.squeeze(0), 'gray', vmin=0, vmax = 1)
    except TypeError: 
        if torch.cuda.is_available():
            plt.imshow(x_c.permute(1,2,0).cpu(), vmin=0, vmax = 1)
        else: 
            plt.imshow(x_c.permute(1,2,0) , vmin=0, vmax = 1)

    plt.colorbar()    
    

def print_dim(measurment_dim, image_dim):
    print('*** Retained {} / {} ({}%) of dimensions'.format(int(measurment_dim), image_dim
                                                   , np.round(int(measurment_dim)/int(image_dim)*100,
                                                              decimals=3) ))    
    
###################################### Inverse problems Tasks ##################################
#############################################################################################
class synthesis:
    def __init__(self):
        super(synthesis, self).__init__()

    def M_T(self, x):
        return torch.zeros_like(x)

    def M(self, x):
        return torch.zeros_like(x)

class inpainting:
    '''
    makes a blanked area in the center
    @x_size : image size, tuple of (n_ch, im_d1,im_d2)
    @x0,y0: center of the blanked area
    @w: width of the blanked area
    @h: height of the blanked area
    '''
    def __init__(self, x_size,x0,y0,h, w):
        super(inpainting, self).__init__()

        n_ch , im_d1, im_d2 = x_size
        self.mask = torch.ones(x_size)
        if torch.cuda.is_available():
            self.mask = self.mask.cuda()
        c1, c2 = int(x0), int(y0)
        h , w= int(h/2), int(w/2)
        self.mask[0:n_ch, c1-h : c1+h , c2-w:c2+w] = 0

    def M_T(self, x):
        return x*self.mask

    def M(self, x):
        return x*self.mask

    
class rand_pixels:
    '''
    @x_size : tuple of (n_ch, im_d1,im_d2)
    @p: fraction of dimensions kept in (0,1)
    '''
    def __init__(self, x_size, p):
        super(rand_pixels, self).__init__()

        self.mask = np.zeros(x_size).flatten()
        self.mask[0:int(p*np.prod(x_size))] = 1
        self.mask = torch.tensor(np.random.choice(self.mask, size = x_size , replace = False).astype('float32').reshape(x_size))
        if torch.cuda.is_available():
            self.mask = self.mask.cuda()        
        
    def M_T(self, x):                                                                                                       
        return x*self.mask

    def M(self, x):
        return x*self.mask

    
class super_resolution:
    '''
    block averaging for super resolution.
    creates a low rank matrix (thin and tall) for down sampling
    @s: downsampling factor, int
    @x_size: tuple of three int  (n_ch, im_d1, im_d2)
    '''

    def __init__(self, x_size, s):
        super(super_resolution, self).__init__()

#         if x_size[1]%2 !=0 or x_size[2]%2 != 0 :
#             raise Exception("image dimensions need to be even")

        self.down_sampling_kernel = torch.ones(x_size[0],1,s,s)
        self.down_sampling_kernel = self.down_sampling_kernel/np.linalg.norm(self.down_sampling_kernel[0,0])
        if torch.cuda.is_available():
            self.down_sampling_kernel = self.down_sampling_kernel.cuda()
        self.x_size = x_size
        self.s = s

    def M_T(self, x):
        down_im = torch.nn.functional.conv2d(x.unsqueeze(0), self.down_sampling_kernel, stride= self.s, groups = self.x_size[0])
        return down_im[0]

    def M(self, x):
        rec_im = torch.nn.functional.conv_transpose2d(x.unsqueeze(0), self.down_sampling_kernel, stride= self.s, groups = self.x_size[0])

        return rec_im[0]


    
class random_basis:
    '''
    @x_size : tuple of (im_d1,im_d2)
    @p: fraction of dimensions kept in (0,1)
    '''
    def __init__(self, x_size, p):
        super(random_basis, self).__init__()
        n_ch , im_d1, im_d2 = x_size
        self.x_size = x_size
        self.U, _ = torch.qr(torch.randn(int(np.prod(x_size)),int(np.prod(x_size)*p) ))
        if torch.cuda.is_available():
            self.U = self.U.cuda()

    def M_T(self, x):
        # gets 2d or 3d image and returns flatten partial measurement(1d)
        # Flatten the image into a single vector
        # x_flat contains all pixel values stacked into a length-N vector
        x_flat = x.flatten()   # shape: (N,)
    
        # Project the image onto each of the K basis vectors
        # This computes K dot products: x_c[i] = <u_i, x>
        x_c = torch.matmul(self.U.T, x_flat)   # shape: (K,)

        # return torch.matmul(self.U.T,x.flatten())
        return x_c
    
    def M(self, x):
        # gets flatten partial measurement (1d), and returns 2d or 3d reconstruction
    
        # Step 1: Linearly combine the K basis vectors using the K coefficients
        # This computes: sum_i x_c[i] * u_i
        # Result is a flattened image vector of length N = C*H*W
        x_proj_flat = torch.matmul(self.U, x_c)   # shape: (N,)
    
        # Reshape the flat vector back into image form
        C, H, W = self.x_size
        x_proj = x_proj_flat.reshape(C, H, W)
        
        # return torch.matmul(self.U,x).reshape(self.x_size[0], self.x_size[1], self.x_size[2])
        return x_proj
    

#### important: when using fftn from torch the reconstruction is more lossy than when fft2 from numpy
#### the difference between reconstruction and clean image in pytorch is of order of e-8, but in numpy is e-16

class spectral_super_resolution:
    '''
    creates a mask for dropping high frequency coefficients
    @im_d: dimension of the input image is (im_d, im_d)
    @p: portion of coefficients to keep
    '''
    def __init__(self, x_size, p):
        super(spectral_super_resolution, self).__init__()

        self.x_size = x_size
        
        xf = int(round(x_size[1]*np.sqrt(p) )/2)
        yf = int(round(x_size[1]*x_size[2]*p/xf )/4)
                
        mask = torch.ones((x_size[1],x_size[2]))

        mask[xf:x_size[1]-xf,:]=0
        mask[:, yf:x_size[2]-yf]=0        
        self.mask = mask
        if torch.cuda.is_available():
            self.mask = self.mask.cuda()

    def M_T(self, x):
        # returns fft of each of the three color channels independently
        return self.mask*torch.fft.fftn(x, norm= 'ortho', dim = (1,2),s = (self.x_size[1], self.x_size[2]) )

    def M(self, x):
        return torch.real(torch.fft.ifftn(x, norm= 'ortho',  dim = (1,2), s = (self.x_size[1], self.x_size[2]) ))


class dichromat:
    """
    Dichromat inverse problem in LMS space.
    The measurement is produced by projecting LMS contrast onto a dichromat plane.
    """
    def __init__(self, x_size, dichromat_type, grayLMS, wls, T_cones, constraint_wl=None):
        # Store relevent info 
        self.x_size         = x_size  # expected LMS image size (C, H, W)
        self.dichromat_type = dichromat_type # dichromat type (Deuteranopia, Protanopia, Tritanopia)
        self.grayLMS        = torch.as_tensor(grayLMS, dtype=torch.float32).reshape(3, 1, 1)
        self.wls            = torch.as_tensor(wls, dtype=torch.float32)
        self.T_cones        = torch.as_tensor(T_cones, dtype=torch.float32) # cone fundamentals (3 x N)

        # Choose default constraint wavelength if not provided
        if constraint_wl is None:
            if dichromat_type == 'Deuteranopia':
                constraint_wl = 575  # or 475 (blue)
            elif dichromat_type == 'Protanopia':
                constraint_wl = 575  # or 475 (blue)
            elif dichromat_type == 'Tritanopia':
                constraint_wl = 660  # or 485 (blue-green)
            else:
                raise ValueError('dichromat_type must be one of: Deuteranopia, Protanopia, Tritanopia')

        # Set missing/available cone indices for LMS ordering (L=0, M=1, S=2)
        if dichromat_type == 'Deuteranopia':
            self.missing_idx = 1
            self.available_idx = [0, 2]
        elif dichromat_type == 'Protanopia':
            self.missing_idx = 0
            self.available_idx = [1, 2]
        elif dichromat_type == 'Tritanopia':
            self.missing_idx = 2
            self.available_idx = [0, 1]

        # Find the closest wavelength index for the monochromatic constraint
        wl_idx = torch.argmin(torch.abs(self.wls - float(constraint_wl))).item()
        # Extract the cone response vector at that wavelength.
        constraint2LMS = self.T_cones[:, wl_idx].reshape(3)

        # Achromatic constraint direction in LMS contrast
        constraint1LMScontrast = torch.tensor([1.0, 1.0, 1.0], dtype=torch.float32)
        # Define monochromatic constraint direction in LMS contrast
        constraint2LMScontrast = (constraint2LMS - self.grayLMS.reshape(3)) / self.grayLMS.reshape(3)

        # Stack the two constraint directions column-wise to form a 3x2 matrix
        # Each column is a direction in LMS contrast space
        # The column span defines the 2D plane of contrasts "accessible" to a dichromat
        constraint_matrix = torch.stack([constraint1LMScontrast, constraint2LMScontrast], dim=1)

        A = constraint_matrix[self.available_idx, :]              # available-cone rows (2x2)
        B = constraint_matrix[self.missing_idx, :].reshape(1, 2)  # missing-cone row (1x2)
        # Solve for the missing-cone linear combination coefficients
        # Compute coefficients that reconstruct the missing cone contrast
        # from the two available cone contrasts...
        # Any dichromat color must lie in the 2D plane spanned by the
        # achromatic and monochromatic constraint directions
        # Solving row * A = B finds the unique linear rule that makes the
        # missing cone consistent with those same plane coordinates
        
        # A = plane -> visible cones
        # B = plane -> missing cone
        # row = visible cones -> missing cone

        
        row = torch.linalg.solve(A.T, B.T).T.reshape(2)
        # Build the trichromat-to-dichromat contrast projection matrix
        self.M_triToDi                      = torch.eye(3, dtype=torch.float32)
        # Zero out the row corresponding to the missing cone
        # This clears any direct contribution from the original missing-cone contrast
        # so its value will be fully recomputed from the remaining cones
        self.M_triToDi[self.missing_idx, :] = 0.0  # Zero the missing-cone row
        # Fill missing-cone row with the linear combination of available cones
        self.M_triToDi[self.missing_idx, self.available_idx] = row

    # Trichromat -> Dichromat (always returns full 3-channel projected LMS image)
    def M_T(self, x):
        # x is LMS image with shape (3, H, W).
        M           = self.M_triToDi.to(device=x.device, dtype=x.dtype)
        gray        = self.grayLMS.to(device=x.device, dtype=x.dtype)
        x_contrast  = (x - gray) / gray
        # Apply the dichromat projection in contrast space
        di_contrast = torch.matmul(M, x_contrast.reshape(3, -1)).reshape_as(x_contrast)
        di_lms      = (di_contrast * gray) + gray # Convert contrast back to LMS excitations
        return di_lms

    # Dichromat -> Trichromat (This just returns itself...)
    def M(self, x):
        # The goal here is to just make sure you have a matrix that is 3xnPixels, which the dichromat version will already be
        return x




