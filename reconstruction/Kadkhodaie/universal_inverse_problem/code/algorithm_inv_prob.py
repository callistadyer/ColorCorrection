import numpy as np
import torch
import time
import os



### Takes a tensor of size (n_ch, im_d1, im_d2)
### and returns a tensor of size (n_ch, im_d1, im_d2)
def univ_inv_sol(model, x_c ,task ,sig_0=1, sig_L=.01, h0=.01 , beta=.01 , freq=5):
    '''
    @x_c:  M^T.x)
    @task: the specific linear inverse problem
    @sig_0: initial sigma (largest)
    @sig_L: final sigma (smallest)
    @h0: 1st step size
    @beta:controls added noise in each iteration (0,1]. if 1, no noise is added. As it decreases more noise added.
    '''

    # M_T is the function that simulates how the image was degraded or sampled to produce the observed data.
    M_T = task.M_T #low rank measurement matrix - in function form
    M = task.M #inverse of M_T
    
    n_ch, im_d1,im_d2 = M(x_c).size()
    N = n_ch* im_d1*im_d2
    intermed_Ys=[]

    # initialize y

    # Create an image-shaped thing filled with 1s (white image)
    e =  torch.ones_like(M(x_c), requires_grad= False )

    # LETS BREAK THIS UP INTO MORE READABLE LINES BELOW
    # y = torch.normal((e - M(M_T(e)))*.5 + M(x_c), sig_0)
    # y = y.unsqueeze(0)
    # y.requires_grad = False

    # Where pixels survive measurement (1 = observed, 0 = missing)
    observed_mask = M_T(e)
    
    # Map mask back to image space
    observed_mask_image = M(observed_mask)
    
    # flip it: 1 = missing, 0 = observed
    missing_mask = e - observed_mask_image
        
    # Build a reasonable starting image:
    # use the observed pixels from the corrupted image
    # fill missing pixels with gray (0.5 = mid-gray)
    # mean_init = M(x_c) + 0.5 * missing_mask
    
    measured_component = M(x_c)

    # Neutral initialization for the unmeasured degrees of freedom
    gray_level = 0.5
    unmeasured_component = missing_mask
    
    # Combine measured data with neutral fill for unconstrained components
    #“Use the correct values wherever the measurements tell us what to do, 
    # and use a neutral placeholder everywhere else.”

    # gray_level * unmeasured_component turns missing part to gray
    # 
    mean_init = measured_component + gray_level * unmeasured_component
    
    # Add Gaussian noise to the initial image estimate
    # sig_0 controls how noisy the starting point is
    y = torch.normal(mean_init, sig_0)
    
    # Add a batch dimension (model expects shape [1, C, H, W])
    y = y.unsqueeze(0)
    
    # We are not optimizing y with backprop. updates are done manually
    y.requires_grad = False

    if freq > 0:
        intermed_Ys.append(y.squeeze(0))


    if torch.cuda.is_available():
        y = y.cuda()

    f_y = model(y)


    sigma = torch.norm(f_y)/np.sqrt(N)


    t=1
    start_time_total = time.time()
    # Keep iterating until the sigma is small enough
    while sigma > sig_L:

        # Compute the step size h for this iteration.
        # Early iterations: larger steps (coarse exploration)
        # Later iterations: smaller steps (fine refinement)
        h = h0*t/(1+ (h0*(t-1)) )

        # Apply the denoiser (image prior)
        # Run the denoiser on the current image estimate y
        # torch.no_grad(): this is not backprop, the denoiser is a fixed operator
        with torch.no_grad():
            f_y = model(y)

        
        ####### Compute the update direction d #######
        # This combines:
        #   (1) a denoising / prior-driven term
        #   (2) a data-consistency correction term
        #
        # f_y:
        #   what the denoiser thinks is noise / residual in y
        #
        # M(M_T(f_y[0])):
        #   projection of the denoiser output onto the measured subspace
        #
        # f_y - M(M_T(f_y[0])):
        #   keep only the denoiser's suggestion in the unmeasured directions
        #
        # M(M_T(y[0])) - M(x_c):
        #   how much the current estimate violates the measurements
        #   (forces y back toward measurement consistency)
        
        d = f_y - M(M_T(f_y[0])) + ( M(M_T(y[0]))  - M(x_c) )


        # Measure how big the update is
        # sigma is the RMS magnitude of the update direction.
        # It acts like a "noise level" indicator
        sigma = torch.norm(d)/np.sqrt(N)

        # gamma controls how much random noise is injected
        #
        # Early iterations:
        #   sigma is large -> gamma is large -> more randomness (exploration)
        #
        # Late iterations:
        #   sigma is small -> gamma is small -> mostly deterministic refinement
        #
        # beta controls how much noise is allowed:
        #   beta = 1   -> no extra randomness
        #   beta < 1   -> more noise
        gamma = sigma*np.sqrt(((1 - (beta*h))**2 - (1-h)**2 ))

        # Sample random noise in image space
        noise = torch.randn(n_ch, im_d1,im_d2)

        if torch.cuda.is_available():
            noise = noise.cuda()

        # Update the image estimate
        #   y             : current estimate
        #   - h * d       : step size times update direction
        #   + gamma*noise : some randomness
        y = y -  h*d + gamma*noise

        if freq > 0 and t%freq== 0:
            print('-----------------------------', t)
            print('sigma ' , sigma.item() )
            intermed_Ys.append(y.squeeze(0))


        t +=1


    print("-------- total number of iterations, " , t )
    print("-------- average time per iteration (s), " , np.round((time.time() - start_time_total)/(t-1)  ,4) )

    denoised_y = y - model(y)


    return denoised_y.squeeze(0), intermed_Ys


