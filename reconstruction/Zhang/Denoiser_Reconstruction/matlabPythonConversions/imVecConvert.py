import torch
import numpy as np

def ImToVec(x_image):
    """
    Converts an image into one long vector using MATLAB ordering

    Syntax:
        x_vec = ImageToVec(x_image)

    Inputs:
        x_image:    Input image of size (3, H, W)

    Outputs:
        x_vec:      Output vector of size (3*H*W,)

    Description:
        This function converts an image into one long vector using the
        MATLAB ordering convention

        The input image is assumed to have shape (3, H, W). The function
        first reorders it to (H, W, 3), and then vectorizes using
        column-major ordering so that the result matches MATLAB xImage(:).

        Thus this function should output a vector with ordering that matches 
        the MATLAB routine ImageToVec.

    History:
        03/10/2026  cmd    Wrote it.

    Examples:
    #{
    import torch

    x_image = torch.rand(3, 4, 5)
    x_vec = ImageToVec(x_image)
    print(x_vec.shape)
    #}
    """
    # Reorder image from (3, H, W) to (H, W, 3), then convert to NumPy
    x_hwc = x_image.permute(1, 2, 0).cpu().numpy()
    # reshape with (-1,) collapses all dimensions together
    # so shape changes from (H, W, 3) to (H*W*3,)
    # order='F' changes the order entries are read out, but not the final shape 
    x_vec = np.reshape(x_hwc, (-1,), order='F')
    # Convert back to torch and move to the original device
    return torch.from_numpy(x_vec).to(x_image.device)


def VecToIm(x_vec, image_size):
    """
    Converts a vector back into an image using MATLAB ordering

    Syntax:
        x_image = VecToImage(x_vec, image_size)

    Inputs:
        x_vec:        Input vector of size (3*H*W,)
        image_size:   Desired output image size (H, W, 3)

    Outputs:
        x_image:      Output image of size (3, H, W)

    Description:
        This function reshapes a vector back into an image using the
        MATLAB ordering convention

        This function inverts ImageToVec. It reshapes the vector into an
        image of shape (H, W, 3) using MATLAB column-major ordering, and
        then reorders dimensions to return a Python image of shape
        (3, H, W).
    """
    H, W, C = image_size
    assert C == 3, "Expected image_size = (H, W, 3)"

    # x_vec starts with shape (H*W*3,)
    # reshape(..., (H, W, C), order='F') changes shape from (H*W*3,) to (H, W, 3)
    # order='F' makes the entries go back into the array using MATLAB column-major ordering
    x_hwc = np.reshape(x_vec.cpu().numpy(), (H, W, C), order='F')

    # torch.from_numpy keeps the same shape, so it stays (H, W, 3)
    # permute(2, 0, 1) changes shape from (H, W, 3) to (3, H, W)        
    x_image = torch.from_numpy(x_hwc).permute(2, 0, 1)

    return x_image.to(x_vec.device)



def CalFormatToVec(cal_format):
    """
    Converts a CalFormat matrix into one long vector using MATLAB ordering

    Syntax:
        x_vec = CalFormatToVec(cal_format)

    Inputs:
        cal_format:   Input CalFormat matrix of size (3, N)

    Outputs:
        x_vec:        Output vector of size (3*N,)

    Description:
        This function converts a CalFormat matrix into one long vector using
        MATLAB ordering.

        The input Python CalFormat matrix is assumed to have shape (3, N).
        To match MATLAB vectorization of a 3 x N matrix, this function
        reshapes using column-major ordering.

        Thus if the columns of cal_format correspond to pixels, the output
        vector is ordered as
            [R1, G1, B1, R2, G2, B2, ...].

    History:
        03/10/2026  cmd    Wrote it.
    """

    # cal_format starts with shape (3, N)
    cal_np = cal_format.cpu().numpy()

    # reshape with (-1,) collapses all dimensions together
    # so shape changes from (3, N) to (3*N,)
    # order='F' makes NumPy read entries in MATLAB / column-major order
    x_vec = np.reshape(cal_np, (-1,), order='F')

    # torch.from_numpy keeps the same shape, so shape stays (3*N,)
    # .to(cal_format.device) changes device only, not shape
    return torch.from_numpy(x_vec).to(cal_format.device)


def VecToCalFormat(x_vec, cal_format_size):
    """
    Converts a vector back into a CalFormat matrix using MATLAB ordering

    Syntax:
        cal_format = VecToCalFormat(x_vec, cal_format_size)

    Inputs:
        x_vec:              Input vector of size (3*N,)
        cal_format_size:    Desired output CalFormat size (3, N)

    Outputs:
        cal_format:         Output CalFormat matrix of size (3, N)

    Description:
        This function reshapes a CalFormat vector back into a CalFormat
        matrix using MATLAB ordering.

        The input vector is assumed to have been created from a MATLAB-style
        3 x N CalFormat matrix, where each column contains the 3 values for
        one pixel. Thus the vector is assumed to be ordered as

            [R1, G1, B1, R2, G2, B2, ..., RN, GN, BN]

        where pixel 1, pixel 2, etc. refer to the columns of the original
        CalFormat matrix.

        This function reshapes that vector into shape (3, N) using MATLAB
        column-major ordering.

        It inverts CalFormatToVec, provided that cal_format_size matches
        the size of the original CalFormat matrix.

        !!!!!!!Note that this function does not expect an image vector produced by
        ImToVec. Image vectors and CalFormat vectors are different
        conventions, even though both may have size (3*N,).

    History:
        03/10/2026  cmd    Wrote it.
    """

    # cal_format_size should be (3, N)
    C, N = cal_format_size

    # checks that the first dimension is 3
    # assert C == 3, "Expected cal_format_size = (3, N)"

    # x_vec starts with shape (3*N,)
    # reshape(..., (3, N), order='F') changes shape from (3*N,) to (3, N)
    # order='F' puts entries back using MATLAB / column-major ordering
    cal_np = np.reshape(x_vec.cpu().numpy(), (C, N), order='F')

    return torch.from_numpy(cal_np).to(x_vec.device)

def ImToCalFormat(x_image):
    """
    Converts an image into CalFormat using MATLAB pixel ordering

    Syntax:
        cal_format = ImToCalFormat(x_image)

    Inputs:
        x_image:      Input image of size (3, H, W)

    Outputs:
        cal_format:   Output CalFormat matrix of size (3, N), where N = H*W

    Description:
        This function converts a Python image into CalFormat using MATLAB
        pixel ordering.

        The input image is assumed to have shape (3, H, W). The output
        CalFormat matrix has shape (3, N), where each column contains the
        3 values for one pixel.

        Pixel order matches MATLAB image ordering, meaning pixels are taken
        column-by-column across the image.

    History:
        03/10/2026  cmd    Wrote it.
    """

    # x_image starts with shape (3, H, W)

    # permute changes shape from (3, H, W) to (W, H, 3)
    # this swaps the spatial dimensions so that when reshape is applied
    # the pixels are read out in MATLAB column-major image order
    x_whc = x_image.permute(2, 1, 0)

    # reshape changes shape from (W, H, 3) to (N, 3), where N = H*W
    # now each row contains the RGB values for one pixel
    xN3 = x_whc.reshape(-1, 3)

    # transpose changes shape from (N, 3) to (3, N)
    # now each column contains one pixel's RGB values
    cal_format = xN3.T

    return cal_format

# 1) What CalFormatToIm expects
#
# CalFormatToIm expects a true CalFormat matrix of shape (3, N), where
# each COLUMN is one pixel, and pixels are ordered in MATLAB pixel order.
#
# EXPECTED:
#
#   cal_format =
#   [ R1   R2   R3   ...  RN
#     G1   G2   G3   ...  GN
#     B1   B2   B3   ...  BN ]
#
# So CalFormatToIm assumes:
#   - rows = channels
#   - columns = pixels


def CalFormatToIm(cal_format, image_size):
    """
    Converts a CalFormat matrix back into an image using MATLAB pixel ordering

    Syntax:
        x_image = CalFormatToIm(cal_format, image_size)

    Inputs:
        cal_format:   Input CalFormat matrix of size (3, N)
        image_size:   Desired output image size (H, W, 3)

    Outputs:
        x_image:      Output image of size (3, H, W)

    Description:
        This function converts a CalFormat matrix back into a Python image
        using MATLAB pixel ordering.

        The input CalFormat matrix is assumed to have shape (3, N), where
        each column contains the 3 values for one pixel, ordered according
        to MATLAB image pixel ordering.

        It exactly inverts ImToCalFormat, provided that image_size matches
        the original image size.

    History:
        03/10/2026  cmd    Wrote it.
    """

    H, W, C = image_size
    assert C == 3, "Expected image_size = (H, W, 3)"

    # cal_format starts with shape (3, N)

    # transpose changes shape from (3, N) to (N, 3)
    # now each row contains one pixel's RGB values
    xN3 = cal_format.T

    # reshape changes shape from (N, 3) to (W, H, 3)
    # this reverses the reshape used in ImToCalFormat
    x_whc = xN3.reshape(W, H, 3)

    # permute changes shape from (W, H, 3) to (3, H, W)
    # this returns the image to standard Python channel-first format
    x_image = x_whc.permute(2, 1, 0)

    return x_image