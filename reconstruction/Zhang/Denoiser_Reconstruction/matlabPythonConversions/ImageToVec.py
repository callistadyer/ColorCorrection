import torch


def ImageToVec(x_image):
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

        The input Python image is assumed to have shape (3, H, W), with
        color channel first. To match MATLAB ordering the
        function first reorders the image to shape (H, W, 3) and then
        flattens it.

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
    return x_image.permute(1, 2, 0).flatten()


def VecToImage(x_vec, image_size):
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
        MATLAB ordering convention.

        The input vector is assumed to have been created from a MATLAB-style
        image ordering, where the corresponding image would have shape
        (H, W, 3). This function first reshapes the vector to (H, W, 3),
        and then reorders dimensions to return a Python image of shape
        (3, H, W).

        It exactly inverts ImageToVec, provided that image_size matches
        the size of the original image.
    """
    H, W, C = image_size
    assert C == 3, "Expected image_size = (H, W, 3)"
    return x_vec.reshape(H, W, C).permute(2, 0, 1)


# ---------------------------------
# Small script to show it working
# ---------------------------------

# Make a tiny test image with easy-to-track values
# Shape is (3, H, W) = (3, 2, 3)
x_image = torch.tensor([
    [[1, 2, 3],
     [4, 5, 6]],

    [[10, 20, 30],
     [40, 50, 60]],

    [[100, 200, 300],
     [400, 500, 600]]
], dtype=torch.float32)

print("Original image x_image (shape = 3 x H x W):")
print(x_image)
print("shape:", x_image.shape)

# Convert image to vector
x_vec = ImageToVec(x_image)

print("\nVectorized image x_vec:")
print(x_vec)
print("shape:", x_vec.shape)

# Convert vector back to image
x_image_recovered = VecToImage(x_vec, (2, 3, 3))

print("\nRecovered image x_image_recovered:")
print(x_image_recovered)
print("shape:", x_image_recovered.shape)

# Check exact equality
print("\nDoes VecToImage(ImageToVec(x_image)) recover the original?")
print(torch.equal(x_image, x_image_recovered))

# Also show the intermediate H x W x 3 arrangement
x_hw3 = x_image.permute(1, 2, 0)
print("\nIntermediate image after permute to (H, W, 3):")
print(x_hw3)
print("shape:", x_hw3.shape)