def ReadCalFormatFromMatlab(rgb_cal, image_size, dtype=None, device=None):
    """
    Convert a MATLAB CalFormat RGB matrix into a Python image tensor.

    Syntax
    ------
    x_image = ReadCalFormatFromMatlab(rgb_cal, image_size, dtype=None, device=None)

    Inputs
    ------
    rgb_cal : numpy.ndarray or torch.Tensor
        CalFormat matrix of shape (3, N), where N = H*W.
        Expected structure:

            [ R1   R2   R3   ...   RN
              G1   G2   G3   ...   GN
              B1   B2   B3   ...   BN ]

        with pixels ordered in MATLAB column-major image order.

    image_size : tuple
        Either (H, W) or (H, W, 3)

    dtype : torch.dtype, optional
        Output dtype. If None, defaults to torch.float32.

    device : torch.device or str, optional
        Output device. If None, leaves tensor on CPU.

    Output
    ------
    x_image : torch.Tensor
        Image tensor of shape (3, H, W)

    Description
    -----------
    MATLAB CalFormat stores each color channel as one row, flattened in
    column-major order. To reconstruct the image correctly in Python, each
    channel must be reshaped separately using MATLAB/Fortran ordering, then
    stacked into an (H, W, 3) image, and finally permuted to channel-first
    format (3, H, W).

    History
    -------
    03/11/2026  cmd   Wrote it.
    """
    import numpy as np
    import torch

    if len(image_size) == 2:
        H, W = image_size
    elif len(image_size) == 3:
        H, W, C = image_size
        if C != 3:
            raise ValueError("Expected image_size = (H, W, 3)")
    else:
        raise ValueError("image_size must be (H, W) or (H, W, 3)")

    if isinstance(rgb_cal, torch.Tensor):
        rgb_cal_np = rgb_cal.detach().cpu().numpy()
    else:
        rgb_cal_np = np.asarray(rgb_cal)

    if rgb_cal_np.shape != (3, H * W):
        raise ValueError(f"Expected rgb_cal shape (3, {H*W}), got {rgb_cal_np.shape}")

    img_hw3 = np.stack([
        rgb_cal_np[0, :].reshape((H, W), order='F'),
        rgb_cal_np[1, :].reshape((H, W), order='F'),
        rgb_cal_np[2, :].reshape((H, W), order='F')
    ], axis=2)  # (H, W, 3)

    if dtype is None:
        dtype = torch.float32

    x_image = torch.from_numpy(img_hw3).permute(2, 0, 1).to(dtype=dtype)

    if device is not None:
        x_image = x_image.to(device)

    return x_image