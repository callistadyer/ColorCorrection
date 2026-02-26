# Denoiser Reconstruction Example

This example demonstrates image reconstruction from compressed measurements using a learned denoiser prior.

## Setup

### 1. Create a Python virtual environment

```bash
cd example_code
python3 -m venv venv
```

### 2. Activate the virtual environment

```bash
source venv/bin/activate
```

### 3. Install dependencies

```bash
pip install -r requirements.txt
```

Note: The `cupy-cuda12x` package requires CUDA 12.x. If you have a different CUDA version, install the appropriate CuPy package:
- CUDA 11.x: `pip install cupy-cuda11x`
- CUDA 10.x: `pip install cupy-cuda10x`

### 4. Register Jupyter kernel (optional, for notebook use)

```bash
python -m ipykernel install --user --name=denoiser_example --display-name="Denoiser Example"
```

## Running the Code

### Jupyter Notebook

1. Open `recon_visualize.ipynb`

2. Select the "Denoiser Example" kernel (or your venv kernel)

3. Run the cells sequentially

## Notebook Contents

The notebook demonstrates:

1. **Dataset Loading**: Load and visualize ILSVRC test images

2. **Denoising Demo**: Add Gaussian noise to images and denoise using the trained CNN denoiser

3. **Reconstruction Demo**:
   - Create a random orthogonal measurement matrix (5% of original dimension)
   - Simulate compressed measurements
   - Reconstruct images using iterative denoising

4. **Multi-image Reconstruction**: Visualize reconstruction results on multiple test images

## File Structure

```
example_code/
├── README.md
├── requirements.txt
├── main.py                 # Configuration module
├── test_run.py             # Test script
├── recon_visualize.ipynb   # Main notebook
├── assets/
│   ├── conv3_ln.pt         # Trained denoiser weights
│   └── gamma.mat           # Gamma correction table
├── utils/
│   ├── dataset.py          # Dataset loading utilities
│   ├── helper.py           # Helper functions
│   └── dataset/islvrc/     # ILSVRC dataset
│       ├── train/
│       └── test/
├── models/
│   └── denoiser.py         # CNN denoiser model
└── inverse/
    ├── solver.py           # Linear inverse solver
    └── sampler.py          # Prior sampler
```

