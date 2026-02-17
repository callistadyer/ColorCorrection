"""Test script to verify the example code works."""
import os

# Change to example_code directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

print("Testing example_code setup...")
print("=" * 50)

# Test 1: Import modules
print("\n1. Testing imports...")
import main
import torch
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

from utils.dataset import DataSet, gamma_correct
from models.denoiser import Denoiser
from inverse.solver import RenderMatrix, linear_inverse
print("   All imports successful!")

# Test 2: Load dataset
print("\n2. Loading dataset...")
main.args.data_path = 'islvrc'
main.args.test_scale = [0.385]
test_set = DataSet.load_dataset(main.args, test_mode=True).test_set()[:40]
print(f"   Loaded {len(test_set)} test images, shape: {test_set[0].shape}")

# Test 3: Load denoiser
print("\n3. Loading denoiser model...")
main.args.model_path = './assets/conv3_ln.pt'
model = Denoiser(main.args)
model.load_state_dict(torch.load(main.args.model_path, weights_only=True))
model = model.eval()
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"   Model loaded, using device: {device}")

# Test 4: Create random sampling matrix
print("\n4. Creating random sampling matrix (5% of original dimension)...")
im_dim = 128 * 128 * 3  # 49152
msmt_dim = im_dim // 20  # 5%

np.random.seed(42)
random_mtx = np.random.randn(msmt_dim, im_dim).astype(np.float32)

# Orthogonalize using SVD
U, S, Vt = np.linalg.svd(random_mtx, full_matrices=False)
ortho_mtx = (U @ Vt).astype(np.float32)
print(f"   Matrix shape: {ortho_mtx.shape}")

# Test 5: Run reconstruction
print("\n5. Running reconstruction...")
test_set_torch = torch.tensor(test_set).permute([0, 3, 1, 2]).to(device)
im_size = test_set_torch.size()[1:]

render = RenderMatrix(torch.tensor(ortho_mtx), im_size, device)

test_idx = 5
test_img = test_set_torch[test_idx]
print(f"   Test image shape: {test_img.shape}")

msmt = render.measure(test_img)
print(f"   Measurement shape: {msmt.shape}")

print("   Running linear_inverse...")
recon_result = linear_inverse(model, render, msmt, stride=50)

# Test 6: Calculate error and save result
print("\n6. Results:")
original = test_img.permute([1, 2, 0]).cpu().numpy()
reconstructed = recon_result[-1]
mse = np.mean((original - reconstructed) ** 2)
print(f"   MSE: {mse:.6f}")

# Save visualization
fig, axs = plt.subplots(1, 2, figsize=(10, 5))
axs[0].imshow(gamma_correct(np.clip(original, 0, 1)))
axs[0].set_title('Original')
axs[0].axis('off')
axs[1].imshow(gamma_correct(np.clip(reconstructed, 0, 1)))
axs[1].set_title(f'Reconstructed ({msmt_dim/im_dim*100:.0f}% measurements)')
axs[1].axis('off')
plt.tight_layout()
plt.savefig('test_result.png', dpi=150)
print("   Saved visualization to test_result.png")

print("\n" + "=" * 50)
print("All tests passed!")
