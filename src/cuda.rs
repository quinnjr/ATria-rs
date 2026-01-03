//! CUDA-accelerated Floyd-Warshall implementation using cudarc
//!
//! This module provides NVIDIA CUDA acceleration for the Floyd-Warshall
//! algorithm used in ATria. Enable with the `cuda` feature flag.
//!
//! Requires NVIDIA GPU and CUDA toolkit installed.

use cudarc::driver::{CudaDevice, CudaSlice, DeviceRepr, LaunchAsync, LaunchConfig};
use cudarc::nvrtc::Ptx;
use std::sync::Arc;

/// CUDA kernel source code for Floyd-Warshall
const FLOYD_WARSHALL_KERNEL: &str = r#"
extern "C" __global__ void floyd_warshall_kernel(
    float* matrix,
    unsigned int n,
    unsigned int k
) {
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    if (i >= n || j >= n) return;
    if (i == j || j == k || i == k) return;
    
    unsigned int curloc = i * n + j;
    unsigned int loca = i * n + k;
    unsigned int locb = k * n + j;
    
    float g_i_k = matrix[loca];
    float g_k_j = matrix[locb];
    float product = g_i_k * g_k_j;
    float current = matrix[curloc];
    
    unsigned int evenodd = i + j;
    
    if ((evenodd & 1u) == 0u) {
        if (current < product) {
            matrix[curloc] = product;
        }
    } else {
        if (current > product) {
            matrix[curloc] = product;
        }
    }
}
"#;

/// CUDA context for running Floyd-Warshall computations
pub struct CudaContext {
    device: Arc<CudaDevice>,
}

// Manual Debug implementation since CudaDevice doesn't implement Debug
impl std::fmt::Debug for CudaContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CudaContext")
            .field("device", &"<CudaDevice>")
            .finish()
    }
}

impl CudaContext {
    /// Create a new CUDA context, initializing the device and compiling kernel
    pub fn new() -> Option<Self> {
        // Try to get the first CUDA device
        let device = CudaDevice::new(0).ok()?;
        
        log::info!("CUDA device initialized: device 0");
        
        // Compile the kernel using NVRTC
        let ptx = cudarc::nvrtc::compile_ptx(FLOYD_WARSHALL_KERNEL).ok()?;
        
        // Load the PTX module
        device.load_ptx(ptx, "floyd_warshall", &["floyd_warshall_kernel"]).ok()?;
        
        log::info!("CUDA kernel compiled and loaded successfully");
        
        Some(Self { device })
    }

    /// Run Floyd-Warshall algorithm on CUDA GPU
    pub fn floyd_warshall(&self, matrix: &mut [f32], n: usize) {
        // Copy matrix to device
        let mut d_matrix: CudaSlice<f32> = self.device
            .htod_sync_copy(matrix)
            .expect("Failed to copy matrix to CUDA device");

        // Get the kernel function
        let kernel = self.device
            .get_func("floyd_warshall", "floyd_warshall_kernel")
            .expect("Failed to get CUDA kernel function");

        // Configure launch parameters (16x16 thread blocks)
        let block_size = 16u32;
        let grid_x = n.div_ceil(block_size as usize) as u32;
        let grid_y = n.div_ceil(block_size as usize) as u32;
        
        let config = LaunchConfig {
            grid_dim: (grid_x, grid_y, 1),
            block_dim: (block_size, block_size, 1),
            shared_mem_bytes: 0,
        };

        // Run Floyd-Warshall iterations
        for k in 0..n {
            let n_u32 = n as u32;
            let k_u32 = k as u32;
            
            // Launch kernel
            // Safety: kernel parameters match the CUDA function signature
            unsafe {
                kernel.clone().launch(config, (&mut d_matrix, n_u32, k_u32))
                    .expect("Failed to launch CUDA kernel");
            }
        }

        // Synchronize and copy results back
        self.device.synchronize().expect("Failed to synchronize CUDA device");
        self.device
            .dtoh_sync_copy_into(&d_matrix, matrix)
            .expect("Failed to copy results from CUDA device");
    }
}

/// Check if CUDA acceleration is available
pub fn is_cuda_available() -> bool {
    CudaContext::new().is_some()
}
