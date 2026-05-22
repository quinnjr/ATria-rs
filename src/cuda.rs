//! CUDA-accelerated Floyd-Warshall implementation using cudarc 0.19+.
//!
//! Provides NVIDIA CUDA acceleration for the modified Floyd-Warshall
//! step at the heart of the ATria algorithm. Enable with the `cuda`
//! Cargo feature. Requires an NVIDIA GPU exposed to the host and the
//! CUDA toolkit's `nvcc` reachable through `$CUDA_PATH` / `$CUDA_ROOT`
//! (cudarc compiles the kernel with NVRTC at runtime, and the build
//! script picks up the toolkit version automatically via the
//! `cuda-version-from-build-system` feature that atria-rs's `cuda`
//! feature propagates).

use cudarc::driver::{CudaContext as CudarcContext, CudaSlice, LaunchConfig, PushKernelArg};
use cudarc::nvrtc::Ptx;
use std::sync::Arc;

/// CUDA kernel source for the modified Floyd-Warshall step. One k
/// iteration per launch; the host loop drives k = 0..n.
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
    unsigned int loca   = i * n + k;
    unsigned int locb   = k * n + j;

    float g_i_k   = matrix[loca];
    float g_k_j   = matrix[locb];
    float product = g_i_k * g_k_j;
    float current = matrix[curloc];

    unsigned int evenodd = i + j;
    if ((evenodd & 1u) == 0u) {
        if (current < product) matrix[curloc] = product;
    } else {
        if (current > product) matrix[curloc] = product;
    }
}
"#;

/// CUDA context for the ATria Floyd-Warshall step.
///
/// Owns the cudarc context + a compiled module containing the
/// `floyd_warshall_kernel`. Construction is fallible (no GPU, no
/// driver, or NVRTC compile failure all return `None`); callers should
/// fall back to the CPU implementation in that case.
pub struct CudaContext {
    ctx: Arc<CudarcContext>,
    module: Arc<cudarc::driver::CudaModule>,
}

impl std::fmt::Debug for CudaContext {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("CudaContext")
            .field("device_ordinal", &self.ctx.ordinal())
            .finish()
    }
}

impl CudaContext {
    /// Initialize CUDA device 0 and compile the kernel.
    ///
    /// Returns `None` if any step fails — call `is_cuda_available()` from
    /// a higher layer to surface the failure to the user.
    pub fn new() -> Option<Self> {
        let ctx = CudarcContext::new(0).ok()?;
        log::info!(
            "CUDA device initialized: ordinal {}, compute capability {:?}",
            ctx.ordinal(),
            ctx.compute_capability().ok()
        );

        let ptx: Ptx = cudarc::nvrtc::compile_ptx(FLOYD_WARSHALL_KERNEL).ok()?;
        let module = ctx.load_module(ptx).ok()?;
        log::info!("CUDA kernel compiled and loaded");

        Some(Self { ctx, module })
    }

    /// Run the iterative Floyd-Warshall step on the GPU. Mutates
    /// `matrix` in place; assumes it's the flat row-major representation
    /// of an `n × n` matrix.
    pub fn floyd_warshall(&self, matrix: &mut [f32], n: usize) {
        debug_assert_eq!(matrix.len(), n * n, "matrix size must be n*n");

        let stream = self.ctx.default_stream();

        // Copy matrix to device.
        let mut d_matrix: CudaSlice<f32> = stream
            .clone_htod(matrix)
            .expect("htod: failed to copy matrix to device");

        // Load the kernel function once and reuse it across iterations.
        let kernel = self
            .module
            .load_function("floyd_warshall_kernel")
            .expect("load_function: floyd_warshall_kernel");

        // 16 x 16 thread blocks; one block tile per (i / 16, j / 16).
        let block_size = 16u32;
        let grid_x = n.div_ceil(block_size as usize) as u32;
        let grid_y = n.div_ceil(block_size as usize) as u32;
        let cfg = LaunchConfig {
            grid_dim: (grid_x, grid_y, 1),
            block_dim: (block_size, block_size, 1),
            shared_mem_bytes: 0,
        };

        let n_u32 = n as u32;
        for k in 0..n {
            let k_u32 = k as u32;
            let mut builder = stream.launch_builder(&kernel);
            builder.arg(&mut d_matrix);
            builder.arg(&n_u32);
            builder.arg(&k_u32);
            // Safety: kernel signature matches the (float*, uint, uint)
            // args pushed above; the device buffer is not aliased
            // outside this scope for the duration of the launch.
            unsafe { builder.launch(cfg) }.expect("kernel launch failed");
        }

        stream.synchronize().expect("CUDA stream synchronize");
        stream
            .memcpy_dtoh(&d_matrix, matrix)
            .expect("dtoh: failed to copy results from device");
    }
}

/// Quick is-CUDA-available probe. Allocates a context and immediately
/// drops it; intended for one-shot "should I use Cuda or fall back to
/// Cpu" decisions during `effective_backend()` resolution.
pub fn is_cuda_available() -> bool {
    CudaContext::new().is_some()
}
