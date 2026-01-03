// Floyd-Warshall CUDA kernel for ATria algorithm
// This kernel performs one iteration of the k-loop on the GPU

extern "C" __global__ void floyd_warshall_kernel(
    float* matrix,
    unsigned int n,
    unsigned int k
) {
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    // Bounds check
    if (i >= n || j >= n) {
        return;
    }
    
    // Skip diagonal and when indices equal k
    if (i == j || j == k || i == k) {
        return;
    }
    
    unsigned int curloc = i * n + j;
    unsigned int loca = i * n + k;
    unsigned int locb = k * n + j;
    
    float g_i_k = matrix[loca];
    float g_k_j = matrix[locb];
    float product = g_i_k * g_k_j;
    float current = matrix[curloc];
    
    unsigned int evenodd = i + j;
    
    // Even indices: maximize, Odd indices: minimize
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
