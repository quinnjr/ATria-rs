// Floyd-Warshall compute shader for ATria algorithm
// This shader performs one iteration of the k-loop on the GPU

struct Params {
    n: u32,      // Matrix dimension
    k: u32,      // Current k iteration
    _pad0: u32,
    _pad1: u32,
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read_write> matrix: array<f32>;

// Workgroup size - 16x16 threads per workgroup
@compute @workgroup_size(16, 16)
fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let i = global_id.x;
    let j = global_id.y;
    let n = params.n;
    let k = params.k;

    // Bounds check
    if (i >= n || j >= n) {
        return;
    }

    // Skip diagonal and when indices equal k
    if (i == j || j == k || i == k) {
        return;
    }

    let curloc = i * n + j;
    let loca = i * n + k;
    let locb = k * n + j;

    let g_i_k = matrix[loca];
    let g_k_j = matrix[locb];
    let product = g_i_k * g_k_j;
    let current = matrix[curloc];

    let evenodd = i + j;

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
