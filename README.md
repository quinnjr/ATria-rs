# ATria-rs

A reimplementation of the ATria algorithm in Rust.

[![Rust](https://img.shields.io/badge/rust-2024_edition-orange.svg)](https://www.rust-lang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE.md)

## Overview

Library for the Ablatio Triadum (ATria) centrality algorithm (Cickovski et al, 2015, 2017).

ATria can run on signed and weighted networks and produces a list of central nodes as both screen output and as a NOde Attribute (NOA) file for Cytoscape. The NOA file can subsequently be imported into Cytoscape resulting in centrality values becoming node attributes and enabling further analysis and visualization based on these values.

## Features

- **High Performance**: Optimized CPU implementation with unsafe pointer arithmetic (~17% faster than baseline)
- **GPU Acceleration**: Optional cross-platform GPU support via wgpu (Vulkan/Metal/DX12)
- **CUDA Support**: Optional NVIDIA CUDA acceleration for maximum performance
- **Cross-Platform**: Works on Linux, macOS, and Windows
- **Configurable Backend**: Runtime selection between CPU, GPU, and CUDA computation
- **Modern Rust**: Built with Rust 2024 edition

## Requirements

- Rust 1.85+ (2024 edition)
- For GPU support: Compatible graphics driver (Vulkan, Metal, or DX12)
- For CUDA support: NVIDIA GPU with CUDA toolkit installed

## Input Format

The input network should be specified in CSV format with nodes as rows and columns and entry (i, j) representing the weight of the edge from node i to node j.

## Output Format

The output is a NOA file, with both centrality value and rank as attributes. Larger magnitude values indicate higher centrality for both centrality and rank. This is typically more convenient for visualization, etc.

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
atria-rs = { git = "https://github.com/quinnjr/ATria-rs" }
```

### With GPU Support (Cross-Platform)

To enable wgpu GPU acceleration (Vulkan/Metal/DX12):

```toml
[dependencies]
atria-rs = { git = "https://github.com/quinnjr/ATria-rs", features = ["gpu"] }
```

### With CUDA Support (NVIDIA)

To enable NVIDIA CUDA acceleration:

```toml
[dependencies]
atria-rs = { git = "https://github.com/quinnjr/ATria-rs", features = ["cuda"] }
```

### With Both GPU and CUDA

```toml
[dependencies]
atria-rs = { git = "https://github.com/quinnjr/ATria-rs", features = ["gpu", "cuda"] }
```

## Usage

### Basic Usage (CPU)

```rust
use ATriaPlugin::ATriaPlugin;
use pluma_plugin_trait::PluMAPlugin;

fn main() {
    let mut plugin = ATriaPlugin::default();
    
    // Load input CSV
    plugin.input("path/to/network.csv".to_string()).unwrap();
    
    // Run ATria algorithm
    plugin.run().unwrap();
    
    // Write output NOA file
    plugin.output("path/to/output.noa".to_string()).unwrap();
}
```

### With GPU/CUDA Acceleration

```rust
use ATriaPlugin::{ATriaPlugin, ComputeBackend};
use pluma_plugin_trait::PluMAPlugin;

fn main() {
    // Create plugin with specific backend
    let mut plugin = ATriaPlugin::with_backend(ComputeBackend::Cuda);
    
    // Or use auto-detection (prefers CUDA > GPU > CPU)
    // let mut plugin = ATriaPlugin::with_backend(ComputeBackend::Auto);
    
    // Check available backends
    println!("CUDA available: {}", plugin.is_cuda_available());
    println!("GPU available: {}", plugin.is_gpu_available());
    println!("Effective backend: {:?}", plugin.effective_backend());
    
    plugin.input("path/to/network.csv".to_string()).unwrap();
    plugin.run().unwrap();
    plugin.output("path/to/output.noa".to_string()).unwrap();
}
```

### Compute Backends

| Backend | Description | Feature |
|---------|-------------|---------|
| `ComputeBackend::Cpu` | CPU-only computation (default) | - |
| `ComputeBackend::Gpu` | wgpu GPU acceleration (Vulkan/Metal/DX12) | `gpu` |
| `ComputeBackend::Cuda` | NVIDIA CUDA acceleration | `cuda` |
| `ComputeBackend::Auto` | Best available (CUDA > GPU > CPU) | - |

## Building

```bash
# CPU-only build
cargo build --release

# With wgpu GPU support
cargo build --release --features gpu

# With CUDA support
cargo build --release --features cuda

# With all GPU backends
cargo build --release --features "gpu cuda"
```

## Testing

```bash
cargo test
```

## Benchmarking

This project uses [Criterion](https://github.com/bheisler/criterion.rs) v0.8 for benchmarking. Run benchmarks with:

```bash
cargo bench
```

### Performance

The implementation includes several optimizations for the core Floyd-Warshall algorithm:

- Unsafe pointer arithmetic to eliminate bounds checking in hot loops
- Pre-computed index offsets to reduce redundant calculations
- Bitwise operations for parity checks
- Buffered I/O for file operations
- Optional GPU/CUDA acceleration for large matrices

Benchmark results on a 126-bacteria network (252×252 matrix):

| Operation | CPU Time |
|-----------|----------|
| Floyd-Warshall (252×252) | ~12.9 ms |
| Full ATria run | ~425 ms |
| Input parsing | ~183 µs |

## GPU Requirements

### wgpu (Cross-Platform)

When using the `gpu` feature, the following backends are supported:

- **Vulkan** (Linux, Windows)
- **Metal** (macOS)
- **DX12** (Windows)

### CUDA (NVIDIA)

When using the `cuda` feature:

- NVIDIA GPU (Compute Capability 3.5+)
- CUDA Toolkit 11.0+ installed
- NVIDIA driver 450.80.02+ (Linux) or 452.39+ (Windows)

## Dependencies

### Core Dependencies

| Crate | Version | Purpose |
|-------|---------|---------|
| csv | 1.1 | CSV file parsing |
| log | 0.4 | Logging framework |
| rayon | 1.4 | Parallel processing |

### GPU Dependencies (gpu feature)

| Crate | Version | Purpose |
|-------|---------|---------|
| wgpu | 0.20 | Cross-platform GPU compute |
| pollster | 0.3 | Async runtime |
| bytemuck | 1.16 | Safe buffer casting |

### CUDA Dependencies (cuda feature)

| Crate | Version | Purpose |
|-------|---------|---------|
| cudarc | 0.12 | CUDA runtime bindings |

## Project Structure

```
ATria-rs/
├── src/
│   ├── lib.rs              # Main library with ATriaPlugin
│   ├── gpu.rs              # wgpu GPU module (optional)
│   ├── cuda.rs             # CUDA module (optional)
│   ├── shaders/
│   │   └── floyd_warshall.wgsl  # WGSL compute shader
│   └── kernels/
│       └── floyd_warshall.cu    # CUDA kernel
├── benches/
│   └── atria_benchmark.rs  # Criterion benchmarks
├── tests/
│   ├── corrP.never.csv     # Test input data
│   └── corrP.never.noa.expected  # Expected output
├── Cargo.toml
├── CHANGELOG.md
├── LICENSE.md
└── README.md
```

## References

- Original C++ implementation: [movingpictures83/ATria](https://github.com/movingpictures83/ATria)
- Cickovski, T., et al. (2015, 2017) - ATria centrality algorithm papers

## License

MIT
