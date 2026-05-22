# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.4.1] - 2026-05-21

### Fixed

- **`cuda` feature now builds end-to-end on CUDA 13.x.** `cudarc 0.12`
  (the prior pin) capped at CUDA 12.6; on hosts with newer toolkits
  the build script panicked with "Unsupported cuda toolkit version".
  Bumped to `cudarc 0.19`, which supports CUDA 11.4–13.2.
- **`cuda` feature now actually picks a CUDA version.** Previously
  the feature only enabled the cudarc dependency, leaving downstream
  builds to fail with "Must specify one of the following features:
  [cuda-version-from-build-system, …]". The feature now propagates
  `cudarc/cuda-version-from-build-system` so the toolkit version is
  auto-detected from the active `nvcc`.

### Changed

- `src/cuda.rs` migrated to the cudarc 0.19 API surface
  (`CudaContext` + `Stream::launch_builder()` + `clone_htod` /
  `memcpy_dtoh` / `synchronize`); old `LaunchAsync` / `load_ptx`
  path is retired. Module + kernel are compiled at runtime with
  NVRTC as before; algorithmic behaviour is unchanged.

## [1.4.0] - 2026-05-21

### Fixed

- **Output format now matches the canonical C++ ATria byte-for-byte.**
  The previous output emitted `<name>\t<centrality_float>\t\t<rank>` with
  the rank convention inverted (`rank == size` = most central) and the
  centrality value in the second column. This produced a NOA file that
  was incompatible with PluMA's `testPluMA.py` diff harness when run
  against the upstream `corrP.never.ATria.noa.expected` reference. The
  Cytoscape `noa` format upstream tooling expects is now emitted:

  ```
  Name<TAB>Centrality<TAB>Rank
  <name>   #<rank> <name>   <rank>     // for ranked nodes (rank 1 = most central)
  <name>   <name>            NR          // for nodes never selected as max-pay
  ```

- **Tied nodes now correctly share a rank.** When two or more nodes have
  the same `|pay|` at iteration step *k*, the C++ original assigns all
  of them rank *k* and advances `currentrank` by the tie count, so the
  next-ranked node skips ahead. The previous Rust implementation
  recorded only the first tied node's pay value and dropped the rest to
  "unranked" — losing ranking information for any iteration with ties.
  ATria-rs now matches the upstream tie semantics exactly: in the
  reference corrP.never network, the three-way tie at iteration 4
  (Bifidobacterium.01 / Peptoniphilus.02 / Peptoniphilus.03) is now
  recorded as `#4`, `#4`, `#4` with the next-ranked node at `#7`.

### Added

- New `ranks: Vec<u32>` field on `ATriaPlugin` to track iteration ranks
  separately from the `output: Vec<f32>` pay values (which are still
  populated for diagnostic purposes).
- `run()` is now defensive about hand-constructed `ATriaPlugin` instances
  that skip `input()`: it ensures `output` and `ranks` are sized to
  `gsize` before assigning to them.

### Changed

- The bundled `tests/corrP.never.noa.expected` is now the canonical C++
  ATria reference output (same fixture used by the upstream
  `movingpictures83/ATria` example).
- `tests::it_works` therefore now asserts byte-for-byte equality
  against the upstream output rather than against ATria-rs's previous
  divergent format.

### Verification

- `cargo test --release --lib` → 3/3 pass.
- `python3 testPluMA.py ATria` (PluMA's official diff harness, against
  the upstream `corrP.never.ATria.noa.expected`) → **Passing Rate: 100.0%**.
- Standalone diff vs upstream reference is empty (sort-identical).

## [1.3.0] - 2026-01-03

### Added

- **PluMA FFI Exports**: New `pluma_ffi` module for direct plugin loading
  - `plugin_create()` - Create new ATriaPlugin instance
  - `plugin_destroy()` - Clean up plugin instance  
  - `plugin_input()` - Read CSV adjacency matrix
  - `plugin_run()` - Execute ATria centrality algorithm
  - `plugin_output()` - Write NOA file for Cytoscape
- Plugin can now be loaded directly by PluMA via dlopen/dlsym
- Uses `#[unsafe(no_mangle)]` for Rust 2024 edition compatibility

### Changed

- Added `pluma` keyword to crate metadata
- Added repository URL to Cargo.toml

## [1.2.0] - 2026-01-03

### Added

- **NVIDIA CUDA Support**: New `cuda` feature for NVIDIA GPU acceleration
  - CUDA kernel for Floyd-Warshall algorithm (`src/kernels/floyd_warshall.cu`)
  - Runtime PTX compilation via NVRTC
  - `CudaContext` for device management
- `ComputeBackend::Cuda` variant for explicit CUDA backend selection
- `ATriaPlugin::is_cuda_available()` method to check CUDA support
- Auto backend now prefers CUDA > GPU > CPU

### Changed

- Renamed internal GPU context references for clarity (wgpu vs CUDA)
- Updated documentation to cover both GPU and CUDA options

### Dependencies (cuda feature)

- cudarc 0.12 - CUDA runtime and NVRTC bindings

## [1.1.0] - 2026-01-03

### Added

- **GPU Acceleration**: Optional GPU compute shader support via wgpu
  - WGSL compute shader for Floyd-Warshall algorithm (`src/shaders/floyd_warshall.wgsl`)
  - GPU context management with automatic device selection
  - Cross-platform support (Vulkan, Metal, DX12)
- `ComputeBackend` enum for runtime backend selection (Cpu, Gpu, Auto)
- `ATriaPlugin::with_backend()` constructor for specifying compute backend
- `ATriaPlugin::set_backend()` and `ATriaPlugin::set_use_gpu()` methods
- `ATriaPlugin::is_gpu_available()` method to check GPU support
- `ATriaPlugin::effective_backend()` method to get actual backend in use
- New Cargo feature flag: `gpu` for optional GPU dependencies

### Changed

- **Breaking**: Updated to Rust 2024 edition (requires Rust 1.85+)
- Replaced `write!` with `writeln!` macros for cleaner code
- Refactored loops to use idiomatic iterators with `enumerate()`
- Used `div_ceil()` instead of manual ceiling division
- Collapsed nested `if` statements for better readability
- Removed redundant imports flagged by clippy

### Dependencies (gpu feature)

- wgpu 0.20 - Cross-platform GPU compute
- pollster 0.3 - Async runtime for GPU initialization
- bytemuck 1.16 - Safe casting for GPU buffers

## [1.0.0] - 2026-01-03

### Added

- Criterion benchmarking suite (v0.8.0) for performance testing
  - `cpu_floyd` benchmark for Floyd-Warshall algorithm at various matrix sizes
  - `atria_full_run` benchmark for complete algorithm execution
  - `input_parsing` benchmark for CSV parsing performance
  - `pay_calculation` benchmark for pay calculation loop
- Comprehensive README documentation with usage examples and performance metrics

### Changed

- Optimized Floyd-Warshall algorithm with unsafe pointer arithmetic (~17% faster)
- Replaced modulo operation with bitwise AND for parity checks
- Pre-computed index offsets in hot loops to reduce redundant calculations
- Added buffered I/O (BufWriter) for file output operations
- Improved input parsing with pre-allocated vector capacity

### Performance

- Floyd-Warshall (252×252): 15.6 ms → 12.9 ms (−17.6%)
- Full ATria run (126 bacteria): 508 ms → 425 ms (−16.4%)
- Input parsing: 190 µs → 183 µs (−3.4%)

## [0.1.0] - 2020

### Added

- Initial Rust implementation of the ATria centrality algorithm
- CSV input parsing for network data
- NOA file output for Cytoscape compatibility
- PluMAPlugin trait implementation
- Unit tests for algorithm correctness
