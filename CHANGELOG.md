# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
