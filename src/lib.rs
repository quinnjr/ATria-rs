// Copyright (C) 2020 Joseph R. Quinn
// SPDX-License-Identifier: MIT

#![allow(non_snake_case)]

//! # ATria-rs
//!
//! Library for the Ablatio Triadum (ATria) centrality algorithm
//! (Cickovski et al, 2015, 2017).
//!
//! ATria can run on signed and weighted networks and produces a list of
//! central nodes as both screen output and as a NOde Attribute (NOA)
//! file for Cytoscape. The NOA file can subsequently be imported into
//! Cytoscape resulting in centrality values becoming node attributes and
//! enabling further analysis and visualization based on these values.
//!
//! The input network should be specified in CSV format with nodes as rows
//! and columns and entry (i, j) representing the weight of the edge from
//! node i to node j.
//!
//! The output is the NOA file, with both centrality value and rank as
//! attributes. Larger magnitude values indicate higher centrality for
//! both centrality and rank. This is typically more convenient for
//! visualization, etc.
//!
//! ## GPU Acceleration
//!
//! Enable the `gpu` feature for GPU-accelerated Floyd-Warshall computation:
//!
//! ```toml
//! [dependencies]
//! atria-rs = { version = "1.0", features = ["gpu"] }
//! ```
//!
//! Then configure the plugin to use GPU:
//!
//! ```ignore
//! let mut plugin = ATriaPlugin::default();
//! plugin.set_use_gpu(true);
//! ```
//!
//! Original code for the C++ version of this library may be
//! found [here](https://github.com/movingpictures83/ATria).

use std::fs::File;
use std::io::prelude::*;
use std::io::BufWriter;

use log::*;
use pluma_plugin_trait::PluMAPlugin;

// GPU module (conditional compilation)
#[cfg(feature = "gpu")]
pub mod gpu;

/// Standard replacement for crate-level `std::result::Result<(), Box<dyn std::error::Error>>`
type Result<T = ()> = std::result::Result<T, Box<dyn std::error::Error>>;

/// Configuration for compute backend selection
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ComputeBackend {
    /// Use CPU for computation (default)
    #[default]
    Cpu,
    /// Use GPU for computation (requires `gpu` feature)
    Gpu,
    /// Automatically select best available backend
    Auto,
}

#[derive(Debug)]
pub struct ATriaPlugin {
    /// Number of bacteria (GSIZE in C++)
    pub gsize: usize,
    /// Vector of the bacteria types in the CSV file.
    pub bacteria: Vec<String>,
    /// The original matrix being worked on by the ATria algorithm (2N x 2N).
    pub orig_graph: Vec<f32>,
    /// Output centrality values (stores pay values, NOT ranks)
    pub output: Vec<f32>,
    /// Compute backend selection
    backend: ComputeBackend,
    /// GPU context (lazily initialized when gpu feature is enabled)
    #[cfg(feature = "gpu")]
    gpu_context: Option<gpu::GpuContext>,
}

// Manual Default impl required due to conditional #[cfg(feature = "gpu")] field
#[allow(clippy::derivable_impls)]
impl Default for ATriaPlugin {
    fn default() -> Self {
        ATriaPlugin {
            gsize: 0,
            bacteria: Vec::new(),
            orig_graph: Vec::new(),
            output: Vec::new(),
            backend: ComputeBackend::default(),
            #[cfg(feature = "gpu")]
            gpu_context: None,
        }
    }
}

impl ATriaPlugin {
    /// Create a new ATriaPlugin with the specified compute backend
    pub fn with_backend(backend: ComputeBackend) -> Self {
        let mut plugin = Self::default();
        plugin.set_backend(backend);
        plugin
    }

    #[inline(always)]
    fn size(&self) -> usize {
        self.gsize
    }

    /// Set the compute backend to use
    pub fn set_backend(&mut self, backend: ComputeBackend) {
        self.backend = backend;

        #[cfg(feature = "gpu")]
        {
            // Initialize GPU context if needed
            if matches!(backend, ComputeBackend::Gpu | ComputeBackend::Auto)
                && self.gpu_context.is_none()
            {
                self.gpu_context = gpu::GpuContext::new();
                if self.gpu_context.is_some() {
                    info!("GPU context initialized successfully");
                } else {
                    warn!("Failed to initialize GPU context, falling back to CPU");
                }
            }
        }

        #[cfg(not(feature = "gpu"))]
        {
            if matches!(backend, ComputeBackend::Gpu) {
                warn!("GPU backend requested but 'gpu' feature is not enabled, using CPU");
            }
        }
    }

    /// Get the current compute backend
    pub fn backend(&self) -> ComputeBackend {
        self.backend
    }

    /// Enable GPU acceleration (convenience method)
    pub fn set_use_gpu(&mut self, use_gpu: bool) {
        self.set_backend(if use_gpu { ComputeBackend::Gpu } else { ComputeBackend::Cpu });
    }

    /// Check if GPU is available for computation
    #[cfg(feature = "gpu")]
    pub fn is_gpu_available(&self) -> bool {
        self.gpu_context.is_some()
    }

    /// Check if GPU is available for computation
    #[cfg(not(feature = "gpu"))]
    pub fn is_gpu_available(&self) -> bool {
        false
    }

    /// Get the effective backend that will be used for computation
    pub fn effective_backend(&self) -> ComputeBackend {
        match self.backend {
            ComputeBackend::Cpu => ComputeBackend::Cpu,
            ComputeBackend::Gpu => {
                if self.is_gpu_available() {
                    ComputeBackend::Gpu
                } else {
                    ComputeBackend::Cpu
                }
            }
            ComputeBackend::Auto => {
                if self.is_gpu_available() {
                    ComputeBackend::Gpu
                } else {
                    ComputeBackend::Cpu
                }
            }
        }
    }
}

/// Modified Floyd-Warshall algorithm for ATria - optimized version
/// Uses unsafe pointer arithmetic to avoid bounds checking in hot loops
#[inline]
pub fn cpu_floyd(g: &mut [f32], n: usize) {
    // Pre-compute k*n once per outer loop iteration
    // Use unsafe to avoid bounds checks in the innermost loop
    let ptr = g.as_mut_ptr();

    for k in 0..n {
        let k_row_offset = k * n;

        for i in 0..n {
            let i_row_offset = i * n;

            // SAFETY: All indices are within bounds since i, j, k < n
            // and the array has n*n elements
            unsafe {
                let g_i_k = *ptr.add(i_row_offset + k);

                for j in 0..n {
                    // Skip diagonal and when j == k
                    if i == j || j == k {
                        continue;
                    }

                    let curloc = i_row_offset + j;
                    let g_k_j = *ptr.add(k_row_offset + j);
                    let product = g_i_k * g_k_j;
                    let current = *ptr.add(curloc);
                    let evenodd = i + j;

                    // Use bitwise AND for parity check (faster than modulo)
                    if (evenodd & 1) == 0 {
                        // Even: maximize
                        if current < product {
                            *ptr.add(curloc) = product;
                        }
                    } else {
                        // Odd: minimize (for negative paths)
                        if current > product {
                            *ptr.add(curloc) = product;
                        }
                    }
                }
            }
        }
    }
}

impl PluMAPlugin for ATriaPlugin {
    /// Create a 2Nx2N adjacency matrix from the input CSV file.
    fn input(&mut self, file_path: String) -> Result {
        let mut reader = csv::Reader::from_path(&file_path).expect("Unable to open CSV file");

        // First pass: count rows to determine GSIZE
        {
            let headers = reader
                .headers()
                .expect("Unable to read CSV headers")
                .clone();

            // Pre-allocate with expected capacity
            self.bacteria.reserve(headers.len() - 1);
            for header in headers.iter().skip(1) {
                self.bacteria.push(header.to_string());
            }
        }

        self.gsize = self.bacteria.len();
        let gsize = self.gsize;

        // Allocate 2N x 2N matrix
        let matrix_size = (gsize * 2) * (gsize * 2);
        self.orig_graph = vec![0.0f32; matrix_size];
        self.output = vec![0.0f32; gsize];

        // Re-read to populate matrix
        let mut reader = csv::Reader::from_path(&file_path).expect("Unable to open CSV file");
        let stride = 2 * gsize;

        for (row_count, result) in reader.records().enumerate() {
            let row = result.expect("Unable to read CSV row");
            let bac1 = row_count;
            let bac1_2 = bac1 * 2;
            let bac1_2_1 = bac1_2 + 1;

            for i in 1..row.len() {
                let bac2 = i - 1;
                let bac2_2 = bac2 * 2;
                let bac2_2_1 = bac2_2 + 1;

                // Pre-compute indices
                let idx_00 = bac1_2 * stride + bac2_2;
                let idx_11 = bac1_2_1 * stride + bac2_2_1;
                let idx_10 = bac1_2_1 * stride + bac2_2;
                let idx_01 = bac1_2 * stride + bac2_2_1;

                if bac1 != bac2 {
                    let weight: f32 = row[i].parse().expect("Unable to parse weight");

                    if weight > 0.0 {
                        self.orig_graph[idx_00] = weight;
                        self.orig_graph[idx_11] = weight;
                        // idx_10 and idx_01 already 0 from initialization
                    } else if weight < 0.0 {
                        self.orig_graph[idx_10] = weight;
                        self.orig_graph[idx_01] = weight;
                        // idx_00 and idx_11 already 0 from initialization
                    }
                    // weight == 0: all already initialized to 0
                } else {
                    // Diagonal: start at 1 because they are starting verts
                    self.orig_graph[idx_00] = 1.0;
                    self.orig_graph[idx_11] = 1.0;
                    // idx_10 and idx_01 already 0 from initialization
                }
            }
        }

        Ok(())
    }

    /// Run the ATria algorithm over the input data.
    fn run(&mut self) -> Result {
        let effective_backend = self.effective_backend();
        info!("Running ATria with {:?} backend", effective_backend);

        let gsize = self.gsize;
        let n = gsize * 2;  // Matrix dimension
        let stride = n;
        let matrix_len = n * n;

        // Working copy of graph for Floyd-Warshall - allocate once
        let mut h_g = vec![0.0f32; matrix_len];
        // Pay values for each bacterium
        let mut h_pay = vec![0.0f32; gsize];
        // Pre-allocate maxnodes vector
        let mut maxnodes = Vec::with_capacity(gsize);

        for _ in 0..gsize {
            // Copy original graph for computation - use fast copy
            h_g.copy_from_slice(&self.orig_graph);

            // Run modified Floyd-Warshall using selected backend
            #[cfg(feature = "gpu")]
            {
                if matches!(effective_backend, ComputeBackend::Gpu) {
                    if let Some(ref gpu_ctx) = self.gpu_context {
                        gpu_ctx.floyd_warshall(&mut h_g, n);
                    } else {
                        cpu_floyd(&mut h_g, n);
                    }
                } else {
                    cpu_floyd(&mut h_g, n);
                }
            }

            #[cfg(not(feature = "gpu"))]
            {
                cpu_floyd(&mut h_g, n);
            }

            // Calculate pay for each bacterium - using iterators for better SIMD
            for (i, pay) in h_pay.iter_mut().enumerate() {
                let row_start = (i * 2) * stride;
                *pay = h_g[row_start..row_start + n].iter().sum::<f32>() - 1.0;
            }

            // Find node(s) with maximum pay - find FIRST maximum (matching original behavior)
            let mut mnode = 0usize;
            let mut maxpay = -1.0f32;
            for (i, &pay) in h_pay.iter().enumerate() {
                if pay.abs() > maxpay {
                    mnode = i;
                    maxpay = pay.abs();
                }
            }

            if maxpay == 0.0 {
                break;
            }

            // Find all nodes with same max pay
            maxnodes.clear();
            maxnodes.push(mnode);
            for (i, &pay) in h_pay.iter().enumerate() {
                if i != mnode && pay.abs() == maxpay {
                    maxnodes.push(i);
                }
            }

            // Process each max node
            for &maxnode in &maxnodes {
                info!("Node with highest pay: {}: {}", self.bacteria[maxnode], h_pay[maxnode]);

                // Only record centrality for the first node (mnode), not ties
                if maxnode == mnode {
                    self.output[maxnode] = h_pay[maxnode];
                }

                let maxnode_2 = maxnode * 2;
                let maxnode_2_1 = maxnode_2 + 1;
                let maxnode_row_even = maxnode_2 * stride;
                let maxnode_row_odd = maxnode_2_1 * stride;

                // Non-GPU Triad Removal
                for i in 0..n {
                    if (i / 2) != maxnode {
                        let edge_even = self.orig_graph[maxnode_row_even + i] != 0.0;
                        let edge_odd = self.orig_graph[maxnode_row_odd + i] != 0.0;

                        if edge_even || edge_odd {
                            let i_row = i * stride;

                            for j in (i + 1)..n {
                                if (j / 2) != maxnode {
                                    let connected_j = self.orig_graph[maxnode_row_even + j] != 0.0
                                        || self.orig_graph[maxnode_row_odd + j] != 0.0;

                                    if connected_j && self.orig_graph[i_row + j] != 0.0 {
                                        self.orig_graph[i_row + j] = 2.0;
                                        self.orig_graph[j * stride + i] = 2.0;
                                    }
                                }
                            }

                            if edge_even {
                                self.orig_graph[maxnode_row_even + i] = 2.0;
                                self.orig_graph[i_row + maxnode_2] = 2.0;
                            }

                            if edge_odd {
                                self.orig_graph[maxnode_row_odd + i] = 2.0;
                                self.orig_graph[i_row + maxnode_2_1] = 2.0;
                            }
                        }
                    }
                }

                // Sweep through and remove marked edges - vectorized
                self.orig_graph.iter_mut().for_each(|v| {
                    if *v == 2.0 {
                        *v = 0.0;
                    }
                });
            }
        }

        Ok(())
    }

    /// Write the results of the ATria calulations to a NOA file.
    fn output(&mut self, file_path: String) -> Result {
        // Use buffered writer for better I/O performance
        let file = File::create(file_path).expect("Unable to open output file location");
        let mut output_file = BufWriter::new(file);

        // Sort by absolute value of output (descending) using bubble sort
        // (maintains compatibility with original algorithm output)
        let size = self.size();
        for i in (0..size).rev() {
            for j in 0..i {
                if self.output[j].abs() < self.output[j + 1].abs() {
                    self.output.swap(j, j + 1);
                    self.bacteria.swap(j, j + 1);
                }
            }
        }

        writeln!(output_file, "Name\tCentrality\tRank")
            .expect("Unable to write headers to output file");

        for i in 0..size {
            self.output[i] = self.output[i].abs();

            writeln!(
                output_file,
                "{}\t{}\t\t{}",
                self.bacteria[i],
                self.output[i],
                size - i
            )
            .expect("Unable to write to output file");
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_should_load_bacteria() {
        let mut plugin = ATriaPlugin::default();

        plugin
            .input("./tests/corrP.never.csv".to_string())
            .expect("Failed to read CSV file");

        assert_eq!(126, plugin.bacteria.len());
    }

    #[test]
    fn it_can_run() {
        let mut plugin = ATriaPlugin::default();

        plugin.gsize = 2;

        plugin.orig_graph = vec![
            1.0, 0.0, 0.5, 0.0,
            0.0, 1.0, 0.0, 0.5,
            0.5, 0.0, 1.0, 0.0,
            0.0, 0.5, 0.0, 1.0,
        ];

        plugin.bacteria = vec![
            String::from("Test Bac 1"),
            String::from("Test Bac 2"),
        ];

        plugin.output.resize(2, 0.0);

        assert!(plugin.run().is_ok());
    }

    #[test]
    fn it_works() {
        let mut plugin = ATriaPlugin::default();

        plugin
            .input("./tests/corrP.never.csv".to_string())
            .expect("Failed to read CSV file...");

        plugin.run().expect("Failed to run ATria...");

        plugin
            .output("./tests/corrP.never.noa".to_string())
            .expect("Failed to write NOA file...");

        let mut expected = File::open("./tests/corrP.never.noa.expected")
            .expect("Failed to open expected output file...");

        let mut expected_content = String::new();
        expected
            .read_to_string(&mut expected_content)
            .expect("Failed to read expected output file.");

        let mut actual =
            File::open("./tests/corrP.never.noa").expect("Failed to open generated output file...");

        let mut actual_content = String::new();
        actual
            .read_to_string(&mut actual_content)
            .expect("Failed to read actual output to file.");

        assert_eq!(expected_content, actual_content);
    }
}
