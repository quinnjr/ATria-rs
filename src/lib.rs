// Copyright (C) 2020 Joseph R. Quinn
// SPDX-License-Identifier: MIT

#![allow(non_snake_case)]

/// Library for the Ablatio Triadum (ATria) centrality algorithm
/// (Cickovski et al, 2015, 2017).
///
/// ATria can run on signed and weighted networks and produces a list of
/// central nodes as both screen output and as a NOde Attribute (NOA)
/// file for Cytoscape. The NOA file can subsequently be imported into
/// Cytoscape resulting in centrality values becoming node attributes and
/// enabling further analysis and visualization based on these values.
///
/// The input network should be specified in CSV format with nodes as rows
/// and columns and entry (i, j) representing the weight of the edge from
/// node i to node j.
///
/// The output is the NOA file, with both centrality value and rank as
/// attributes. Larger magnitude values indicate higher centrality for
/// both centrality and rank. This is typically more convenient for
/// visualization, etc.
///
/// Original code for the C++ version of this library may be
/// found (here)[https://github.com/movingpictures83/ATria].
use std::fs::File;
use std::io::prelude::*;
use std::sync::{Arc, Mutex};

use csv;
use libm::fabsf;
use log::*;
use pluma_plugin_trait::PluMAPlugin;
use rayon::prelude::*;

// pub mod ffi;

const PARTITION_SIZE: usize = 25;

/// Standard replacement for crate-level `std::result::Result<(), Box<dyn std::error::Error>>`
type Result<T = ()> = std::result::Result<T, Box<dyn std::error::Error>>;

/// Calculate memory offset of a value in a 1D array as
/// if it was a 2D array.
#[inline]
fn vec_offset(i: usize, j: usize, sz: usize) -> usize {
    i * sz + j
}

// fn modified_floyd_warshall(mut graph: &Vec, size: usize) {
//     for k in 0..size {
//         for i in 0..size {
//             for j in 0..size {
//                 if i != j && j != k {
//                     let even_odd = i + j;
//                     let current_node = graph[i * size + j];
//                     let node_a = graph[i * size + k];
//                     let node_b = graph[k * size + j];

//                     if (even_odd % 2 == 0) &&
//                         (current_node < (node_a * node_b)) ||
//                         (even_odd % 2 == 1) &&
//                         (current_node > (node_a * node_b)) {
//                         graph[i * size + j] = graph[i * size + k] * graph[k * size + j];
//                     }
//                 }
//             }
//         }
//     }
// }

#[derive(Debug)]
pub struct ATriaPlugin {
    /// Vector of the bacteria types in the CSV file.
    pub bacteria: Vec<String>,
    /// The original matrix being worked on by the ATria algorithm.
    pub orig_graph: Vec<f32>,
    ///
    pub output: Vec<f32>,
}

impl Default for ATriaPlugin {
    fn default() -> Self {
        ATriaPlugin {
            bacteria: Vec::new(),
            orig_graph: Vec::new(),
            output: Vec::new(),
        }
    }
}

impl ATriaPlugin {
    #[inline]
    fn size(&self) -> usize {
        self.bacteria.len()
    }
}

impl PluMAPlugin for ATriaPlugin {
    /// Create a NxN adjacency matrix from the input CSV file.
    /// Space complextiy is expected to be O(n^2).
    fn input(&mut self, file_path: String) -> Result {
        let mut reader = csv::Reader::from_path(file_path).expect("Unable to open CSV file");

        {
            let headers = reader
                .headers()
                .expect("Unable to read CSV headers")
                .clone();
            for header in headers.iter().skip(1) {
                self.bacteria.push(header.to_string());
            }
        }

        let size = self.size();

        // The expectation is that the matrix provided will always be
        // an NxN matrix.
        // Using a 1D vector for performance improvements.
        self.orig_graph.resize(size * size, 0.0);

        self.output.resize(size, 0.0);

        for (row, result) in reader.records().enumerate() {
            let row_value = result.expect("Unable to read CSV row");

            // Insert the row-column weights into the original graph.
            // Note: We skip the first column in each row
            // as that is just the bacteria name.
            for column in 1..row_value.len() {
                let weight = row_value[column]
                    .parse::<f32>()
                    .expect("Unable to parse string into float");
                // Change the column to accurately place it in the vector
                let column = column - 1;

                if row != column && weight != 0.0 {
                    if weight > 0.0 {
                        self.orig_graph
                            .insert(vec_offset(row, column, size), weight);
                        self.orig_graph
                            .insert(vec_offset(row + 1, column + 1, size), weight);
                    // self.orig_graph[(row + 1) * size + real_column] = 0.0;
                    // self.orig_graph[row * size + (real_column + 1)] = 0.0;
                    } else if weight < 0.0 {
                        self.orig_graph
                            .insert(vec_offset(row + 1, column, size), weight);
                        self.orig_graph
                            .insert(vec_offset(row, column + 1, size), weight);
                        // self.orig_graph[row * size + real_column] = 0.0;
                        // self.orig_graph[(row + 1) * size + (real_column + 1)] = 0.0
                    } // else {
                      //     if real_column + 1 > size {
                      //         continue;
                      //     }
                      //     self.orig_graph[row * size + real_column] = 0.0;
                      //     self.orig_graph[(row + 1) * size + (real_column + 1)] = 0.0;
                      //     self.orig_graph[(row + 1) * size + real_column] = 0.0;
                      //     self.orig_graph[row * size + (real_column + 1)] = 0.0;
                      // }
                } else {
                    self.orig_graph.insert(vec_offset(row, column, size), 1.0);

                    if vec_offset(row + 1, column + 1, size) < size * size {
                        self.orig_graph
                            .insert(vec_offset(row + 1, column + 1, size), 1.0);
                    }
                    // self.orig_graph[(row + 1) * size + real_column] = 0.0;
                    // self.orig_graph[row * size + (real_column + 1)] = 0.0;
                }
            }
        }

        Ok(())
    }

    /// Run the ATria algorithm over the input data.
    fn run(&mut self) -> Result {
        info!("I am running ATria");

        let size = self.size();
        // For parallelization, we divide our work load into smaller, workable partitions.
        let partitions = size / PARTITION_SIZE + 1;

        // ATria iterates over all rows and columns for
        // as large as the graph is.
        for _ in 0..size {
            let mut max_node = 0;
            let mut max_pay = -1.0;
            // Adjacency graph scope. Drops adjacency allocations once computations are completed.
            {
                // Copy the original graph to work on.
                let mut adj_graph = self.orig_graph.clone();
                let adj_graph_guard = Arc::new(Mutex::new(&mut adj_graph));

                // Running a modified Floyd-Warshall Algorithm over each partition.
                (0..partitions).into_par_iter().for_each(|partition| {
                    let adj_graph = Arc::clone(&adj_graph_guard);
                    let mut adj_graph = adj_graph
                        .lock()
                        .expect("Could not lock mutex on adjacency graph");

                    let min = partition * PARTITION_SIZE;
                    let max = if min + PARTITION_SIZE < size {
                        min + PARTITION_SIZE
                    } else {
                        size
                    };

                    for k in min..max {
                        for i in min..max {
                            for j in min..max {
                                if i != j && j != k {
                                    let even_odd = i + j;
                                    let current = vec_offset(i, j, partition);
                                    let current_node = adj_graph[current];
                                    let comperitor = adj_graph[vec_offset(i, k, partition)]
                                        * adj_graph[vec_offset(k, j, partition)];

                                    if (even_odd % 2 == 0) && (current_node < comperitor)
                                        || (even_odd % 2 == 1) && (current_node > comperitor)
                                    {
                                        adj_graph.insert(current, comperitor);
                                    }
                                }
                            }
                        }
                    }
                });

                let mut output_pay = vec![0.0; size];

                for i in 0..size {
                    let mut pay = 0.0;
                    for j in 0..size {
                        pay += adj_graph[vec_offset(i, j, size)];
                    }
                    pay = pay - 1.0;
                    output_pay.insert(i, pay);
                }

                for i in 0..size {
                    if fabsf(output_pay[i]) > max_pay {
                        max_node = i;
                        max_pay = fabsf(output_pay[i]);
                    }
                }

                self.output.insert(max_node, output_pay[max_node]);
            }

            if max_pay == 0.0 {
                break;
            }

            // Non-GPU Triad removal
            let orig_graph = Arc::new(Mutex::new(&mut self.orig_graph));

            (0..partitions).into_par_iter().for_each(|partition| {
                debug!("Working on partition {}", partition);
                let o_g = Arc::clone(&orig_graph);
                let mut o_g = o_g
                    .lock()
                    .expect("Unable to get a runtime lock on the input graph data");
                let size = o_g.len();

                let min = partition * PARTITION_SIZE;
                let max = if min + PARTITION_SIZE > size {
                    size
                } else {
                    min + PARTITION_SIZE
                };

                for i in min..max {
                    if (i / 2) != max_node
                        && (o_g[vec_offset(max_node, i, partition)] != 0.0
                            || o_g[vec_offset(max_node + 1, i, partition)] != 0.0)
                    {
                        for j in min..max {
                            if (j / 2) != max_node
                                && ((o_g[vec_offset(max_node, j, partition)] != 0.0
                                    || o_g[vec_offset(max_node + 1, j, partition)] != 0.0)
                                    && o_g[vec_offset(i, j, partition)] != 0.0)
                            {
                                o_g[vec_offset(i, j, partition)] = f32::MAX;
                                o_g[vec_offset(j, i, partition)] = f32::MAX;
                            }
                        }

                        if o_g[vec_offset(max_node, i, partition)] != 0.0 {
                            o_g[vec_offset(max_node, i, partition)] = f32::MAX;
                            o_g[vec_offset(i, max_node, partition)] = f32::MAX;
                        }

                        if o_g[vec_offset(max_node + 1, i, partition)] != 0.0 {
                            o_g[vec_offset(max_node + 1, i, partition)] = f32::MAX;
                            o_g[vec_offset(i, max_node + 1, partition)] = f32::MAX;
                        }
                    }
                }
            });

            //Sweep through row == column values.
            for i in 0..size {
                if self.orig_graph[vec_offset(i, i, size)] == f32::MAX {
                    self.orig_graph[vec_offset(i, i, size)] = 0.0;
                }
            }
        }

        Ok(())
    }

    /// Write the results of the ATria calulations to a NOA file.
    fn output(&mut self, file_path: String) -> Result {
        let mut output_file = File::create(file_path).expect("Unable to open output file location");

        for i in (0..self.size()).rev() {
            for j in 0..i {
                if fabsf(self.output[j]) < fabsf(self.output[j + 1]) {
                    self.output.swap(j, j + 1);
                    self.bacteria.swap(j, j + 1);
                }
            }
        }

        write!(output_file, "Name\tCentrality\tRank\n")
            .expect("Unable to write headers to output file");

        let mut min = 0.0;
        let mut max = 0.0;

        for i in 0..self.size() - 1 {
            self.output[i] = fabsf(self.output[i]);

            if self.output[i] > max {
                max = self.output[i];
            } else if self.output[i] < min {
                min = self.output[i];
            }

            write!(
                output_file,
                "{}\t{}\t\t{}\n",
                self.bacteria[i],
                self.output[i],
                self.size() - i
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

        let size = 4;

        plugin.orig_graph = vec![
            0.0,
            f32::INFINITY,
            -2.0,
            f32::INFINITY,
            4.0,
            0.0,
            3.0,
            f32::INFINITY,
            f32::INFINITY,
            f32::INFINITY,
            0.0,
            2.0,
            f32::INFINITY,
            -1.0,
            f32::INFINITY,
            0.0,
        ];

        plugin.bacteria = vec![
            String::from("Test Bac 1"),
            String::from("Test Bac 2"),
            String::from("Test Bac 3"),
            String::from("Test Bac 4"),
        ];

        plugin.output.resize(size, 0.0);

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
