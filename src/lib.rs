// Copyright (C) 2020 Joseph R. Quinn
// SPDX-License-Identifier: MIT

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

use csv;
use libm::fabsf;
use pluma_plugin_trait::PluMAPlugin;
use rayon::prelude::*;

/// Crate Result type.
pub type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

mod dijkstras;

/// Matrix type using built-in Rust Vectors to create a dynamically sized 2D array
pub type Matrix = Vec<Vec<f32>>;

#[repr(C)]
#[derive(Debug)]
pub struct ATria {
    /// Vector of the bacteria types in the CSV file.
    pub bacteria: Vec<String>,
    /// The original matrix being worked on by the ATria algorithm.
    pub orig_graph: Matrix,
    ///
    pub output: Vec<f32>,
    ///
    pub output_pay: Vec<f32>,
}

impl Default for ATria {
    fn default() -> Self {
        ATria {
            bacteria: vec![],
            orig_graph: vec![vec![]],
            output: vec![],
            output_pay: vec![],
        }
    }
}

impl ATria {
    /// Get the internal size of the original graph.
    pub fn len(&self) -> usize {
        self.bacteria.len() + 1
    }
}

impl PluMAPlugin for ATria {
    /// Create a NxN adjacency matrix from the input CSV file.
    /// Space complextiy is expected to be O(n^2).
    fn input<'a>(&mut self, file_path: &'a str) -> crate::Result<()> {

        let mut reader = csv::Reader::from_path(file_path)?;

        {
            let headers = reader.headers()?.clone();
            for header in headers.iter().skip(1) {
                self.bacteria.push(header.to_string());
            }
        }

        // The expectation is that the matrix provided will always be
        // an NxN matrix.
        self.orig_graph = vec![vec![0.0; self.len()]; self.len()];

        self.output.resize(self.len(), 0.0);
        self.output_pay.resize(self.len(), 0.0);

        for (row, result) in reader.records().enumerate() {
            let row_value = result?;

            for column in 1..row_value.len() {
                let column = column - 1;
                let weight = row_value[column + 1].parse::<f32>().unwrap();

                if row != column {
                    if weight > 0.0 {
                        self.orig_graph[row][column] = weight;
                        self.orig_graph[row + 1][column + 1] = weight;
                        self.orig_graph[row + 1][column] = 0.0;
                        self.orig_graph[row][column + 1] = 0.0;
                    } else if weight < 0.0 {
                        self.orig_graph[row + 1][column] = weight;
                        self.orig_graph[row][column + 1] = weight;
                        self.orig_graph[row][column] = 0.0;
                        self.orig_graph[row + 1][column + 1] = 0.0
                    } else {
                        self.orig_graph[row][column] = 0.0;
                        self.orig_graph[row + 1][column + 1] = 0.0;
                        self.orig_graph[row + 1][column] = 0.0;
                        self.orig_graph[row][column + 1] = 0.0;
                    }
                } else {
                    self.orig_graph[row][column] = 1.0;
                    self.orig_graph[row + 1][column + 1] = 1.0;
                    self.orig_graph[row + 1][column] = 0.0;
                    self.orig_graph[row][column + 1] = 0.0;
                }
            }
        }

        Ok(())
    }

    /// Run the ATria algorithm over the input data.
    fn run(&mut self) -> crate::Result<()> {
        println!("I am running ATria");

        let mut output_graph = self.orig_graph.clone();

        // ATria iterates over all rows and columns for
        // as large as the graph is.
        for _count in 0..self.len() {

            for vertex in 0..self.len() {
                let adjacency = dijkstras::dijkstras(
                    &self.orig_graph,
                    vertex
                );

                output_graph[vertex] = adjacency;

                for i in 0..self.len() {
                    let mut pay = 0.0;
                    for j in 0..self.len() {
                        pay = pay + self.output[j];
                    }
                    pay = pay - 1.0;
                    self.output_pay[i] = pay;
                }

                let mut max_node = 0;
                let mut max_pay = 0.0;

                for i in 0..self.len() {
                    if fabsf(self.output_pay[i]) > max_pay {
                        max_node = i;
                        max_pay = fabsf(self.output_pay[i]);
                    }
                }

                self.output[max_node] = self.output_pay[max_node];

                if max_pay == 0.0 {
                    break;
                }

                // Non-GPU Triad removal
                //
                // Divide and conquer by 25-unit chunks
                self.orig_graph.par_chunks_mut(25)
                    .enumerate()
                    .for_each(|(i, slice)| {
                        let step = i + 1;

                        // Offset the step by 25 to account
                        // for chunk size.
                        for j in 0..(step * 25) {
                            if (j/2) != max_node &&
                                (self.orig_graph[max_node][j] != 0.0 ||
                                self.orig_graph[max_node + 1][j] != 0.0) {
                                for k in j + 1..self.len() {
                                    if k / 2 != max_node &&
                                        ((self.orig_graph[max_node][k] != 0.0 ||
                                        self.orig_graph[max_node + 1][k] != 0.0) &&
                                        self.orig_graph[j][k] != 0.0) {
                                        self.orig_graph[j][k] = f32::MAX;
                                        self.orig_graph[k][j] = f32::MAX;
                                    }
                                }
                            }

                            if self.orig_graph[max_node][j] != 0.0 {
                                self.orig_graph[max_node][j] = f32::MAX;
                                self.orig_graph[j][max_node] = f32::MAX;
                            }

                            if self.orig_graph[max_node + 1][j] != 0.0 {
                                self.orig_graph[max_node + 1][j] = f32::MAX;
                                self.orig_graph[j][max_node + 1] = f32::MAX;
                            }
                        }
                    });

                // Sweep through row == column values
                for i in 0..self.len() {
                    if self.orig_graph[i][i] == f32::MAX {
                        self.orig_graph[i][i] = 0.0;
                    }
                }
            }
        }

        Ok(())
    }

    /// Write the results of the ATria calulations to a NOA file.
    fn output<'a>(&mut self, file_path: &'a str) -> crate::Result<()> {
        let mut output_file = File::create(file_path)?;

        for i in (self.len() - 1)..0 {
            for j in 0..i {
                if fabsf(self.output[j]) < fabsf(self.output[j + 1]) {
                    self.output.swap(j, j + 1);
                    self.bacteria.swap(j, j + 1);
                }
            }
        }

        write!(output_file, "Name\tCentrality\tRank\n")?;

        let mut min = 0.0;
        let mut max = 0.0;

        for i in 0..self.len() - 1 {
            self.output[i] = fabsf(self.output[i]);

            if self.output[i] > max {
                max = self.output[i];
            } else if self.output[i] < min {
                min = self.output[i];
            }

            write!(output_file,
                "{}\t{}\t\t{}\n",
                self.bacteria[i],
                self.output[i],
                self.len() - i
            )?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_can_run() {
        let mut plugin = ATria::default();

        plugin.orig_graph = vec![
            vec![0.0, f32::INFINITY, -2.0, f32::INFINITY],
            vec![4.0, 0.0, 3.0, f32::INFINITY],
            vec![f32::INFINITY, f32::INFINITY, 0.0, 2.0],
            vec![f32::INFINITY, -1.0, f32::INFINITY, 0.0]
        ];

        plugin.output.resize(plugin.orig_graph.len(), 0.0);
        plugin.output_pay.resize(plugin.orig_graph.len(), 0.0);

        assert!(plugin.run().is_ok());
    }
}
