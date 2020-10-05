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
use std::marker::PhantomData;

use csv;
use libm::fabsf;
use pluma_plugin_trait::PluMAPlugin;
use generic_floyd_warshall::floyd_warshall;

pub type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

pub type Matrix = Vec<Vec<f32>>;

#[derive(Debug)]
pub struct ATria<'a> {
    pub bacteria: Vec<String>,
    pub graph_size: usize,
    pub orig_graph: Matrix,
    pub output: Vec<f32>,
    pub output_pay: Vec<f32>,
    phantom: PhantomData<&'a str>
}

impl<'a> Default for ATria<'a> {
    fn default() -> Self {
        ATria {
            graph_size: 0,
            orig_graph: vec![vec![]],
            output: vec![],
            bacteria: vec![],
            output_pay: vec![],
            phantom: PhantomData
        }
    }
}

impl<'a> PluMAPlugin<'a> for ATria<'a> {
    /// Create a NxN adjacency matrix from the input CSV file.
    /// Space complextiy is expected to be O(n^2).
    fn input(&mut self, file_path: &'a str) -> crate::Result<()> {

        let mut reader = csv::Reader::from_path(file_path)?;

        // Our bacteria should be initialized as a Vector of n elements.
        {
            let headers = reader.headers()?.clone();
            for header in headers.iter().skip(1) {
                self.bacteria.push(header.to_string());
            }
        }
        
        self.graph_size = self.bacteria.len();

        // The expectation is that the matrix provided
        // will always be an NxN matrix.
        self.orig_graph = vec![vec![0.0; self.graph_size]; self.graph_size];

        for (row, result) in reader.records().skip(1).enumerate() {
            let row_value = result?;

            for column in 1..row {
                let weight = row_value[column].parse::<f32>()?;

                if row != column {
                    if weight > 0.0 {
                        self.orig_graph[row][column] = weight;
                        self.orig_graph[row + 1][column + 1] = weight;
                    } else if weight < 0.0 {
                        self.orig_graph[row + 1][column] = weight;
                        self.orig_graph[row][column + 1] = weight;
                    } else {
                        self.orig_graph[row][column] = 0.0;
                        self.orig_graph[row + 1][column + 1] = 0.0;
                    }
                } else {
                    self.orig_graph[row][column] = 1.0;
                    self.orig_graph[row + 1][column + 1] = 1.0;
                }
            }
        }

        self.output.resize(self.graph_size, 0.0);
        self.output_pay.resize(self.graph_size, 0.0);

        Ok(())
    }

    /// Run the ATria algorithm over the input data.
    fn run(&mut self) -> crate::Result<()> {
        println!("I am running ATria");

        for _a in 0..self.graph_size {
            let output_graph = floyd_warshall(&mut self.orig_graph.clone(), self.graph_size);

            for i in 0..self.graph_size {
                let mut pay = 0.0;
                for j in 0..self.graph_size {
                    pay = pay + output_graph[i][j];
                }
                pay = pay - 1.0;
                self.output_pay[i] = pay;
            }

            let mut max_node = 0;
            let mut max_pay = -1.0;

            for i in 0..self.graph_size {
                if fabsf(self.output_pay[i]) > max_pay {
                    max_node = i;
                    max_pay = fabsf(self.output_pay[i]);
                }
            }

            self.output[max_node] = self.output_pay[max_node];

            if max_pay == 0.0 {
                break;
            }

            // Non-GPU Triad Removal
            for i in 0..self.graph_size {
                if (i/2 != max_node) &&
                    ((self.orig_graph[max_node][i] != 0.0) ||
                    (self.orig_graph[max_node + 1][i] != 0.0))
                {
                    for j in i + 1..self.graph_size {
                        if (j/2 != max_node) && 
                            ((self.orig_graph[max_node][j] != 0.0 || 
                            self.orig_graph[max_node + 1][j] != 0.0) && 
                            (self.orig_graph[i][j] != 0.0))
                        {
                            self.orig_graph[i][j] = 2.0;
                            self.orig_graph[j][i] = 2.0;
                        }
                    }

                    if self.orig_graph[max_node][i] != 0.0 {
                        self.orig_graph[max_node][i] = 2.0;
                        self.orig_graph[i][max_node] = 2.0;
                    }

                    if self.orig_graph[max_node + 1][i] != 0.0 {
                        self.orig_graph[max_node + 1][i] = 2.0;
                        self.orig_graph[i][max_node + 1] = 2.0;
                    }
                }
            }

            // Sweep through row == column values
            for i in 0..self.graph_size {
                if self.orig_graph[i][i] == 2.0 {
                    self.orig_graph[i][i] = 0.0;
                }
            }
        }

        Ok(())
    }

    /// Write the results of the ATria calulations to a NOA file.
    fn output(&mut self, file_path: &'a str) -> crate::Result<()> {
        let mut output_file = File::create(file_path)?;

        for i in (self.graph_size - 1)..0 {
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

        for i in 0..self.graph_size {
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
                self.graph_size-i
            )?;
        }

        Ok(())
    }
}