// Copyright (C) 2020 Joseph R. Quinn
// SPDX-License-Identifier: MIT

use atria_rs::ATria;
use pluma_plugin_trait::PluMAPlugin;
use csv;

#[test]
#[ignore]
fn it_should_correctly_load_matrix() {

    let mut reader = csv::Reader::from_path("./tests/corrP.never.csv")
        .expect("Failed to read input CSV file");

    let _headers = reader.headers()
        .expect("Failed to read headers");

    let mut plugin = ATria::default();

    plugin.input("./tests/corrP.never.csv")
        .expect("Failed to read CSV file");

    let size = plugin.orig_graph.len();

    for i in 0..size {
        let row = &plugin.orig_graph[i];
        assert_eq!(size, row.len());

        for j in 0..size {
            let column = row[j];
            assert!(column >= f32::NEG_INFINITY);
        }
    }

    // @TODO: Build out test to include reading the CSV file again
    // and iterating through it to check correctness.

    for (row, result) in reader.records().enumerate() {
        let row_value = result.expect("Failed to load row_value");

        for column in 1..row_value.len() {
            let weight = row_value[column].parse::<f32>()
                .expect("Failed to parse value as a float");

            assert_eq!(
                plugin.orig_graph[row][column],
                weight,
                "\nRow {} Column {}: Weight expected: {} - Weight actual: {}",
                row,
                column,
                plugin.orig_graph[row][column],
                weight
            );
        }
    }
}
