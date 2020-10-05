// Copyright (C) 2020 Joseph R. Quinn
// SPDX-License-Identifier: MIT

use std::fs::File;
use std::io::prelude::*;

use atria_rs::ATria;
use pluma_plugin_trait::PluMAPlugin;

#[test]
fn it_works() {
    let mut plugin = ATria::default();

    plugin.input("./tests/corrP.never.csv")
        .expect("Failed to read CSV file...");
    plugin.run()
        .expect("Failed to run ATria...");
    plugin.output("./tests/corrP.never.noa")
        .expect("Failed to write NOA file...");

    let mut expected = File::open("./tests/corrP.never.noa.expected")
        .expect("Failed to open expected output file...");
    let mut expected_content = String::new();
    expected.read_to_string(&mut expected_content)
        .expect("Failed to read expected output file.");

    let mut actual = File::open("./tests/corrP.never.noa")
        .expect("Failed to open generated output file...");
    let mut actual_content = String::new();
    actual.read_to_string(&mut actual_content)
        .expect("Failed to read actual output to file.");

    assert_eq!(expected_content, actual_content);
}
