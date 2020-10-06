// Copyright (C) 2020 Joseph R. Quinn
// SPDX-License-Identifier: MIT

use atria_rs::ATria;
use pluma_plugin_trait::PluMAPlugin;

#[test]
fn it_should_load_bacteria() {
    let mut plugin = ATria::default();

    plugin.input("./tests/corrP.never.csv")
        .expect("Failed to read CSV file");
    
    assert_eq!(126, plugin.bacteria.len());
}