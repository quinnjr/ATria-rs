// Copyright (C) Joseph R. Quinn
// SPDX-License-Identifier: MIT

use cpp::cpp;

use crate::ATriaPlugin;

cpp!{{
    #include "atria_rs.hxx"

    ATriaPlugin::ATriaPlugin() {
        this->internal = rust!(ATriaPlugin_constructor [] -> *mut ATriaPlugin as "void *" {
            let b = Box::new(ATriaPlugin::default());
            Box::into_raw(b)
        });
    }

    ATriaPlugin::~ATriaPlugin() {
        rust!(ATriaPlugin_destructor [internal: *mut ATriaPlugin as "void *"] {
            let _b = unsafe {
                Box::from_raw(internal)
            };
        });
    }

    void ATriaPlugin::input(std::string file) {

    }

    void ATriaPlugin::run() {

    }

    void ATriaPlugin::output(std::string file) {

    }
}}
