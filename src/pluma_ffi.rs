// Copyright (C) Joseph R. Quinn
// SPDX-License-Identifier: MIT
//
// PluMA FFI exports for the ATria plugin
// These functions are called by the PluMA C++ runtime to execute the plugin

use crate::ATriaPlugin;
use pluma_plugin_trait::PluMAPlugin;
use std::ffi::CStr;
use std::os::raw::c_char;

/// Create a new ATria plugin instance
#[unsafe(no_mangle)]
pub extern "C" fn plugin_create() -> *mut std::ffi::c_void {
    let plugin = Box::new(ATriaPlugin::default());
    Box::into_raw(plugin) as *mut std::ffi::c_void
}

/// Destroy an ATria plugin instance
#[unsafe(no_mangle)]
pub extern "C" fn plugin_destroy(ptr: *mut std::ffi::c_void) {
    if !ptr.is_null() {
        unsafe {
            let _ = Box::from_raw(ptr as *mut ATriaPlugin);
        }
    }
}

/// Execute the input phase
#[unsafe(no_mangle)]
pub extern "C" fn plugin_input(ptr: *mut std::ffi::c_void, filename: *const c_char) {
    if ptr.is_null() || filename.is_null() {
        eprintln!("[ATria] Error: null pointer in plugin_input");
        return;
    }
    unsafe {
        let plugin = &mut *(ptr as *mut ATriaPlugin);
        let filename_str = match CStr::from_ptr(filename).to_str() {
            Ok(s) => s.to_string(),
            Err(e) => {
                eprintln!("[ATria] Error converting filename: {}", e);
                return;
            }
        };
        if let Err(e) = plugin.input(filename_str) {
            eprintln!("[ATria] Error in input: {}", e);
        }
    }
}

/// Execute the run phase
#[unsafe(no_mangle)]
pub extern "C" fn plugin_run(ptr: *mut std::ffi::c_void) {
    if ptr.is_null() {
        eprintln!("[ATria] Error: null pointer in plugin_run");
        return;
    }
    unsafe {
        let plugin = &mut *(ptr as *mut ATriaPlugin);
        if let Err(e) = plugin.run() {
            eprintln!("[ATria] Error in run: {}", e);
        }
    }
}

/// Execute the output phase
#[unsafe(no_mangle)]
pub extern "C" fn plugin_output(ptr: *mut std::ffi::c_void, filename: *const c_char) {
    if ptr.is_null() || filename.is_null() {
        eprintln!("[ATria] Error: null pointer in plugin_output");
        return;
    }
    unsafe {
        let plugin = &mut *(ptr as *mut ATriaPlugin);
        let filename_str = match CStr::from_ptr(filename).to_str() {
            Ok(s) => s.to_string(),
            Err(e) => {
                eprintln!("[ATria] Error converting filename: {}", e);
                return;
            }
        };
        if let Err(e) = plugin.output(filename_str) {
            eprintln!("[ATria] Error in output: {}", e);
        }
    }
}
