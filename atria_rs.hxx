// Copyright (C) Joseph R. Quinn
// SPDX-License-Identifier: MIT

#ifndef _ATRIAPLUGIN_H
#define _ATRIAPLUGIN_H

#include <string>

/*
#include "Plugin.h"
#include "PluginProxy.h"
*/

class ATriaPlugin/*: public Plugin */ {
public:
    ATriaPlugin();
    ~ATriaPlugin();
    void input(std::string file);
    void run();
    void output(std::string file);
private:
    void *internal;
};

#endif
