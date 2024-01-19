#pragma once

#include <vector>

#ifndef FTETWLD_EXPORT
#  if defined(_MSC_VER) || defined(__CYGWIN__)
#    ifdef FTetWild_EXPORT
#      define FTETWLD_EXPORT __declspec(dllexport)
#    else
#      define FTETWLD_EXPORT __declspec(dllimport)
#    endif
#  elif defined(__clang__) || defined(__GNUC__)
#    define FTETWLD_EXPORT __attribute__((visibility("default")))
#  endif
#endif


FTETWLD_EXPORT int runFTetWild(std::vector<std::vector<int>> &faces,
                       std::vector<std::vector<float>> &verts,
                       int argc, char **argv);
