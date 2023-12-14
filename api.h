#pragma once

#include <vector>

#ifndef EXPORT
#  if defined(_MSC_VER) || defined(__CYGWIN__)
#    ifdef FTetWild_EXPORT
#      define EXPORT __declspec(dllexport)
#    else
#      define EXPORT __declspec(dllimport)
#    endif
#  elif defined(__clang__) || defined(__GNUC__)
#    define EXPORT __attribute__((visibility("default")))
#  endif
#endif


EXPORT int runFTetWild(int argc, char **argv);
