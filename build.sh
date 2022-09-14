#!/usr/bin/env bash
build_dir="build"

cmake "$@" -B "${build_dir}"
cmake --build "${build_dir}"
