#!/bin/bash

julia "$@" --project -e 'using Pkg; Pkg.test()'
