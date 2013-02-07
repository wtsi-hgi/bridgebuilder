#!/bin/sh

if [[ ! -d gnulib ]]; then 
    ./bootstrap
fi

autoreconf -v --install || exit 1
