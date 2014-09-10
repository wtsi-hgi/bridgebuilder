#!/bin/sh

if [ ! -d ./gnulib ] ; then 
    ./bootstrap || (echo "bootstrap failed" && rm -rf ./gnulib && exit 1)
else
    autoreconf --verbose --install || (echo "autoreconf failed" && exit 1)
fi

