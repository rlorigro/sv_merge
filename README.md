# Hapestry

For merging, annotating, and filtering population datasets, using graph alignment, multiobjective optimization, and machine learning.

# Installation

## Installing from source

You will need a C++20 compiler (for example [GCC 13](https://gcc.gnu.org)), [CMake](https://cmake.org/download/) version 3.22 or above,
and the following packages:
- `autoconf`
- `automake`
- `pkg-config`
- `zlib1g-dev`
- `libbz2-dev`
- `libssl-dev`
- `liblzma-dev`
- `libjansson-dev`

If you are on Ubuntu 22, you will also need:
- `libnghttp2-dev`, for curl HTTP2 support.

If you are on macOS, you will also need:
- `libcurl-dev` (generally `libcurl4-openssl-dev`).

Finally, you need to install [GraphAligner](https://github.com/maickrau/GraphAligner) v1.0.19 and add it to you `PATH`.

To install Hapestry:
```
git clone https://github.com/rlorigro/sv_merge.git \
&& cd sv_merge \
&& git checkout dev \
&& mkdir build \
&& cd build \
&& cmake .. \
&& make -j $N_THREADS
```

## Workflows and docker images

#### Hapestry merge WDL
https://dockstore.org/workflows/github.com/rlorigro/sv_merge/hapestry_merge:fc_fetching_tests?tab=info

#### Hapestry merge docker image
https://hub.docker.com/layers/fcunial/hapestry/merge/images/sha256-8fa681256cfe2f459b6d3fdce9e160f880ecc49f2b9cb39fa84699b49c890dd3

