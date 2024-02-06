# sv_merge
C++ implementation of rlorigro/hapslap



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

To install Hapestry, just do the following:
```
git clone https://github.com/rlorigro/sv_merge.git \
&& cd sv_merge \
&& git checkout dev \
&& mkdir build \
&& cd build \
&& cmake .. \
&& make -j $N_THREADS
```

## Docker

We also provide a Docker image at...







<!-- REMOVED BECAUSE GRAPHALIGNER API IS BAD
### GraphAligner

- protobuf-compiler
- libsparsehash-dev
- libsdsl-dev
- jemalloc (source installation required on Ubuntu 22.04)

```
cd ~/software/
git clone https://github.com/jemalloc/jemalloc.git
cd jemalloc/
./autogen.sh
make -j [n_threads]
sudo make install
```
-->
