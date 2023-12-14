# sv_merge
C++ implementation of rlorigro/hapslap


## Dependencies

- c++17
- gcc or clang
- cmake
- zlib1g-dev
- libbz2-dev
- lzma
- autoconf
- automake
- libssl-dev
- pkg-config
- libnghttp2-dev (For curl http2 support, **only required for Ubuntu22**)
- libcurl-dev (generally libcurl4-openssl-dev, currently **only required for MacOS**)

## Installation

First verify that the system level dependencies are installed and then run:
```
git clone https://github.com/rlorigro/sv_merge
cd sv_merge
mkdir build
cmake ..
make -j [n_threads]
```

(if any dependencies are not installed you should see an error in the `cmake` step)

