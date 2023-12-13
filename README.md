# sv_merge
C++ implementation of rlorigro/hapslap


## Dependencies

- c++17
- gcc
- cmake
- zlib1g-dev
- libbz2-dev
- lzma
- autoconf
- automake
- libssl-dev
- pkg-config
- libnghttp2-dev (For Ubuntu22 curl http2 support)
- libcurl-dev (generally libcurl4-openssl-dev, currently only required for MacOS)
- libjansson4

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
