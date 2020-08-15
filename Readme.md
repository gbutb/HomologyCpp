# HomologyCpp
[![Actions Status](https://github.com/gbutb/HomologyCpp/workflows/GTest/badge.svg)](https://github.com/gbutb/HomologyCpp/actions)

HomologyCpp is a header-only library for computing simplicial homology of meshes.

## Usage:
Given VCGLIB's mesh: `mesh`, initialize its chain complex using:
```cpp
homology::ChainComplex<MeshType, FaceType, EdgeType, VertexType> complex(mesh);
```
To compute `i`-th betti number, use:
```cpp
complex.getBetti(i);
```

## Building examples:
In order to build examples, run the following commands from the root directory of this repository:
```zsh
> mkdir build && cd build
> cmake -DBUILD_TESTS=OFF -DBUILD_EXAMPLES=ON ..
> make -j
```
