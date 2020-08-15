/**
 * main.cpp
 * --------
 */

#include <stdio.h>
#include <stdlib.h>

#include <homology/mesh.hpp>
#include <homology/ChainComplex.hpp>
using homology::ChainComplex;

#include <iostream>

int main() {
    Mesh mesh;
    vcg::tri::Torus(mesh, 10, 3);
    ChainComplex<Mesh, Face, Edge, Vertex> complex(mesh);
    Eigen::MatrixXf mat = complex.getd(1)*complex.getd(2);
    assert(mat.isZero(0));

    for (int i = 0; i < 3; ++i) {
        std::cout << "getBetti(" << i << ") = " << complex.getBetti(i) << std::endl;
    }

    std::cout << "chiE = " << complex.getEulerCharacteristic() << std::endl;

    return EXIT_SUCCESS;
}
