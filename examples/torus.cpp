/**
 * torus.cpp
 * ---------
 * 
 * Computes betti numbers
 * of torus.
 */

#include <stdio.h>
#include <stdlib.h>

#include <homology/mesh.hpp>
#include <homology/ChainComplex.hpp>
using homology::ChainComplex;

int main(int argc, char** argv) {
    if (argc != 3) {
        printf("Usage: ./torus [horizontal radius] [vertical radius]\n");
        return EXIT_FAILURE;
    }

    // Parse args
    float r_h = atof(argv[1]);
    float r_v = atof(argv[2]);

    // Initialize torus
    Mesh mesh;
    vcg::tri::Torus(mesh, r_h, r_v);

    // Initialize chain complex
    ChainComplex<Mesh, Face, Edge, Vertex> complex(mesh);

    // Compute and print betti numbers
    for (int i = 0; i < 4; ++i)
        printf("betti%d = %d\n", i, complex.getBetti(i));

    // Print Euler characteristic
    printf("chiE = %d\n", complex.getEulerCharacteristic());

    return EXIT_SUCCESS;
}
