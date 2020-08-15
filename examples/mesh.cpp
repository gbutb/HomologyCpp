/**
 * mesh.cpp
 * --------
 * 
 * Computes betti numbers
 * of a mesh.
 */

#include <stdio.h>
#include <stdlib.h>

#include <homology/mesh.hpp>
#include <homology/ChainComplex.hpp>
using homology::ChainComplex;

#include <wrap/io_trimesh/import_obj.h>

int main(int argc, char** argv) {
    if (argc != 2) {
        printf("Usage: ./mesh [path_to_obj]");
        return EXIT_FAILURE;
    }
    // Parse args
    char* filename = argv[1];

    // Load mesh
    Mesh mesh;
    vcg::tri::io::ImporterOBJ<Mesh>::Info info;
    vcg::tri::io::ImporterOBJ<Mesh>::Open(mesh, filename, info);
    printf("Loaded mesh with: %d vertices and %d faces\n",
        info.numVertices, info.numFaces);

    // Initialize chain complex
    ChainComplex<Mesh, Face, Edge, Vertex> complex(mesh);

    // Compute and print betti numbers
    for (int i = 0; i < 4; ++i)
        printf("betti%d = %d\n", i, complex.getBetti(i));

    // Print Euler characteristic
    printf("chiE = %d\n", complex.getEulerCharacteristic());

    return EXIT_SUCCESS;
}
