/**
 * main.cpp
 * --------
 */

#include <stdio.h>
#include <stdlib.h>

#include <homology/ChainComplex.hpp>

class Face;
class Vertex;
class Edge;

struct UsedTypes: public vcg::UsedTypes<
    vcg::Use<Vertex>::AsVertexType,
    vcg::Use<Face>::AsFaceType,
    vcg::Use<Edge>::AsEdgeType> {};

class Vertex : public vcg::Vertex<
    UsedTypes,
    vcg::vertex::Coord3f,
    vcg::vertex::Normal3f,
    vcg::vertex::Color4b,
    vcg::vertex::VFAdj,
    vcg::vertex::BitFlags> {};

class Face : public vcg::Face<
    UsedTypes,
    vcg::face::VertexRef,
    vcg::face::Normal3f,
    vcg::face::FFAdj,
    vcg::face::FEAdj,
    vcg::face::VFAdj,
    vcg::face::Mark,
    vcg::face::BitFlags> {};

class Edge : public vcg::Edge<
    UsedTypes,
    vcg::edge::VertexRef,
    vcg::edge::EEAdj,
    vcg::edge::EFAdj,
    vcg::edge::Mark,
    vcg::edge::BitFlags> {};

class Mesh : public vcg::tri::TriMesh<
    std::vector<Vertex>,
    std::vector<Face>,
    std::vector<Edge>> {};

#include <iostream>

int main() {
    Mesh mesh;
    vcg::tri::Tetrahedron(mesh);
    ChainComplex<Mesh, Face, Edge, Vertex> complex(mesh);

    std::cout << complex.getd(1) << std::endl;
    std::cout << std::endl;
    std::cout << complex.getd(2) << std::endl;

    Eigen::MatrixXf mat = complex.getd(1)*complex.getd(2);
    assert(mat.isZero(0));

    for (int i = 0; i < 3; ++i) {
        std::cout << "getBetti(" << i << ") " << complex.getBetti(i) << std::endl;
    }

    std::cout << "Euler Char. = " << complex.getEulerCharacteristic() << std::endl;

    return EXIT_SUCCESS;
}
