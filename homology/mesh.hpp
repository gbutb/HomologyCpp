/**
 * mesh.hpp
 * --------
 * 
 * Defines minimal Mesh
 * and its components s.t.
 * all requirements for ChainComplex
 * hold.
 */

#ifndef HOMOLOGY_MESH_HPP_
#define HOMOLOGY_MESH_HPP_

// VCG
#include <vcg/complex/complex.h>

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

#endif  // HOMOLOGY_MESH_HPP_