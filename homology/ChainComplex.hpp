/**
 * ChainComplex.hpp
 * ----------------
 * 
 * Contains declaration of
 * the ChainComplex for a mesh.
 */

#ifndef HOMOLOGY_CHAINCOMPLEX_HPP_
#define HOMOLOGY_CHAINCOMPLEX_HPP_

// VCG
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/topology.h>

// Other
#include <wrap/io_trimesh/export_obj.h>
#include <eigenlib/Eigen/Core>

template<class MeshType, class FaceType, class EdgeType, class VertexType>
class ChainComplex {
 private:
    // Boundary maps
    Eigen::MatrixXf _d1;
    Eigen::MatrixXf _d2;

    // Properties of the mesh
    bool _orientable;
    int b0 = -1, b1 = -1, b2 = -1;
    int _num_faces;
    int _num_edges;
    int _num_vertices;

    /**
     * Computes boundary maps.
     * @param mesh input mesh
     * @param d1 output d1 boundary map.
     * @param d2 output d2 boundary map.
     */
    static void computeBoundaryMaps(
            MeshType& mesh,
            Eigen::MatrixXf& d1,
            Eigen::MatrixXf& d2) {
        // Initialize boundary maps
        d1 = Eigen::MatrixXf::Zero(mesh.vert.size(), mesh.edge.size());
        d2 = Eigen::MatrixXf::Zero(mesh.edge.size(), mesh.face.size());

        for (int face_i = 0; face_i < mesh.face.size(); ++face_i) {
            FaceType face = mesh.face.at(face_i);

            // Iterate over each edge of face_i
            for (int edge_i = 0; edge_i < 3; ++edge_i) {
                EdgeType* edge = face.FEp(edge_i);
                int edge_orientation = (edge->V(0) == face.V(0) ? 1 : -1);
                if (!edge->IsV()) {
                    // Update d1 boundary map
                    d1(edge->V(0) - &mesh.vert.at(0), edge - &mesh.edge.at(0)) += 1;
                    d1(edge->V(1) - &mesh.vert.at(0), edge - &mesh.edge.at(0)) -= 1;
                }
                d2(edge - &mesh.edge.at(0), face_i) += edge_orientation;
                edge->SetV();
            }
        }
    }

    /**
     * Prepares the mesh
     * @param mesh the mesh.
     * @param orientable writes true if the mesh is orientable,
     *      otherwise, writes false.
     */
    static void prepareMesh(MeshType& mesh, bool& orientable) {
        bool isOriented = false;
        vcg::tri::UpdateTopology<MeshType>::FaceFace(mesh);
        vcg::tri::Clean<MeshType>::OrientCoherentlyMesh(mesh, isOriented, orientable);
        vcg::tri::UpdateTopology<MeshType>::AllocateEdge(mesh);
    }

 public:
    /**
     * Initializes ChainComplex from mesh.
     * @param mesh A mesh from which the
     *      chain complex will be constructed.
     * @note requires FE adjacency.
     */
    ChainComplex(MeshType& mesh) {
        vcg::tri::RequireFEAdjacency(mesh);
        prepareMesh(mesh, _orientable);

        _num_faces = mesh.face.size();
        _num_edges = mesh.edge.size();
        _num_vertices = mesh.vert.size();


        computeBoundaryMaps(mesh, _d1, _d2);
    }

    /**
     * @param i index of the boundary map.
     * @return a matrix representation of
     *      the i-th boundary map.
     */
    Eigen::MatrixXf getd(int i) const {
        if (i == 1) return _d1;
        else if (i == 2) return _d2;
        return Eigen::MatrixXf::Zero(1, 1);
    }

    /**
     * @param i index of betti number.
     * @return i-th betti number. -1 if invalid index.
     */
    int getBetti(int i) const {
        switch (i) {
            case 0:
                break;
            case 1:
                break;
            case 2:
                break;
        }
        return -1;
    }

    /**
     * @return Euler Characteristic
     *      of the mesh.
     */
    int eulerCharacteristic() const;
};

#endif  // HOMOLOGY_CHAINCOMPLEX_HPP_
