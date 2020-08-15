/**
 * main_test.cpp
 * -------------
 * Contains unit tests.
 */

#include "gtest/gtest.h"
#include <homology/mesh.hpp>
#include <homology/ChainComplex.hpp>


TEST(testHomology, torus) {
    Mesh mesh;
    vcg::tri::Torus(mesh, 10, 3);
    ChainComplex<Mesh, Face, Edge, Vertex> complex(mesh);
    Eigen::MatrixXf mat = complex.getd(1)*complex.getd(2);
    EXPECT_TRUE(mat.isZero(0));

    EXPECT_EQ(1, complex.getBetti(0));
    EXPECT_EQ(2, complex.getBetti(1));
    EXPECT_EQ(1, complex.getBetti(2));
    EXPECT_EQ(0, complex.getBetti(3));

    EXPECT_EQ(0, complex.getEulerCharacteristic());
}

TEST(testHomology, tetrahedron) {
    Mesh mesh;
    vcg::tri::Tetrahedron(mesh);
    ChainComplex<Mesh, Face, Edge, Vertex> complex(mesh);
    Eigen::MatrixXf mat = complex.getd(1)*complex.getd(2);
    EXPECT_TRUE(mat.isZero(0));

    EXPECT_EQ(1, complex.getBetti(0));
    EXPECT_EQ(0, complex.getBetti(1));
    EXPECT_EQ(1, complex.getBetti(2));
    EXPECT_EQ(0, complex.getBetti(3));

    EXPECT_EQ(2, complex.getEulerCharacteristic());
}

int main(int argc, char**argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
