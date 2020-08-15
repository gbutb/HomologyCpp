/**
 * main_test.cpp
 * -------------
 * Contains unit tests.
 */

#include "gtest/gtest.h"

#include <homology/mesh.hpp>
#include <homology/ChainComplex.hpp>
using homology::ChainComplex;

TEST(testHomology, torus) {
    Mesh mesh;
    vcg::tri::Torus(mesh, 10, 3);
    ChainComplex<Mesh, Face, Edge, Vertex> complex(mesh);
    Eigen::MatrixXf mat = complex.getd(1)*complex.getd(2);
    EXPECT_TRUE(mat.isZero(0));

    int expected[] = {1, 2, 1, 0};
    for (int i = 0; i < 4; ++i)
        EXPECT_EQ(expected[i], complex.getBetti(i));

    EXPECT_EQ(0, complex.getEulerCharacteristic());
}

TEST(testHomology, tetrahedron) {
    Mesh mesh;
    vcg::tri::Tetrahedron(mesh);
    ChainComplex<Mesh, Face, Edge, Vertex> complex(mesh);
    Eigen::MatrixXf mat = complex.getd(1)*complex.getd(2);
    EXPECT_TRUE(mat.isZero(0));

    int expected[] = {1, 0, 1, 0};
    for (int i = 0; i < 4; ++i)
        EXPECT_EQ(expected[i], complex.getBetti(i));

    EXPECT_EQ(2, complex.getEulerCharacteristic());
}

TEST(testChainComplex, repeatedCalls) {
    int num_repeats = 5;
    Mesh mesh;
    vcg::tri::Tetrahedron(mesh);
    ChainComplex<Mesh, Face, Edge, Vertex> complex(mesh);

    for (int repeat = 0; repeat < num_repeats; ++repeat) {
        int expected[] = {1, 0, 1, 0};
        for (int i = 0; i < 4; ++i)
            EXPECT_EQ(expected[i], complex.getBetti(i));

        EXPECT_EQ(2, complex.getEulerCharacteristic());
    }
}

int main(int argc, char**argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
