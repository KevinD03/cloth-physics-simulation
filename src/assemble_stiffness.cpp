#include <assemble_stiffness.h>
#include <iostream>

//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  dX - an mx9 matrix which stores the flattened dphi/dX matrices for each tetrahedron.
//       Convert this values back to 3x3 matrices using the following code (NOTE YOU MUST USE THE TEMPORARY VARIABLE tmp_row):
//       Eigen::Matrix<double, 1,9> tmp_row
//       tmp_row = dX.row(ei); //ei is the triangle index. 
//       Eigen::Map<const Eigen::Matrix3d>(tmp_row.data())
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  F - the mx3 triangle connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  a0 - the mx1 vector of undeformed triangle areas
//  mu,lambda - material parameters for the cloth material model
//Output:
//  K- the 3n by 3n sparse stiffness matrix. 

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0, 
                     double mu, double lambda) { 

    K.resize(q.rows(), q.rows());
    K.setZero();
    typedef Eigen::Triplet<double> Triplet;
    std::vector<Triplet> coordinate_triplets;
    coordinate_triplets.reserve(q.rows() * q.rows());

    for (int tri_i = 0; tri_i < F.rows(); tri_i++) {
        Eigen::Matrix99d d2V;
        Eigen::RowVector3i element = F.row(tri_i);
        Eigen::Matrix<double, 1, 9> tmp_row;
        tmp_row = dX.row(tri_i);
        // 1*9 to a 3*3
        Eigen::Map<const Eigen::Matrix3d> tmp_dX(tmp_row.data());

        d2V_membrane_corotational_dq2(d2V, q, tmp_dX, V, element, a0(tri_i), mu, lambda);

        // describe the relation between a selected vertice and other 3 vertices
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {

                coordinate_triplets.push_back(Triplet(3 * F(tri_i, i), 3 * F(tri_i, j), -d2V(3 * i, 3 * j)));
                coordinate_triplets.push_back(Triplet(3 * F(tri_i, i), 3 * F(tri_i, j) + 1, -d2V(3 * i, 3 * j + 1)));
                coordinate_triplets.push_back(Triplet(3 * F(tri_i, i), 3 * F(tri_i, j) + 2, -d2V(3 * i, 3 * j + 2)));

                coordinate_triplets.push_back(Triplet(3 * F(tri_i, i) + 1, 3 * F(tri_i, j), -d2V(3 * i + 1, 3 * j)));
                coordinate_triplets.push_back(Triplet(3 * F(tri_i, i) + 1, 3 * F(tri_i, j) + 1, -d2V(3 * i + 1, 3 * j + 1)));
                coordinate_triplets.push_back(Triplet(3 * F(tri_i, i) + 1, 3 * F(tri_i, j) + 2, -d2V(3 * i + 1, 3 * j + 2)));

                coordinate_triplets.push_back(Triplet(3 * F(tri_i, i) + 2, 3 * F(tri_i, j), -d2V(3 * i + 2, 3 * j)));
                coordinate_triplets.push_back(Triplet(3 * F(tri_i, i) + 2, 3 * F(tri_i, j) + 1, -d2V(3 * i + 2, 3 * j + 1)));
                coordinate_triplets.push_back(Triplet(3 * F(tri_i, i) + 2, 3 * F(tri_i, j) + 2, -d2V(3 * i + 2, 3 * j + 2)));

            }
        }
    }
    K.setFromTriplets(coordinate_triplets.begin(), coordinate_triplets.end());
};
