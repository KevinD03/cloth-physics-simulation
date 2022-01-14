#include <assemble_forces.h>
#include <iostream>

//Input:
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
//  f - the vector 3xn vector of forces acting on each node of the mass-spring system

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0,
                     double mu, double lambda) { 

    f.resize(q.rows());
    f.setZero();

    for (int i = 0; i < F.rows(); i++) {

        Eigen::Vector9d dV;
        Eigen::RowVector3i element = F.row(i);
        Eigen::Matrix<double, 1, 9> tmp_row;
        tmp_row = dX.row(i);
        // 1*9 to a 3*3
        Eigen::Map<const Eigen::Matrix3d> tmp_dX(tmp_row.data());

        dV_membrane_corotational_dq(dV, q, tmp_dX, V, element, a0(i), mu, lambda);

        // set up for each vertices, x0,x1,x2,x3
        // 
        // x0  0   ......
        // y0  0   ......
        // z0  0   ......
        //  0  x1  ......
        //  0  y1  ......
        //  0  z1  ......

        // setup x0
        f(3 * F(i, 0)) -= dV(0);
        f(3 * F(i, 0) + 1) -= dV(1);
        f(3 * F(i, 0) + 2) -= dV(2);

        // setup x1
        f(3 * F(i, 1)) -= dV(3);
        f(3 * F(i, 1) + 1) -= dV(4);
        f(3 * F(i, 1) + 2) -= dV(5);

        // set up x2
        f(3 * F(i, 2)) -= dV(6);
        f(3 * F(i, 2) + 1) -= dV(7);
        f(3 * F(i, 2) + 2) -= dV(8);
    }
};
