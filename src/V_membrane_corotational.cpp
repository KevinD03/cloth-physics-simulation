#include <V_membrane_corotational.h>

//  q - generalized coordinates for the FEM system
//  dX - the 3x3 matrix containing dphi/dX
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the vertex indices of this triangle
//  area - the area of this triangle
//  mu,lambda - material parameters for the cloth material model
//Output:
//  energy- the per-triangle potential energy (the linear model described in the README).


//Allowed to use libigl SVD or Eigen SVD for this part
void V_membrane_corotational(double &energy, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {

    // see lecture 6 slide 73 for reference
    Eigen::Vector3d x0 = q.segment<3>(3 * element(0));
    Eigen::Vector3d x1 = q.segment<3>(3 * element(1));
    Eigen::Vector3d x2 = q.segment<3>(3 * element(2));
  
    Eigen::Vector3d X0 = V.row(element(0));
    Eigen::Vector3d X1 = V.row(element(1));
    Eigen::Vector3d X2 = V.row(element(2));

    Eigen::Vector3d n = (x1 - x0).cross(x2 - x0).normalized();
    Eigen::Vector3d N = (X1 - X0).cross(X2 - X0).normalized();

    // see slide 68 for the composition of F matrix
    Eigen::Matrix34d x;
    x.setZero();
    x.block(0, 0, 3, 1) = x0;
    x.block(0, 1, 3, 1) = x1;
    x.block(0, 2, 3, 1) = x2;
    x.block(0, 3, 3, 1) = n;

    Eigen::Matrix43d X;
    X.setZero();
    X.block(0, 0, 3, 3) = dX;
    X.block(3, 0, 1, 3) = N.transpose();


    Eigen::Matrix3d F = x * X;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Vector3d S = svd.singularValues();

    // see formula in readMe < Linear Elasticity without the Pesky Rotations >
    double psi = mu * (std::pow(S[0] - 1.0, 2) + std::pow(S[1] - 1.0, 2) + std::pow(S[2] - 1.0, 2)) +
        0.5 * lambda * std::pow(S[0] + S[1] + S[2] - 3.0, 2);
    
    energy = area * psi;

}
