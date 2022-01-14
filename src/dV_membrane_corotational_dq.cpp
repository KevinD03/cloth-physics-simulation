#include <dV_membrane_corotational_dq.h>
#include <iostream>

//  q - generalized coordinates for the FEM system
//  dX - the 3x3 matrix containing dphi/dX
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the vertex indices of this triangle
//  area - the area of this triangle
//  mu,lambda - material parameters for the cloth material model
//Output:
//  f - the per-triangle gradient of the membrane potential energy (the linear model described in the README).

void dV_membrane_corotational_dq(Eigen::Vector9d &dV, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {

    //Deformation Gradient
    Eigen::Matrix3d dx; //deformed tangent matrix 
    Eigen::Matrix3d U;
    Eigen::Vector3d S; 
    Eigen::Matrix3d W; 

    //TODO: SVD Here
    
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
    U = svd.matrixU();
    W = svd.matrixV();
    S = svd.singularValues();

    //Fix for inverted elements (thanks to Danny Kaufman)
    double det = S[0]*S[1];
    
     if(det <= -1e-10)
    {
        if(S[0] < 0) S[0] *= -1;
        if(S[1] < 0) S[1] *= -1;
        if(S[2] < 0) S[2] *= -1;
    }
    
    if(U.determinant() <= 0)
    {
        U(0, 2) *= -1;
        U(1, 2) *= -1;
        U(2, 2) *= -1;
    }
    
    if(W.determinant() <= 0)
    {
        W(0, 2) *= -1;
        W(1, 2) *= -1;
        W(2, 2) *= -1;
    }
    
    //TODO: energy model gradient 


    // see readMe < The Gradient of Principal Stretch Models >
    // lecture 6 slide 63
    Eigen::Vector3d dpsiDs;
    dpsiDs.setZero();
    dpsiDs[0] = 2.0 * mu * (S[0] - 1.0) + lambda * (S[0] + S[1] + S[2] - 3.0);
    dpsiDs[1] = 2.0 * mu * (S[1] - 1.0) + lambda * (S[0] + S[1] + S[2] - 3.0);
    dpsiDs[2] = 2.0 * mu * (S[2] - 1.0) + lambda * (S[0] + S[1] + S[2] - 3.0);
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> dS = U * dpsiDs.asDiagonal() * W.transpose();
    Eigen::Vector9d flatten_dS = Eigen::Map<const Eigen::Vector9d>(dS.data(), dS.size());


    //constructure B and N
    // lecture 6 slide 83
    Eigen::Matrix99d B_j;
    B_j << dX(0, 0), 0, 0, dX(1, 0), 0, 0, dX(2, 0), 0, 0,
           dX(0, 1), 0, 0, dX(1, 1), 0, 0, dX(2, 1), 0, 0,
           dX(0, 2), 0, 0, dX(1, 2), 0, 0, dX(2, 2), 0, 0,
           0, dX(0, 0), 0, 0, dX(1, 0), 0, 0, dX(2, 0), 0,
           0, dX(0, 1), 0, 0, dX(1, 1), 0, 0, dX(2, 1), 0,
           0, dX(0, 2), 0, 0, dX(1, 2), 0, 0, dX(2, 2), 0,
           0, 0, dX(0, 0), 0, 0, dX(1, 0), 0, 0, dX(2, 0),
           0, 0, dX(0, 1), 0, 0, dX(1, 1), 0, 0, dX(2, 1),
           0, 0, dX(0, 2), 0, 0, dX(1, 2), 0, 0, dX(2, 2);

    Eigen::Matrix93d N_j;
    N_j.setZero();
    N_j.block(0, 0, 3, 1) = N;
    N_j.block(3, 1, 3, 1) = N;
    N_j.block(6, 2, 3, 1) = N;

    Eigen::Vector3d delta_x1 = x1 - x0;
    Eigen::Vector3d delta_x2 = x2 - x0;
    double n_til = delta_x1.cross(delta_x2).norm();

    // cross product
    // see slide 81 for a clear visual
    Eigen::Matrix3d change_delta_x1;
    change_delta_x1 <<       0,               -1.0 * delta_x1(2),       delta_x1[1],
                            delta_x1[2],              0,             -1.0 * delta_x1[0],
                       -1.0 * delta_x1[1],        delta_x1[0],               0;

    Eigen::Matrix3d change_delta_x2;
    change_delta_x2 <<          0,            -1.0 * delta_x2[2],         delta_x2[1],
                            delta_x2[2],             0,              -1.0 * delta_x2[0],
                       -1.0 * delta_x2[1],       delta_x2[0],                0;

    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d O = Eigen::Matrix3d::Zero();
    Eigen::Matrix39d tmp1, tmp2;
    tmp1 << -I, O, I;
    tmp2 << -I, I, O;

    Eigen::Matrix39d final_delta_x1, final_delta_x2;
    final_delta_x1 = change_delta_x1 * tmp1;
    final_delta_x2 = change_delta_x2 * tmp2;

    // formula in slide 82
    Eigen::Matrix39d dndq = 1.0 / n_til * (I - n * n.transpose()) * (final_delta_x1 - final_delta_x2);
    Eigen::MatrixXd dFdq = B_j + N_j * dndq;

    dV = 1.0 * area * dFdq.transpose() * flatten_dS;
}
