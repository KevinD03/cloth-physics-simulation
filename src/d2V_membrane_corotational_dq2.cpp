#include <d2V_membrane_corotational_dq2.h>
#include <iostream>

//  q - generalized coordinates for the FEM system
//  dX - the 3x3 matrix containing dphi/dX
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the vertex indices of this triangle
//  area - the area of this triangle
//  mu,lambda - material parameters for the cloth material model
//Output:
//  H - the per-triangle Hessian of the potential energy (the linear model described in the README).

void d2V_membrane_corotational_dq2(Eigen::Matrix99d &H, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {
    

    //SVD = USW^T
    Eigen::Matrix3d U;
    Eigen::Vector3d S; 
    Eigen::Matrix3d W; 
    Eigen::Matrix3d F; //deformation gradient
    
    double tol = 1e-5;
    
    //Compute SVD of F here

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

    F = x * X;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = svd.matrixU();
    W = svd.matrixV();
    S = svd.singularValues();

    //deal with singularity in the svd gradient
    if(std::fabs(S[0] - S[1]) < tol || std::fabs(S[1] - S[2]) < tol || std::fabs(S[0] - S[2]) < tol) {
        F += Eigen::Matrix3d::Random()*tol;
        Eigen::JacobiSVD<Eigen::Matrix3d> svd2(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        U = svd2.matrixU();
        W = svd2.matrixV();
        S = svd2.singularValues();
    }
    
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

    //TODO: compute H, the hessian of the corotational energy

    // see readMe < The Hessian of Principal Stretch Models >
    Eigen::Tensor3333d dU;
    Eigen::Tensor3333d dW;
    Eigen::Tensor333d dS;
    dsvd(dU, dS, dW, F);

    // lecture 6 slide 78
    Eigen::Vector3d dpsiDs;
    dpsiDs.setZero();
    dpsiDs[0] = 2.0 * mu * (S[0] - 1.0) + lambda * (S[0] + S[1] + S[2] - 3.0);
    dpsiDs[1] = 2.0 * mu * (S[1] - 1.0) + lambda * (S[0] + S[1] + S[2] - 3.0);
    dpsiDs[2] = 2.0 * mu * (S[2] - 1.0) + lambda * (S[0] + S[1] + S[2] - 3.0);
  
    // see readMe < The Hessian of Principal Stretch Models >   ds_ij part
    Eigen::Matrix3d dpsiDs_M;
    dpsiDs_M.setZero();
    Eigen::Matrix3d d2psiDs;
    d2psiDs.setZero();

    // lecture 6 slide 75
    dpsiDs_M(0, 0) = dpsiDs[0];
    dpsiDs_M(1, 1) = dpsiDs[1];
    dpsiDs_M(2, 2) = dpsiDs[2];

    d2psiDs(0, 0) = lambda + 2 * mu;
    d2psiDs(0, 1) = lambda;
    d2psiDs(0, 2) = lambda;
    d2psiDs(1, 0) = lambda;
    d2psiDs(1, 1) = lambda + 2 * mu;
    d2psiDs(1, 2) = lambda;
    d2psiDs(2, 0) = lambda;
    d2psiDs(2, 1) = lambda;
    d2psiDs(2, 2) = lambda + 2 * mu;
 
    // lecture 6 slide 75
    Eigen::MatrixXd d2psi_dF;
    d2psi_dF.resize(9, 9);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Eigen::Vector3d S_ij = d2psiDs * dS[i][j];
            
            Eigen::Matrix3d d2psi_dF_ij;
            d2psi_dF_ij.setZero();
            // slide 75 for formula
            d2psi_dF_ij = dU[i][j] * dpsiDs_M * W.transpose() + U * S_ij.asDiagonal() * W.transpose() + U * dpsiDs_M * dW[i][j].transpose();
          
            for (int m = 0; m < 3; m++) {
                d2psi_dF.block(3 * i + j, m*3, 1, 3) = d2psi_dF_ij.row(m);
            }
        }
    }

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
    change_delta_x1 << 0, -1.0 * delta_x1(2), delta_x1[1],
        delta_x1[2], 0, -1.0 * delta_x1[0],
        -1.0 * delta_x1[1], delta_x1[0], 0;

    Eigen::Matrix3d change_delta_x2;
    change_delta_x2 << 0, -1.0 * delta_x2[2], delta_x2[1],
        delta_x2[2], 0, -1.0 * delta_x2[0],
        -1.0 * delta_x2[1], delta_x2[0], 0;

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

    H = 1.0 * area * dFdq.transpose() * d2psi_dF * dFdq;


    //fix errant eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix99d> es(H);
    
    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();
    
    for (int i = 0; i < 9; ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }
    
    H = Evec * DiagEval * Evec.transpose();
    
}
