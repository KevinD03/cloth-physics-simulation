#include <dV_cloth_gravity_dq.h>

//  M - sparse mass matrix for the entire mesh
//  g - the acceleration due to gravity
//Output:
//  fg - the gradient of the gravitational potential for the entire mesh

void dV_cloth_gravity_dq(Eigen::VectorXd &fg, Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::Vector3d> g) {
    
    fg.setZero();
    
    // apply g on all x,y,z of all vertices
    Eigen::VectorXd g_mesh;
    g_mesh.resize(M.rows());

    for (int i = 0; i < M.rows()/3 ; i++) {
        g_mesh.segment(i * 3, 3) = g;
    }
    
    fg = -M * g_mesh;
}
