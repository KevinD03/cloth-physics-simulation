#include <V_spring_particle_particle.h>

void V_spring_particle_particle(double &V, Eigen ::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
    // code from A3
    
    // see the lagrangian formula in lecture 3 page 33
    // 
    // Set up the Matrix B formed by combine -I and I
    // No need to set up B since it's just q1-q0 3D distance

    /*Eigen::MatrixXd B(3, 6);
    B << -1., 0., 0., 1., 0., 0.,
        0., -1., 0., 0., 1., 0.,
        0., 0., -1., 0., 0., 1.;
    Eigen::VectorXd q(6, 1);
    q << q0, q1;
    V = 0.5 * stiffness * pow((sqrt(q.transpose() * B.transpose() * B * q)), 2);*/

    double deformed_length;
    deformed_length = (q1 - q0).norm();

    V = 0.5 * stiffness * std::pow(deformed_length - l0, 2);
}