#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {
    // code from A3

    f(0) = stiffness * (q0(0) * 2.0 - q1(0) * 2.0) * (l0 - sqrt(q0(0) * (q0(0) - q1(0)) + q0(1) * (q0(1) - q1(1)) - q1(0) * (q0(0) - q1(0)) + q0(2) * (q0(2) - q1(2)) - q1(1) * (q0(1) - q1(1)) - q1(2) * (q0(2) - q1(2)))) * 1.0 / sqrt(q0(0) * (q0(0) - q1(0)) + q0(1) * (q0(1) - q1(1)) - q1(0) * (q0(0) - q1(0)) + q0(2) * (q0(2) - q1(2)) - q1(1) * (q0(1) - q1(1)) - q1(2) * (q0(2) - q1(2))) * (-1.0 / 2.0);
    f(1) = stiffness * (q0(1) * 2.0 - q1(1) * 2.0) * (l0 - sqrt(q0(0) * (q0(0) - q1(0)) + q0(1) * (q0(1) - q1(1)) - q1(0) * (q0(0) - q1(0)) + q0(2) * (q0(2) - q1(2)) - q1(1) * (q0(1) - q1(1)) - q1(2) * (q0(2) - q1(2)))) * 1.0 / sqrt(q0(0) * (q0(0) - q1(0)) + q0(1) * (q0(1) - q1(1)) - q1(0) * (q0(0) - q1(0)) + q0(2) * (q0(2) - q1(2)) - q1(1) * (q0(1) - q1(1)) - q1(2) * (q0(2) - q1(2))) * (-1.0 / 2.0);
    f(2) = stiffness * (q0(2) * 2.0 - q1(2) * 2.0) * (l0 - sqrt(q0(0) * (q0(0) - q1(0)) + q0(1) * (q0(1) - q1(1)) - q1(0) * (q0(0) - q1(0)) + q0(2) * (q0(2) - q1(2)) - q1(1) * (q0(1) - q1(1)) - q1(2) * (q0(2) - q1(2)))) * 1.0 / sqrt(q0(0) * (q0(0) - q1(0)) + q0(1) * (q0(1) - q1(1)) - q1(0) * (q0(0) - q1(0)) + q0(2) * (q0(2) - q1(2)) - q1(1) * (q0(1) - q1(1)) - q1(2) * (q0(2) - q1(2))) * (-1.0 / 2.0);
    f(3) = (stiffness * (q0(0) * 2.0 - q1(0) * 2.0) * (l0 - sqrt(q0(0) * (q0(0) - q1(0)) + q0(1) * (q0(1) - q1(1)) - q1(0) * (q0(0) - q1(0)) + q0(2) * (q0(2) - q1(2)) - q1(1) * (q0(1) - q1(1)) - q1(2) * (q0(2) - q1(2)))) * 1.0 / sqrt(q0(0) * (q0(0) - q1(0)) + q0(1) * (q0(1) - q1(1)) - q1(0) * (q0(0) - q1(0)) + q0(2) * (q0(2) - q1(2)) - q1(1) * (q0(1) - q1(1)) - q1(2) * (q0(2) - q1(2)))) / 2.0;
    f(4) = (stiffness * (q0(1) * 2.0 - q1(1) * 2.0) * (l0 - sqrt(q0(0) * (q0(0) - q1(0)) + q0(1) * (q0(1) - q1(1)) - q1(0) * (q0(0) - q1(0)) + q0(2) * (q0(2) - q1(2)) - q1(1) * (q0(1) - q1(1)) - q1(2) * (q0(2) - q1(2)))) * 1.0 / sqrt(q0(0) * (q0(0) - q1(0)) + q0(1) * (q0(1) - q1(1)) - q1(0) * (q0(0) - q1(0)) + q0(2) * (q0(2) - q1(2)) - q1(1) * (q0(1) - q1(1)) - q1(2) * (q0(2) - q1(2)))) / 2.0;
    f(5) = (stiffness * (q0(2) * 2.0 - q1(2) * 2.0) * (l0 - sqrt(q0(0) * (q0(0) - q1(0)) + q0(1) * (q0(1) - q1(1)) - q1(0) * (q0(0) - q1(0)) + q0(2) * (q0(2) - q1(2)) - q1(1) * (q0(1) - q1(1)) - q1(2) * (q0(2) - q1(2)))) * 1.0 / sqrt(q0(0) * (q0(0) - q1(0)) + q0(1) * (q0(1) - q1(1)) - q1(0) * (q0(0) - q1(0)) + q0(2) * (q0(2) - q1(2)) - q1(1) * (q0(1) - q1(1)) - q1(2) * (q0(2) - q1(2)))) / 2.0;
}