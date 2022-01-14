#include <collision_detection_cloth_sphere.h>
#include <iostream>

//  q - generalized coordinates for the FEM system
//  center - the position of the sphere center in the world space
//  radius - the radius of the sphere in the world space 
//Output:
//  cloth_index - the indices of all vertices currently in contact with the sphere
//  normals - the outward facing contact normals for each contacting vertex. 

void collision_detection_cloth_sphere(std::vector<unsigned int> &cloth_index, std::vector<Eigen::Vector3d> &normals, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Vector3d> center, double radius) {

    cloth_index.clear();
    normals.clear();

    Eigen::Vector3d cur_q;

    // see lecture 6 slide 93 for formula
    for (int i = 0; i < q.size() / 3; i++) {

        cur_q[0] = q[3.0 * i];
        cur_q[1] = q[3.0 * i + 1.0];
        cur_q[2] = q[3.0 * i + 2.0];
        
        double x_center_dis;
        x_center_dis = std::pow((cur_q[0] - center[0]), 2.0) + std::pow((cur_q[1] - center[1]), 2.0) + std::pow((cur_q[2] - center[2]), 2.0);
        // collision happened
        if (x_center_dis <= radius*radius) {
            cloth_index.push_back(i);
            Eigen::Vector3d tmp_q;
            tmp_q = (cur_q - center).normalized();
            normals.push_back(tmp_q);
        }
    }
}