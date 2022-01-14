#include <velocity_filter_cloth_sphere.h>

//Input:
//  qdot - the 3nx1 generalized velocities of the cloth mesh
//  index - a list of collision vertex indices from the collision detector
//  normals - a list of collision normals from the collision detector
//Output:
//  qdot- the filtered 3nx1 generalized velocities

void velocity_filter_cloth_sphere(Eigen::VectorXd &qdot, const std::vector<unsigned int> &indices, 
                                  const std::vector<Eigen::Vector3d> &normals) {

    // lecture 6 slide 95 for reference
    for (int i = 0; i < indices.size(); i++) {
        int index = indices[i];
        Eigen::Vector3d normal = normals[i];
        Eigen::Vector3d v = qdot.segment(3 * index, 3);
        double a = 0.0;
        if (normal.transpose() * v < 0) {
            a = -1.0 * normal.transpose() * v;
        }
        qdot.segment(3 * index, 3) = v + a * normal;
    }

}