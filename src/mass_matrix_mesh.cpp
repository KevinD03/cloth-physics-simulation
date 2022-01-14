#include <mass_matrix_mesh.h>

//Input: 
//  q - generalized coordinates for the FEM system
//  V - the nx3 matrix of undeformed vertex positions
//  F - the mx3 matrix of triangle-vertex indices
//  density - the density of the cloth material
//  areas - the mx1 vector of undeformed triangle areas
//Output:
//  M - sparse mass matrix for the entire mesh

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, 
                         Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F,
                         double density, Eigen::Ref<const Eigen::VectorXd> areas) {
    // Assume h=1

    M.resize(q.rows(), q.rows());
    M.setZero();

    typedef Eigen::Triplet<double> Triplet;
    std::vector<Triplet> tripletList;


    // this is before the integration
    // After integration ahould be 9*9. see slide 83 Bj for visual reference
    Eigen::Matrix3d m_i;
    m_i << 1.0 / 6.0, 1.0 / 12.0, 1.0 / 12.0,
        1.0 / 12.0, 1.0 / 6.0, 1.0 / 12.0,
        1.0 / 12.0, 1.0 / 12.0, 1.0 / 6.0;

    // mgh but h is 1 here, no need to mulitple
    m_i = m_i * density;

    // then we need to get the M_0 
    // see slide 32

    for (int i = 0; i < F.rows(); i++) {
        //Eigen::MatrixXd M_for_each_tri;
        //M_for_each_tri.resize(9, 9);

        // now doing it for each M_i
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                tripletList.push_back(Triplet(3 * F(i, j), 3 * F(i, k), m_i(j, k) * areas(i)));
                tripletList.push_back(Triplet(3 * F(i, j), 3 * F(i, k), m_i(j, k) * areas(i)));
                tripletList.push_back(Triplet(3 * F(i, j), 3 * F(i, k), m_i(j, k) * areas(i)));

                tripletList.push_back(Triplet(3 * F(i, j) + 1, 3 * F(i, k) + 1, m_i(j, k) * areas(i)));
                tripletList.push_back(Triplet(3 * F(i, j) + 1, 3 * F(i, k) + 1, m_i(j, k) * areas(i)));
                tripletList.push_back(Triplet(3 * F(i, j) + 1, 3 * F(i, k) + 1, m_i(j, k) * areas(i)));

                tripletList.push_back(Triplet(3 * F(i, j) + 2, 3 * F(i, k) + 2, m_i(j, k) * areas(i)));
                tripletList.push_back(Triplet(3 * F(i, j) + 2, 3 * F(i, k) + 2, m_i(j, k) * areas(i)));
                tripletList.push_back(Triplet(3 * F(i, j) + 2, 3 * F(i, k) + 2, m_i(j, k) * areas(i)));

                /*tripletList.push_back(Triplet(3 * F(i, j), 3 * F(i, k), m_i(j, k) * areas(i)));
                tripletList.push_back(Triplet(3 * F(i, j), 3 * F(i, k) + 1, m_i(j, k) * areas(i)));
                tripletList.push_back(Triplet(3 * F(i, j), 3 * F(i, k) + 2, m_i(j, k) * areas(i)));

                tripletList.push_back(Triplet(3 * F(i, j) + 1, 3 * F(i, k), m_i(j, k) * areas(i)));
                tripletList.push_back(Triplet(3 * F(i, j) + 1, 3 * F(i, k) + 1, m_i(j, k) * areas(i)));
                tripletList.push_back(Triplet(3 * F(i, j) + 1, 3 * F(i, k) + 2, m_i(j, k) * areas(i)));

                tripletList.push_back(Triplet(3 * F(i, j) + 2, 3 * F(i, k), m_i(j, k) * areas(i)));
                tripletList.push_back(Triplet(3 * F(i, j) + 2, 3 * F(i, k) + 1, m_i(j, k) * areas(i)));
                tripletList.push_back(Triplet(3 * F(i, j) + 2, 3 * F(i, k) + 2, m_i(j, k) * areas(i)));*/
            }
        }
    }
    
    M.setFromTriplets(tripletList.begin(), tripletList.end());
   
}
 