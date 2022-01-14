#include <fixed_point_constraints.h>
#include <algorithm>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
	// code from A3
	P.resize(q_size - 3 * indices.size(), q_size);
	P.setZero();

	int p_row = 0;

	for (int i = 0; i < q_size; i += 3) {
		if (std::count(indices.begin(), indices.end(), (i / 3)) == 0) {
			// update corresponding x,y,z to 1
			// 1 0 0 0 ......
			// 0 1 0 0 ......
			// 0 0 1 0 ......
			for (int j = 0; j < 3; j++) {
				P.insert(p_row, i + j) = 1;
				p_row++;
			}
		}
	}
}