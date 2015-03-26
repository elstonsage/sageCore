#ifndef LINALG_H

#include "numerics/matrix.h"

namespace SAGE {

double cnorm(const Matrix2D<double>& m, size_t col, size_t row = 0);

void cholesky(const Matrix2D<double>& A, Matrix2D<double>& G);

Matrix2D<double>& msqrt(const Matrix2D<double>& XX, Matrix2D<double>& X);
Matrix2D<double>  msqrt(const Matrix2D<double>& XX);

void qr(const Matrix2D<double>& A, Matrix2D<double>& Q, Matrix2D<double>& R);

}

#endif
