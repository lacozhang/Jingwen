#include <lapacke.h>
#include "common.h"

/*
void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
    double *w, double *work, int *lwork, int *info);
*/

void GaussHermite(int n, Eigen::VectorXd & w, Eigen::VectorXd & x) {
    Eigen::VectorXd a(n - 1);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> cm(n, n);
    x.resize(n);
    w.resize(n);

    x.setZero();
    w.setZero();

    for (int i = 0; i < n - 1; ++i) {
        a[i] = std::sqrt(static_cast<double>(i+1)*0.5);
    }

#ifdef _DEBUG
    std::cout << "a coeff " << std::endl;
    std::cout << a << std::endl;
#endif // _DEBUG

    cm.setZero();
    for (int i = 0; i < n - 1; ++i) {
        cm.coeffRef(i, i + 1) = a[i];
    }

#ifdef _DEBUG
    std::cout << "cm matrix " << std::endl;
    std::cout << cm << std::endl;
#endif // _DEBUG

    /*
    char jobz = 'V', uplo = 'U';
    int cols = n, lda = n, lwork = -1, info = 0;
    double *work = static_cast<double*>(std::malloc(lwork * sizeof(double)));
    dsyev_(&jobz, &uplo, &cols, cm.data(), &lda, x.data(), work, &lwork, &info);    
    std::cout << "info " << info << std::endl;
    std::cout << "lwork " << lwork << std::endl;
    std::free(work);
    */

    lapack_int retinfo=0, cols = n, lda = n;

#ifdef _DEBUG
    std::cout << "cols " << cols << std::endl;
    std::cout << "lda  " << lda << std::endl;
#endif // _DEBUG

    retinfo = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', cols, cm.data(), lda, x.data());
    if (retinfo != 0) {
        std::cerr << "Error when call dsyev" << std::endl;
        std::cerr << "ret info " << retinfo << std::endl;
    }

    w = cm.row(0).array().pow(2.0);
    w *= std::sqrt(4.0*std::atan(1.0));

#ifdef _DEBUG
    std::cout << "eigenvalues " << std::endl;
    std::cout << x << std::endl;
    std::cout << "res matrix " << std::endl;
    std::cout << cm << std::endl;
#endif // _DEBUG

}
