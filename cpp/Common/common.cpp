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
        std::exit(1);
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

void RangeOfX(double x0, int k, int timesteps, double & xrange) {
    Eigen::VectorXd w, a;
    GaussHermite(k, w, a);

    double timesteplength = 1.0 / (timesteps*1.0);
    double xsteplength = std::pow(timesteplength, 0.75);
    int xsteps = std::floor (xrange / xsteplength);

    Eigen::VectorXd time;
    time.resize(timesteps + 1);
    for (int i = 0; i <= timesteps; ++i) {
        time[i] = i;
    }
    time *= timesteplength;

    int INDEX_OFFSET = xsteps;
    Eigen::VectorXd x;
    x.resize(xsteps * 2 + 1);
    x.setZero();
    for (int i = -xsteps; i <= xsteps; ++i) {
        x[i + INDEX_OFFSET] = x0 + i*xsteplength;
    }

    Eigen::VectorXd x_need_base = std::sqrt(2.0*timesteplength) * a;
    int xneedindex = 0;
    for (int tidx = 0; tidx <= timesteps; ++tidx) {
        double x_max = (x_need_base.array() + x[xneedindex + INDEX_OFFSET]).array().abs().maxCoeff();

        xneedindex = (int)std::floor((x_max - x0) / xsteplength);
        if (xneedindex > xsteps) {
            std::cerr << "current range too small, please extend" << std::endl;
            return;
        }
    }

    xrange = x0 + xneedindex + xsteplength;
}

void NewtonInterpolation(const Eigen::VectorXd & x, const Eigen::VectorXd & y, const double & x_aim, double & y_aim) {
    Eigen::MatrixXd z;
    Eigen::VectorXd x_;
    size_t sizex = x.size();
    z.resize(sizex, sizex);
    x_.resize(sizex);
    x_ = x;
    z.setZero();
    double n = 0.0;
    z.col(0) = y;

#ifdef _DEBUG
    std::cout << "x_ : " << x_ << std::endl;
    std::cout << "z  : " << z << std::endl;
#endif // _DEBUG

    for (int temp2 = 1; temp2 < sizex; ++temp2) {
        for (int temp1 = temp2; temp1 < sizex; ++temp1) {
            z(temp1, temp2) = (z(temp1 - 1, temp2 - 1) - z(temp1, temp2 - 1)) / (x_(temp1 - temp2) - x_(temp1));
#ifdef _DEBUG
            std::cout << "z  : " << std::endl << z << std::endl;
            assert(z.allFinite());
#endif // _DEBUG
        }
    }

#ifdef _DEBUG
    std::cout << "x_ : " << x_ << std::endl;
    std::cout << "z  : " << std::endl << z << std::endl;
    assert(z.allFinite());
#endif // _DEBUG

    for (int i = sizex - 1; i > 0; --i)
        x_[i] = x_[i - 1];

#ifdef _DEBUG
    std::cout << "after move : x_ " << x_ << std::endl;
#endif // _DEBUG


#ifdef _DEBUG
    std::cout << "after assign x_ " << std::endl;
    std::cout << "x_ : " << x_ << std::endl;
    std::cout << "z  : " << std::endl << z << std::endl;
#endif // _DEBUG
    x_[0] = x_aim - 1.0;
    x_ = x_aim - x_.array();

#ifdef _DEBUG
    std::cout << "after transform x_ : " << x_ << std::endl;
#endif // _DEBUG


    for (int temp1 = 1; temp1 < sizex; ++temp1) {
        x_[temp1] = x_[temp1 - 1] * x_[temp1];
    }
#ifdef _DEBUG
    assert(!x_.hasNaN());
#endif // _DEBUG


    for (int temp1 = 0; temp1 < sizex; ++temp1) {
        n += z(temp1, temp1) * x_[temp1];
    }

#ifdef _DEBUG
    assert(!std::isnan(n));
#endif // _DEBUG


    y_aim = n;
}

void LeastSquare(Eigen::VectorXi & timestep, Eigen::VectorXd & error, double & order) {
    int sizeerror = error.size();
    Eigen::VectorXd timesteplength;
    timesteplength.resize(sizeerror);

    timesteplength.setZero();
    for (int i = 0; i < sizeerror; ++i) {
        timesteplength[i] = std::log(1.0 / (timestep[i] * 1.0));
    }

    for (int i = 0; i < sizeerror; ++i) {
        error[i] = std::log(error[i]);
    }

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> a;
    a.resize(sizeerror, 2);
    a.setOnes();
    a.col(1) = timesteplength;
#ifdef _DEBUG
    std::cout << "matrix " << std::endl << a << std::endl;
    std::cout << "error " << std::endl << error << std::endl;
#endif // _DEBUG

    LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', sizeerror, 2, 1, a.data(), sizeerror, error.data(), sizeerror);
#ifdef _DEBUG
    std::cout << "matrix " << a << std::endl;
    std::cout << "result " << std::endl << error << std::endl;
#endif // _DEBUG

    // LAPACKE_dgelss(LAPACK_ROW_MAJOR, sizeerror, 2, 1, a.data(), sizeerror, error.data(), 1, s.data(), 0.01, &rank);
    order = error[1];
}

void InterpolationIndex(Eigen::VectorXd & xvalues, Eigen::MatrixXd & yvalues,
    double xpoint, double xsteplength, const int INDEX_OFFSET,
    int timesteps, int xsteps, int timeindex,
    Eigen::Vector4d & xintvalues, Eigen::Vector4d & yintvalues) {
    bool overflow = false;

    Eigen::Vector4i xindex;
    xindex[2] = (int)std::floor(xpoint / xsteplength);
    xindex[3] = xindex[2] + 1;
    xindex[1] = xindex[2] - 1;
    xindex[0] = xindex[2] - 2;

    if (xindex[0] < -xsteps) {
        xindex[0] = -xsteps;
        xindex[1] = xindex[0] + 1;
        xindex[2] = xindex[0] + 2;
        xindex[3] = xindex[0] + 3;
        overflow = true;
    }
    else if (xindex[3] > xsteps) {
        xindex[3] = xsteps;
        xindex[2] = xindex[3] - 1;
        xindex[1] = xindex[3] - 2;
        xindex[0] = xindex[3] - 3;
        overflow = true;
    }

    xintvalues[0] = xvalues[xindex[0] + INDEX_OFFSET];
    xintvalues[1] = xvalues[xindex[1] + INDEX_OFFSET];
    xintvalues[2] = xvalues[xindex[2] + INDEX_OFFSET];
    xintvalues[3] = xvalues[xindex[3] + INDEX_OFFSET];

    yintvalues[0] = yvalues(timeindex, xindex[0] + INDEX_OFFSET);
    yintvalues[1] = yvalues(timeindex, xindex[1] + INDEX_OFFSET);
    yintvalues[2] = yvalues(timeindex, xindex[2] + INDEX_OFFSET);
    yintvalues[3] = yvalues(timeindex, xindex[3] + INDEX_OFFSET);

    if (xvalues[xindex[0] + INDEX_OFFSET] > xpoint && !overflow) {
        std::cout << "error when compute lower boundary" << std::endl;
    }

    if (xvalues[xindex[3] + INDEX_OFFSET] < xpoint && !overflow) {
        std::cout << "error when compute upper boundary" << std::endl;
    }
}