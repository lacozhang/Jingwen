#include "easyoneleg.h"
#include "../Common/common.h"

/// func: function34
/// yfunc: finalvaluey
/// zfunc: valuez
void solve(int timestep, double xrange, double & yerror, double & zerror, double theta,
    std::function<void(double, double&)> func,
    std::function<void(double, double, double&)> yfunc,
    std::function<void(double, double, double&)> zfunc) {
    Eigen::MatrixXd y, z;
    Eigen::VectorXd t, x;
    Eigen::VectorXd w, a, yb_cache, zb_cache;
    Eigen::Vector4d yinterpl, xinterpl, zinterpl;
    double timesteplength, xsteplength, yb, zb;
    double f, expecty, expectf, expectyw, expectfw, expectz, x_point, y_i;
    int tindex = 0;
    const double yture = 0.5e0, zture = 0.25e0, x0 = 0.0e0;

    const double pi = 4.0e0 * std::atan(1.0e0);
    std::cout << "#steps    : " << timestep << std::endl;
    timesteplength = 1.0e0 / (double)timestep;
    std::cout << "step size : " << timesteplength << std::endl;

    const int k = 10;
    std::cout << "gauss hermite : " << k << std::endl;

    GaussHermite(k, w, a);
    std::cout << "w : " << std::endl << w << std::endl;
    std::cout << "a : " << std::endl << a << std::endl;
    yb_cache.resize(k);
    zb_cache.resize(k);
    yb_cache.setZero();
    zb_cache.setZero();

    if (std::fabs(theta - 0.5e0) < 1e-15) {
        xsteplength = std::pow(timesteplength, 0.75e0);
    }
    else {
        xsteplength = std::pow(timesteplength, 0.5e0);
    }

    // RangeOfX(x0, k, timestep, xrange);
    int xstep = std::ceil(xrange / xsteplength);
    const int IndexOffset = xstep;
    std::cout << "actual x range : " << xrange << std::endl;
    t.resize(timestep + 1);
    x.resize(2 * xstep + 1);
    x.setZero();
    t.setZero();

    t[0] = 0.0e0;
    for (int i = 1; i <= timestep; ++i) {
        t[i] = t[i - 1] + timesteplength;
    }

    if (t.hasNaN()) {
        std::cout << "time array has NaN" << std::endl;
        std::abort();
    }

    x[0 + IndexOffset] = 0.0e0;
    for (int xidx = -xstep; xidx <= xstep; ++xidx) {
        x[xidx + IndexOffset] = xidx;
    }
    x *= xsteplength;
    if (x.hasNaN()) {
        std::cout << "x array has NaN" << std::endl;
        std::abort();
    }

    y.resize(timestep + 1, 2 * xstep + 1);
    z.resize(timestep + 1, 2 * xstep + 1);
    y.setZero();
    z.setZero();

    for (int xidx = -xstep; xidx <= xstep; ++xidx) {
        yfunc(t[timestep], x[xidx + IndexOffset], y.coeffRef(timestep, xidx + IndexOffset));
    }

    if (y.hasNaN()) {
        std::cout << "y has NaN" << std::endl;
        std::abort();
    }

    for (int xidx = -xstep; xidx <= xstep; ++xidx) {
        zfunc(t[timestep], x[xidx + IndexOffset], z.coeffRef(timestep, xidx + IndexOffset));
    }

    if (z.hasNaN()) {
        std::cout << "z has NaN" << std::endl;
        std::abort();
    }

    const double rtpi = std::sqrt(pi);
    const double rt2det = std::sqrt(2.0e0*timesteplength);

    for (int timeindex = timestep - 1; timeindex >= 0; --timeindex) {
//        std::cout << "time step : " << timeindex << std::endl;
// #pragma omp parallel for
        for (int xindex = -xstep; xindex <= xstep; ++xindex) {
            expecty = 0.0e0;
            expecty = 0.0e0;
            expectyw = 0.0e0;
            expectfw = 0.0e0;
            expectz = 0.0e0;
            yb_cache.setZero();
            zb_cache.setZero();
#ifdef _DEBUG
            std::cout << "yb cache : " << yb_cache << std::endl;
            std::cout << "zb cache : " << zb_cache << std::endl;
#endif // _DEBUG
            tindex = timeindex + 1;
            for (int i = 0; i < k; ++i) {
                x_point = x[xindex + IndexOffset] + rt2det * a[i];
#ifdef _DEBUG
                std::cout << "x_point : " << x_point << std::endl;
#endif // _DEBUG

                InterpolationIndex(x, y, x_point, xsteplength, IndexOffset, timestep, xstep, tindex, xinterpl, yinterpl);
#ifdef _DEBUG
                std::cout << "x int : " << std::endl << xinterpl << std::endl;
                std::cout << "y int : " << std::endl << yinterpl << std::endl;
#endif // _DEBUG

                NewtonInterpolation(xinterpl, yinterpl, x_point, yb);
#ifdef _DEBUG
                std::cout << "y int val : " << yb << std::endl;
#endif // _DEBUG

                InterpolationIndex(x, z, x_point, xsteplength, IndexOffset, timestep, xstep, tindex, xinterpl, zinterpl);
#ifdef _DEBUG
                std::cout << "x int : " << std::endl << xinterpl << std::endl;
                std::cout << "z int : " << std::endl << zinterpl << std::endl;
#endif // _DEBUG
                NewtonInterpolation(xinterpl, zinterpl, x_point, zb);
#ifdef _DEBUG
                std::cout << "z int val : " << zb << std::endl;
#endif // _DEBUG
                yb_cache[i] = yb;
                zb_cache[i] = zb;
            }

#ifdef _DEBUG
            std::cout << "yb cache : " << yb_cache << std::endl;
            std::cout << "zb cache : " << zb_cache << std::endl;
#endif // _DEBUG


            for (int i = 0; i < k; ++i) {
                func(yb_cache[i], f);
                expecty += w[i] * yb_cache[i];
                expectyw += w[i] * yb_cache[i] * rt2det * a[i];
                expectfw += w[i] * f * rt2det * a[i];
                expectz += w[i] * zb_cache[i];
            }

            expecty /= std::sqrt(pi);
            expectyw /= std::sqrt(pi);
            expectfw /= std::sqrt(pi);
            expectz /= std::sqrt(pi);

            z(timeindex, xindex + IndexOffset) = (expectyw + (1.0e0 - theta)*timesteplength*expectfw - (1.0e0 - theta)*timesteplength*expectz) / (theta*timesteplength);

            y_i = y(tindex, xindex + IndexOffset);
#ifdef _DEBUG
            std::cout << "y_i : " << y_i << std::endl;
#endif // _DEBUG

            int cycidx = 0;
            do {
                expectf = 0.0e0;
                cycidx++;

                for (int i = 0; i < k; ++i) {
                    func((1.0e0 - theta)*yb_cache(i) + theta*y_i, f);
                    expectf += w[i] * f;
                }
                expectf /= rtpi;
                y(timeindex, xindex + IndexOffset) = expecty + timesteplength*expectf;

#ifdef _DEBUG
                std::cout << "new y_i " << y(timeindex, xindex + IndexOffset) << std::endl;
#endif // _DEBUG

                if (std::fabs(y_i - y(timeindex, xindex + IndexOffset)) <= 1e-14) {
                    y(timeindex, xindex + IndexOffset) = y_i;
                    break;
                }
                else if (cycidx > 1000000) {
                    std::cout << "time " << timeindex << std::endl;
                    std::cout << "xindex " << xindex << std::endl;
                    std::cout << "abs : " << std::fabs(y_i - y(timeindex, xindex + IndexOffset)) << std::endl;
                    break;
                }
                y_i = y(timeindex, xindex + IndexOffset);
            } while (true);

            if (timeindex > 0) {
                yfunc(t[timeindex], x[xindex + IndexOffset], y(timeindex, xindex + IndexOffset));
            }

        }
    }

    zerror = std::fabs(zture - z(0, IndexOffset));
    yerror = std::fabs(yture - y(0, IndexOffset));
}
