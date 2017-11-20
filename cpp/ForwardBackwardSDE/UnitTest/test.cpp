#define BOOST_TEST_MODULE FBSDE
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <Eigen/Dense>
#include "../../Common/common.h"

BOOST_AUTO_TEST_CASE(GHRun)
{
    Eigen::VectorXd w, x;
    GaussHermite(10, w, x);
    std::cout << "w value " << w << std::endl;
    std::cout << "x value " << x << std::endl;
}

BOOST_AUTO_TEST_CASE(Newton) {
    Eigen::VectorXd x, y;
    x.resize(5);
    y.resize(5);
    x.setZero();
    y.setZero();
    x[0] = -1.5, y[0] = -14.1014;
    x[1] = -0.75, y[1] = -0.931596;
    x[2] = 0, y[2] = 0;
    x[3] = 0.75, y[3] = 0.931596;
    x[4] = 1.5, y[4] = 14.1014;
    double x_aim = 2.7, y_aim = 0;
    NewtonInterpolation(x, y, x_aim, y_aim);
    std::cout << std::fabs(y_aim - 91.17478282e0) << std::endl;
    BOOST_ASSERT(std::fabs(y_aim - 91.17478282e0) < 1);
}

BOOST_AUTO_TEST_CASE(LS) {
    Eigen::VectorXi timesteps;
    timesteps.resize(5);
    timesteps[0] = 4;
    for (int i = 1; i < timesteps.size(); ++i) {
        timesteps[i] = timesteps[i - 1] * 2;
    }

    double order = 0.0;
    Eigen::VectorXd error;
    error.resize(timesteps.size());
    error[0] = 1.998e-2;
    error[1] = 9.961e-3;
    error[2] = 4.979e-3;
    error[3] = 2.472e-3;
    error[4] = 1.247e-3;

    std::cout << "start run " << std::endl;
    LeastSquare(timesteps, error, order);
    std::cout << "order " << order << std::endl;
}