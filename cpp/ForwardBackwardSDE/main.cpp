#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "../Common/common.h"
#include "easyoneleg.h"

namespace Function2 {
    void function34(double y, double z, double& f) {
        f = (0.5e0*y - z) / (std::pow(y, 2.0) + std::pow(z, 2.0));
    }

    void finalvaluey(double t, double x, double& y) {
        y = std::sin(x + t);
    }

    void valuez(double t, double x, double& z) {
        z = std::cos(x + t);
    }

    const static double y_true = 0.0, z_true = 1.0;
}

namespace Function3 {
    void function34(double y, double z, double& f) {
        f = -1.0e0 * std::pow(y, 3) + 2.5e0 * std::pow(y, 2.0e0) - 1.5e0 * y;
    }

    void finalvaluey(double t, double x, double& y) {
        double expo = t + x;
        if (expo > 400) {
            y = 1.0;
        }
        else if (expo < -400) {
            y = 0.0;
        }
        else
            y = std::exp(t + x) / (std::exp(t + x) + 1.0e0);
    }

    void valuez(double t, double x, double& z) {
        z = std::exp(x + t) / std::pow(std::exp(x + t) + 1.0e0, 2.0);
    }

    const static double y_true = 0.5e0, z_true = 0.25e0;
}

namespace LinearFunction {

    void function34(double y, double z, double& f) {
        f = 0.5*y - z;
    }

    void finalvaluey(double t, double x, double& y) {
        y = std::sin(x + t);
    }

    void valuez(double t, double x, double& z) {
        z = std::cos(x + t);
    }

    const static double y_true = 0.0, z_true = 1.0;
}


int main(int argc, const char* argv[]) {
    Eigen::VectorXi timesteps;
    timesteps.resize(4);
    timesteps.setZero();
    timesteps[0] = 8;
    for (int i = 1; i < timesteps.size(); ++i) {
        timesteps[i] = timesteps[i - 1] * 2;
    }

    Eigen::VectorXd yerror, zerror;
    yerror.resize(timesteps.size());
    zerror.resize(timesteps.size());
    yerror.setZero();
    zerror.setZero();

    double range = 50.0, theta = 0.5;

    std::cout << "Please input theta : ";
    std::cin >> theta;
    std::cout << std::endl;

    std::cout << "Summary of input " << std::endl;
    std::cout << "range : " << range << " theta : " << theta << std::endl;

    namespace target = Function3;

    for (int i = 0; i < timesteps.size(); ++i) {
        solve(timesteps[i], range, yerror[i], zerror[i], theta,
            target::function34, target::finalvaluey, target::valuez,
            target::y_true, target::z_true,
            true, true, false);
    }

    std::cout << "y errors : " << std::endl << yerror << std::endl;
    std::cout << "z errors : " << std::endl << zerror << std::endl;
    double order = 0.0;
    LeastSquare(timesteps, yerror, order);
    std::cout << "order of y " << order << std::endl;
    order = 0.0;
    LeastSquare(timesteps, zerror, order);
    std::cout << "order of z " << order << std::endl;
}