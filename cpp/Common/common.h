#pragma once
#ifndef __COMMON_H__
#include <iostream>
#include <functional>
#include <Eigen/Dense>

void GaussHermite(int n, Eigen::VectorXd&w, Eigen::VectorXd& x);
void RangeOfX(double x0, int k, int timesteps, double& xrange);
void NewtonInterpolation(const Eigen::VectorXd& x, const Eigen::VectorXd& y, const double& x_aim, double& y_aim);
void LeastSquare(Eigen::VectorXi& timestep, Eigen::VectorXd& error, double& order);
void InterpolationIndex(Eigen::VectorXd& xvalues, Eigen::MatrixXd& yvalues,
    double xpoint, double xsteplength, const int INDEX_OFFSET,
    int timesteps, int xsteps, int timeindex,
    Eigen::Vector4d& xintvalues, Eigen::Vector4d& yintvalues);


#endif // !__COMMON_H__
