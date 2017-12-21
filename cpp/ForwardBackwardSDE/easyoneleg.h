#pragma once

#ifndef __EASY_ONE_LEG_H__
#include <iostream>
#include <functional>
#include <cmath>
#include <Eigen/Dense>

void solve(int timestep, double range, double& yerror, double& zerror, double theta,
    std::function<void(double, double, double&)> func,
    std::function<void(double, double, double&)> yfunc,
    std::function<void(double, double, double&)> zfunc,
    const double ytrue, const double ztrue,
    bool y_approx = true, bool z_approx = true, bool global=false);

#endif // !__EASY_ONE_LEG_H__
