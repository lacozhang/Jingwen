#define BOOST_TEST_MODULE FBSDE
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <Eigen/Dense>
#include "../../Common/common.h"

BOOST_AUTO_TEST_CASE(GHRun)
{
    Eigen::VectorXd w, x;
    GaussHermite(5, w, x);
    std::cout << "w value " << w << std::endl;
    std::cout << "x value " << x << std::endl;
}