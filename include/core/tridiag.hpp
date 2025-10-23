#pragma once

#include <Eigen/Core>

struct Tridiag {
    Eigen::VectorXd a;  // main   (size M)
    Eigen::VectorXd b;  // upper  (size M-1)
    Eigen::VectorXd c;  // lower  (size M-1)
};

