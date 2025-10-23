#pragma once

#include <Eigen/Core>

struct TaylorWorkspace {
    void resize(Eigen::Index n) {
        if (sum_.size() != n) {
            sum_.resize(n);
            vk_.resize(n);
            tmp_.resize(n);
        }
    }

    Eigen::VectorXcd& sum() { return sum_; }
    Eigen::VectorXcd& vk() { return vk_; }
    Eigen::VectorXcd& tmp() { return tmp_; }

private:
    Eigen::VectorXcd sum_;
    Eigen::VectorXcd vk_;
    Eigen::VectorXcd tmp_;
};

