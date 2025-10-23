#include "test_time_reversibility.hpp"

#include <algorithm>

#include "core/math_utils.hpp"
#include "core/workspace.hpp"
#include "evolve/evolve_taylor.hpp"

RevTestResult test_time_reversibility(const Tridiag& T,
                                      const Eigen::VectorXcd& psi_init,
                                      double dx,
                                      double dt,
                                      int nsteps,
                                      int K) {
    RevTestResult r{};
    r.norm0 = l2_norm(psi_init, dx);

    Eigen::VectorXcd psi = psi_init;
    TaylorWorkspace workspace;
    workspace.resize(psi.size());

    auto update_drift = [&](double& cur_max, const Eigen::VectorXcd& v) {
        double n = l2_norm(v, dx);
        cur_max = std::max(cur_max, std::abs(n - 1.0));
    };

    double max_drift = 0.0;
    for (int i = 0; i < nsteps; ++i) {
        taylor_step_tridiag(T, psi, dt, K, workspace);
        if ((i & 0x3FF) == 0) {
            update_drift(max_drift, psi);
        }
    }
    r.norm_fwd = l2_norm(psi, dx);
    update_drift(max_drift, psi);

    for (int i = 0; i < nsteps; ++i) {
        taylor_step_tridiag(T, psi, -dt, K, workspace);
        if ((i & 0x3FF) == 0) {
            update_drift(max_drift, psi);
        }
    }
    r.norm_back = l2_norm(psi, dx);
    update_drift(max_drift, psi);

    r.err_back = (psi - psi_init).norm();
    r.max_drift = max_drift;
    return r;
}

