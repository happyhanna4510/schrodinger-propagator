#include "core/spectral.hpp"

#include <Eigen/Eigenvalues>

#include <complex>
#include <stdexcept>

#include "hamiltonian.hpp"

using cplx = std::complex<double>;

SpectralData make_spectral_data(const Eigen::MatrixXd& H,
                                const Eigen::VectorXcd& psi_init,
                                double dx) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
    if (es.info() != Eigen::Success) {
        throw std::runtime_error("eigensolve failed");
    }

    Eigen::MatrixXd Vreal = es.eigenvectors();
    renormalize_eigenvectors(Vreal, dx);

    SpectralData data;
    data.evals = es.eigenvalues();
    data.eigenvectors = Vreal.cast<cplx>();

    data.phi0 = data.eigenvectors.col(0);
    if (data.phi0(0).real() < 0.0) {
        data.phi0 = -data.phi0;
        data.eigenvectors.col(0) = -data.eigenvectors.col(0);
    }
    data.E0 = data.evals(0);

    data.coeffs0 = dx * data.eigenvectors.adjoint() * psi_init;

    return data;
}

