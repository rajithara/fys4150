#include "utils.hpp"

int main() {

    // Set up the tridiagonal matrix
    int N = 6;
    arma::mat A = get_tridiagonal_matrix_A(N);
    A.save("output_setting_up_A.txt", arma::raw_ascii);

    // Solve the eigenvalue problem using armadillo (and normalize)
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);
    //eigvec = arma::normalise(eigvec);

    // Check armadillo's solution against the analytical
    std::cout << "Comparing Armadillo's solution to the analytical solution:\n";
    compare_numerical_and_analytical(eigvec, eigval, N);
    return 0;
}