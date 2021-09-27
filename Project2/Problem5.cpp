#include "utils.hpp"

// Compare the numerical Jaboci method with the analytical solution
int main(){

    // Set up the tridiagonal matrix
    int N = 6;
    arma::mat A = get_tridiagonal_matrix_A(N);

    // Solve the eigenvalue problem by the Jacobi rotation algorithm
    double eps = 10e-8;
    int maxiter = 1000;
    int iterations = 0;
    bool converged = 0;

    arma::vec eigval;
    arma::mat eigvec;
    jacobi_eigensolver(A, eps, eigval, eigvec, maxiter, iterations, converged);
    std::string message = (converged == 1) ? "Convergence was reached.\n" : "Convergence was not reached. \n";
    std::cout << message;

    // Sort the generated solution pairs
    sort_solution_pairs(eigvec, eigval);

    // Compare results agianst the analytical solution
    std::cout << "Comparing Jacobi's rotation algorithm to the analytical solution:\n";
    compare_numerical_and_analytical(arma::normalise(eigvec), eigval, N);

    return 0;
}