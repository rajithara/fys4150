#include "utils.hpp"

// Produce solution for N=10 and N=100, and save the eigenvectors corrosp. to the 
// three lowest eigenvalues as binary files.
int main(){
    // Initiate
    arma::mat eigvec_9 = arma::mat(9, 3);
    arma::mat eigvec_99 = arma::mat(99, 3);

    for (int N = 9; N <= 99;) {
        // Set up the tridiagonal matrix
        arma::mat A = get_tridiagonal_matrix_A(N);

        // Solve the eigenvalue problem by the Jacobi rotation algorithm
        double eps = 10e-8;
        int maxiter = N*N*N;
        int iterations = 0;
        bool converged = 0;

        arma::vec eigval;
        arma::mat eigvec;
        jacobi_eigensolver(A, eps, eigval, eigvec, maxiter, iterations, converged);

        // Sort the generated solution pairs
        sort_solution_pairs(eigvec, eigval);

        // Extract the eigenvectors corrsp. to the three lowest eigenvalues
        if (N==9) {
            eigvec_9 = eigvec.submat(0, 0, 8, 2);
        }
        else {
            eigvec_99 = eigvec.submat(0, 0, 98, 2);
        }

        // Update the iterator
        N += 90;
    }

    // Save as bianry files
    eigvec_9.save("eigvec_9.bin");
    eigvec_99.save("eigvec_99.bin");

    return 0;
}