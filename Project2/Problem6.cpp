#include "utils.hpp"

// Computing the number of iterations before convergence, I, as a function of N, the size of the matrix. 
int main() {
    // Initialize
    int N_min = 6;
    int N_max = 100;
    arma::mat I_vec = arma::mat(N_max-N_min + 1, 1).fill(0); // storing as matrix as pyarma don't have vec (why tough?)

    int I_idx = 0;
    double eps = 10e-8; // (Keeping this constant???)
    for (int N = 6; N <= N_max; N++){
        // Set up the tridiagonal matrix
        arma::mat A = get_tridiagonal_matrix_A(N);
        
        // Solve the eigenvalue problem by the Jacobi rotation algorithm
        int maxiter = N*N*N;
        int iterations = 0;
        bool converged = 0;
        arma::vec eigval;
        arma::mat eigvec;
        while (converged == 0) {
            jacobi_eigensolver(A, eps, eigval, eigvec, maxiter, iterations, converged);
            maxiter *= N;
        }
        // Store the number of iterations needed before convergence
        I_vec(I_idx, 0) = iterations;
        I_idx++;
    }

    // Save the values of the iterration function I(N) as a binary (Armadillo) file
    // as well information about N_min and N_max
    I_vec.save("iterations_vector.bin");

    arma::mat N_info(2, 1);
    N_info(0,0) = N_min;
    N_info(1,0) = N_max;
    N_info.save("N_information.bin");

    return 0;
}