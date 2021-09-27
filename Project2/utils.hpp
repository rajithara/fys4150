#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <assert.h> 
#include <vector>
#include <armadillo>

// Returns the tridiagonal matrix A for a given integer size N. The values of the diagonal
// elements are assigned in accordance with the problem text.
arma::mat get_tridiagonal_matrix_A(int N) {

    double n = N - 1;
    double h = 1/n;
    double a = -1/(h*h);
    double d = 2/(h*h);

    arma::mat A = arma::mat(N, N).fill(0.);
    int j = 0;
    for (int i = 0; i < N; i++){
        // First row of A
        if (i == 0){
            A(i, 0) = d;
            A(i, 1) = a;
        }
        // Last row of A
        else if(i == N-1){
            A(i, i-1) = a;
            A(i, i) = d;
        }
        // All other rows
        else {
            A(i, j) = a;
            A(i, j+1) = d;
            A(i, j+2) = a;
            j++;
        }
    }
    return A;
}

// Produce the analytical solution for the eigenvalue problem for the index pair (lambda_i, v_i). 
// The results will be stored in arguments $lambda and $v passed by reference. 
void solve_analytical(int i, int N, double &lambda, arma::vec &v){
    assert(v.size() == N);
    assert((i >= 1) && (i <= N));

    double n = N - 1;
    double h = 1/n;
    double a = -1/(h*h);
    double d = 2/(h*h);
    const double pi = 2*acos(0.0);

    lambda = d + 2*a*cos((i*pi)/(N+1));
    for (int j = 1; j <= N; j++) {
        v(j-1) = sin((i*j*pi)/(N+1));
    }
}

// Sort eigenvalues and eigenvector pairs in ascending order according to the eigenvalues
// Adapted from https://stackoverflow.com/questions/46276477/index-a-matrix-using-a-vector-of-indices-in-armadillo-library
void sort_solution_pairs(arma::mat &eigenvectors, arma::vec &eigenvalues) {
    arma::uvec order = arma::sort_index(eigenvalues);
    eigenvectors = eigenvectors.cols(order(arma::span(0, eigenvectors.n_cols-1)));
    eigenvalues = arma::sort(eigenvalues);
}

// Compare a numerical solution against the analytical solution
void compare_numerical_and_analytical(const arma::mat eigenvectors, const arma::vec eigenvalues, int N) {
    double tol = 10e-4;
    bool comparison_failed = 0;
    for (int i = 0; i < N; i++) {
        double lambda = 0;
        arma::vec v(N);
        solve_analytical(i+1, N, lambda, v);
        v = arma::normalise(v);

        // The comparison of the vectors is based on the average deviation of 
        // each element; might be a suboptimal solution. The comparison also allows for oppsite signed vectors.
        // 
        if ((abs(lambda - eigenvalues(i)) > tol) || (arma::accu(arma::abs(eigenvectors.col(i) - v)) > tol*N)) {
            if ((abs(lambda - eigenvalues(i)) > tol) || (arma::accu(arma::abs(-eigenvectors.col(i) - v)) > tol*N)) {
                std::cout << "\nDiscrepancy for i = " << i+1 << ":\n" << "Numerical: " << eigenvalues(i) << " | ";
                std::cout << eigenvectors.col(i).as_row() << "Analytical: " << lambda << " | " << v.as_row() << std::endl;
                comparison_failed = 1;
            }
        }
        //break;
    }
    if (comparison_failed == 0){
        std::cout << "The numerical solution is in accordance with the analytical solution.\n";
    }
    else {
        std::cout << "Comparison failed!\n";
    }
}

// Returns the largest off-diagonal element and also saves its indices by reference to 
// the integer variables $r and $c. 
// Assuming symmetric matrix input A; only checking elements below the diagonal
double find_largest_off_diagonal_element(const arma::mat &A, int &k, int &l) {
    assert((A.is_square()) && (A.n_rows >= 1));

    // Initializing
    double max_elem = A(1,0);
    k = 1;
    l = 0;

    // Iterating through the lower half of the off-diagonal elements
    int col_counter = 0;
    for (int i = 1; i < A.n_rows; i++) {
        for (int j = 0; j <= col_counter; j++) {
            if (abs(A(i,j)) > max_elem) {
                max_elem = abs(A(i,j));
                k = i;
                l = j;
            }
        }
        col_counter++;
    }
    return max_elem;
}

// Unit testing
void test_find_largest_off_diagonal_element() {

    arma::mat B = arma::mat(4,4).eye();
    B(0,3) = 0.5;
    B(1,2) = -0.7;
    B(2,1) = -0.7;
    B(3,0) = 0.5;

    int r;
    int c;
    double max_off_diag = find_largest_off_diagonal_element(B, r, c);
    
    assert ((r == 2) && (c == 1) && (max_off_diag == 0.7));
    return;
}

// Perform a single Jacobi rotation
// (Based on the code from page 220 of the compendium)
void jacobi_rotate(arma::mat &A, arma::mat &R, int k, int l) {
    // Determine sin and cos
    double s, c;
    if (A(k,l) != 0.0) {
        double t, tau;
        tau = (A(l,l) - A(k,k))/(2*A(k,l));
        if (tau > 0) {
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else {
            t = -1.0/(-tau + sqrt(1.0 + tau*tau));
        }
        c = 1/sqrt(1 + t*t);
        s = c*t;
    }
    else {
        c = 1.0;
        s = 0.0;
    }

    // Performing rotation/ updating the elements
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0.0;
    A(l,k) = 0.0;

    for ( int i = 0; i < A.n_cols; i++ ) {
        // Updating A
        if ( i != k && i != l ) {
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);
        }
        // Update the rotation matrix
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }
    return;
}

// const A or not???
void jacobi_eigensolver(arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged) {
    
    // Initialize the rotation matrix
    arma::mat R = arma::mat(A.n_rows, A.n_cols).eye();

    // Perform rotations iteratively
    int k,l;
    double max_offdiag = find_largest_off_diagonal_element(A, k, l);
    while((max_offdiag > eps) && (iterations < maxiter)) {
        max_offdiag = find_largest_off_diagonal_element(A, k, l);
        jacobi_rotate(A, R, k, l);
        iterations++;
    }
    converged = (max_offdiag < eps) ? 1 : 0;

    // Write eigenvalues and eigenvectors
    eigenvalues = A.diag();
    eigenvectors = R;

    return;
}