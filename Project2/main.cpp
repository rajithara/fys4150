#include <iostream>
#include <armadillo>

#include "Matrix.hpp"
#include "utils.hpp"

int main() {

    int N = 6;
    arma::mat A = get_tridiagonal_matrix_A(N);
    A.print("A");

    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);

    eigval.print("\nEigen values:");
    eigvec.print("\nEigen vectors:");

    arma::mat B = arma::mat(4,4).eye();
    B(0,3) = 0.5;
    B(1,2) = -0.7;
    B(2,1) = -0.7;
    B(3,0) = 0.5;

    int r;
    int c;
    double max_off_diag = find_largest_off_diagonal_element(B, r, c);
    std::cout << "\nMax off diagonal element: " << max_off_diag << "\n(r,c) = " << r << "," << c << std::endl;;
    

    return 0;
}