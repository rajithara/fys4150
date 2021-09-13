#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "utilities.hpp"

// Declerations
int u(int x);
std::vector<double> u_vectorize(std::vector<double> x, int N);
std::string vector_to_string(std::vector<double> vec, int N); // For testing purposes

// Main
int main() {

    // Produce results for problem 1
    int N = 100; // Number of x values (spacing)
    run_problem2(N);

    //std::vector<std::vector<double>> A = get_matrix_A(N);
    //std::vector<double> v = solve_general_tridiagonal_matrix_equation(A, N);

    return 0;
}