#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>

// Declerations
double f(double x);
double u(double x);
std::vector<double> get_vector_x(int N);
std::vector<double> u_vectorize(std::vector<double> x, int N);
void run_problem2(int N);
std::string vector_to_string(std::vector<double> vec, int N); // For testing purposes

std::vector<std::vector<double>> get_matrix_A(int N);
std::vector<double> solve_general_tridiagonal_matrix_equation(std::vector<std::vector<double>> A, int N);

// Definitions
void run_problem2(int N) {
    // Compute values for u(x)
    std::vector<double> x = get_vector_x(N);
    std::vector<double> res = u_vectorize(x, N);

    // Store to a .txt file in the format "x u(x)"
    std::string filename = "problem2_output.txt";
    std::ofstream ofile;
    ofile.open(filename);
    for (int i = 0; i < N; i++) {
        ofile << std::setprecision(4) << std::scientific << x.at(i) << " " << res.at(i) << std::endl;
    }
    ofile.close();
}

double u(double x) {
    return 1. - (1. - exp(-10.))*x - exp(-10.*x);
}

double f(double x) {
    return 100*exp(-10*x);
}

std::vector<double> get_vector_x(int N) {

    // Set domain boundaries
    double x1 = 0.0;
    double xN = 1.0;

    // Define x values
    std::vector<double> x(N, 0.0);
    double h = (xN - xN)/(N-1);
    for (int i = 0; i < N; i++) {
        x.at(i) = h*i;
    };

    return x;
}

std::vector<double> u_vectorize(std::vector<double> x, int N) {
    std::vector<double> res(N, 0.0);
    for (int i = 1; i < N-1; i++){ // The boundaries are fixed to zero
        res.at(i) = u(x.at(i));
    }    
    return res;
}

std::string vector_to_string(std::vector<double> vec, int N) { // For testing purposes
    std::string res;
    for (int i = 0; i < N; i++) {
        res.append(std::to_string(vec.at(i)));
        res.append(" ");
    };
    return res;
}

std::vector<std::vector<double>> get_matrix_A(int N) {
    std::vector<std::vector<double>> A(N);
    int j = 0;
    for (int i = 0; i < N; i++){
        // First row of A
        if (i == 0){
            std::vector<double> tmp(N, 0.0);
            tmp.at(0) = 2;
            tmp.at(1) = -1;
            A.at(i) = tmp;
        }
        // Last row of A
        else if(i == N-1){
            std::vector<double> tmp(N, 0.0);
            tmp.at(N-2) = -1;
            tmp.at(N-1) = 2;
            A.at(i) = tmp;
        }
        // All other rows
        else {
            std::vector<double> tmp(N, 0.0);
            tmp.at(j) = -1;
            tmp.at(j + 1) = 2;
            tmp.at(j + 2) = -1;
            A.at(i) = tmp;
            j++;
        }
    }
    return A;
}

std::vector<double> solve_general_tridiagonal_matrix_equation(std::vector<std::vector<double>> A, int N) {

    // Extract diagonal matrices from A
    std::vector<double> a(N-1);
    std::vector<double> b(N);
    std::vector<double> c(N-1);

    int j = 0;
    for (int i = 0; i < N; i++){
        // First row of A
        if (i == 0){
            a.at(0) = A.at(i).at(1);
            b.at(0) = A.at(i).at(0);
        }
        // Last row of A
        else if(i == N-1){
            b.at(i) = A.at(i).at(i);
            c.at(i-1) = A.at(i).at(i-1);
        }
        // All other rows
        else {
            a.at(i) = A.at(i).at(j);
            b.at(i) = A.at(i).at(j+1);
            c.at(i-1) = A.at(i).at(j+2);
            j++;
        }
    }

    std::cout << "starting forward\n";

    // Forward subs. step
    std::vector<double> x = get_vector_x(N);
    std::vector<double> g(N);
    double h = 1.0/(N-1);
    for (int i = 1; i < N-1; i++){
        g.at(i) = (h*h)*f(x.at(i));
    }    

    std::vector<double> b_tilde(N);
    std::vector<double> g_tilde(N);

    b_tilde.at(0) = b.at(0);
    g_tilde.at(0) = g.at(0);
    for (int i = 1; i < N; i++){
        std::cout << "We got to " << i << "\n";
        std::cout << "We got to " << a.at(i) << "\n";
        b_tilde.at(i) = b.at(i) - (a.at(i-1) * c.at(i-1))/b_tilde.at(i-1);
        g_tilde.at(i) = g.at(i) - (a.at(i-1) * g_tilde.at(i-1))/b_tilde.at(i-1);
        
    }

    std::cout << "starting back\n";

    // Back subs. step
    std::vector<double> v(N);
    v.at(N-1) = g_tilde.at(N-1)/b_tilde.at(N-1);
    for (int i = N-2; i >= 0; i--){
        v.at(i) = (g_tilde.at(i) - c.at(i)*v.at(i+1))/b_tilde.at(i);
    }
    
    return v;
}

