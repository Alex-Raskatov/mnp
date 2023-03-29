#include <iostream>
#include <cstdio>
#include <fstream>
#include <type_traits>
#include <cmath>

const int K = 100;
const int T = 100;
const int step = 10;
const double gamm = 0.5;

double initial_f (int k) {
    return std::exp( -((k - K/5.0)*(k - K/5.0))/(2.0*K));
}

void set_initial_values (double* data) {

    for (int k = 0; k < K; k++) {
        data[k] = initial_f(k);
    }
}

void make_iteration( double* current, double* next) {

    next[0] = current[0];
    for (int k = 1; k < K; k++) {
        next[k] = current[k] - gamm * (current[k] - current[k-1]);
    }
}

void save_to_file(double* data, int i) {

    char filename[30];
    sprintf(filename, "./data/task1/out_%03d.dat", i);
    std::ofstream file(filename);
    for (int k = 0; k < K; k++) {
        file << data[k] << std::endl;
    }
    file.close();
}


int main () {

    double* current = new double[K];
    double* next = new double[K];

    set_initial_values(current);

    for (int i = 0; i < T; i++) {

        if (i % step == 0) {
            save_to_file(current, i);
        }

        make_iteration(current, next);
        std::swap(current, next);
    }

    delete [] current;
    delete [] next;

    return 0;
}