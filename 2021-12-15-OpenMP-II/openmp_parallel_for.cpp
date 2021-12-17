#include <omp.h>
#include <iostream>
#include <cmath>
#include <chrono>

int main(int argc, char *argv[]) {
    const int N = 80000000;
    double *a = new double[N];

    auto start = std::chrono::steady_clock::now();
    
#pragma omp parallel for
    for(int i = 0; i < N; i++) {
        a[i] = 2*i*std::sin(std::sqrt(i/56.7)) +
            std::cos(std::pow(i*i, 0.3));
    }
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now()-start;
    std::cout << elapsed.count() << "\n";
    a[1]+=1;
    delete [] a;
    return 0;
}
