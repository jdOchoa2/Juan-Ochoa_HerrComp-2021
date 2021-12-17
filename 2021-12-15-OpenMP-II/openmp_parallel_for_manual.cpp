#include <omp.h>
#include <iostream>
#include <cmath>
#include <chrono>

int main(int argc, char *argv[]) {
    const int N = 80000000;
    double *a = new double[N];
    auto start = std::chrono::steady_clock::now();
#pragma omp parallel
    {
        int thid = omp_get_thread_num();
        int nth = omp_get_num_threads();
        int localsize = N/nth;
        int iimin = thid*localsize;
        int iimax = iimin + localsize;
        for(int i = iimin; i < iimax; i++) {
            a[i] = 2*i*std::sin(std::sqrt(i/56.7)) +
                std::cos(std::pow(i*i, 0.3));
        }
    }
    
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now()-start;
    std::cout << elapsed.count() << "\n";
    a[1]+=1;
    delete [] a;
    return 0;
}
