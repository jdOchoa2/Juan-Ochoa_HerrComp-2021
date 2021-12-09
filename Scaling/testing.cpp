#include <iostream>
#include <cstdlib>
#include <numeric>
#include <algorithm>


int main(int argc, char **argv){
    int N = std::atoi(argv[1]);
    int NTH = std::atoi(argv[2]);

    double localsize =double(N)/NTH;
    std::cout<<localsize;
}
