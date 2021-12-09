N=400000000

g++ threads.cpp -pthread
for NTH in $(seq 1 16); do
    ./a.out ${N} ${NTH} 2>/dev//null
done > data.txt
