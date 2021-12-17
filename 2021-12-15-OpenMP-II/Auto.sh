for NTH in $(seq 1 16); do OMP_NUM_THREADS=${NTH} ./a.out; done > Auto.txt
