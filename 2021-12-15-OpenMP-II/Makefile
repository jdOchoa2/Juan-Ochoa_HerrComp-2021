all: openmp_parallel_for_manual.cpp openmp_parallel_for.cpp speedup.py parallel.py
	g++ -fopenmp $<;\
	bash Man.sh;\
	g++ -fopenmp openmp_parallel_for.cpp;\
	bash Auto.sh;\
	python3 speedup.py;\
	python3 parallel.py
clean:
	rm -f *.x *.png
