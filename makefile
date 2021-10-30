hellomake: psort.cpp
	mpic++ -c psort.cpp -o psort.o -fopenmp
	ar cr libpsort.a psort.o

my:psort.cpp a4.cpp
	mpic++ -c psort.cpp -o psort.o -fopenmp
	ar cr libpsort.a psort.o
	mpic++ -c a4.cpp -o a4.o -fopenmp
	mpic++ a4.o libpsort.a -o a.out -fopenmp
	rm -rf input_dir/*.txt output_dir/*.txt
	# time mpirun -np 4 ./a.out
	#mpirun -np 6 ./a.out
	#cat output_dir/inp_0.txt output_dir/inp_1.txt output_dir/inp_2.txt output_dir/inp_3.txt
	#cat output_dir/out_0.txt output_dir/out_1.txt output_dir/out_2.txt output_dir/out_3.txt
