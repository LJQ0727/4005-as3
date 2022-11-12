n_body := 1000
n_iterations := 5000

n_thds := 4
n_omp_threads := 4
n_mpi_nodes := 4

checkpoint_folder := ./pthread_1000_20221111092509/

C_FLAG := -O2

default: run_mpiopenmpg

run_seq: seq
	./seq $(n_body) $(n_iterations)
run_seqg: seqg
	./seqg $(n_body) $(n_iterations)
run_pthread: pthread
	./pthread $(n_body) $(n_iterations) $(n_thds)
run_pthreadg: pthreadg
	./pthreadg $(n_body) $(n_iterations) $(n_thds)
run_openmp: openmp
	./openmp $(n_body) $(n_iterations) $(n_omp_threads)
run_openmpg: openmpg
	./openmpg $(n_body) $(n_iterations) $(n_omp_threads)
run_cuda: cuda
	./cuda $(n_body) $(n_iterations)
run_cudag: cudag
	./cudag $(n_body) $(n_iterations)
run_mpi: mpi
	mpirun -np $(n_mpi_nodes) ./mpi $(n_body) $(n_iterations)
run_mpig: mpig
	mpirun -np $(n_mpi_nodes) ./mpig $(n_body) $(n_iterations)
run_mpiopenmp: mpiopenmp
	mpirun -np $(n_mpi_nodes) ./mpiopenmp $(n_body) $(n_iterations) $(n_omp_threads)
run_mpiopenmpg: mpiopenmpg
	mpirun -np $(n_mpi_nodes) ./mpiopenmpg $(n_body) $(n_iterations) $(n_omp_threads)


run_video: video
	./video $(checkpoint_folder)

seq:
	g++ $(C_FLAG) ./src/sequential.cpp -o seq -std=c++11
mpi:	
	mpic++ ./src/mpi.cpp -o mpi -std=c++11
pthread:
	g++ $(C_FLAG) ./src/pthread.cpp -o pthread -lpthread -std=c++11
cuda:
	nvcc ./src/cuda.cu -o cuda -O2 --std=c++11
seqg:
	g++ $(C_FLAG) ./src/sequential.cpp -o seqg -I/usr/include -L/usr/local/lib -L/usr/lib -lglut -lGLU -lGL -lm -DGUI -std=c++11
mpig:
	mpic++ ./src/mpi.cpp -o mpig -I/usr/include -L/usr/local/lib -L/usr/lib -lglut -lGLU -lGL -lm -DGUI -std=c++11
pthreadg:
	g++ $(C_FLAG) ./src/pthread.cpp -o pthreadg -I/usr/include -L/usr/local/lib -L/usr/lib -lglut -lGLU -lGL -lm -lpthread -DGUI -std=c++11
cudag:
	nvcc ./src/cuda.cu -o cudag -I/usr/include -L/usr/local/lib -L/usr/lib -lglut -lGLU -lGL -lm -O2 -DGUI --std=c++11
openmp:
	g++ $(C_FLAG) ./src/openmp.cpp -o openmp -fopenmp -std=c++11
openmpg:
	g++ $(C_FLAG) ./src/openmp.cpp -o openmpg -fopenmp -I/usr/include -L/usr/local/lib -L/usr/lib -lglut -lGLU -lGL -lm -DGUI -std=c++11
mpiopenmp:
	mpic++ ./src/mpi.cpp -o mpiopenmp -fopenmp -std=c++11
mpiopenmpg:
	mpic++ ./src/mpi.cpp -o mpiopenmpg -I/usr/include -fopenmp -L/usr/local/lib -L/usr/lib -lglut -lGLU -lGL -lm -DGUI -std=c++11

video:
	g++ ./src/video.cpp -o video -I/usr/include -L/usr/local/lib -L/usr/lib -lglut -lGLU -lGL -lm -DGUI -O2 -std=c++11
all:
	make seq
	make mpi
	make pthread
	make cuda
	make seqg
	make mpig
	make pthreadg
	make cudag
	make openmp
	make openmpg
	make mpiopenmp
	make mpiopenmpg
	make video
clean:
	rm -f seq mpi pthread seqg mpig pthreadg cuda cudag openmp openmpg video mpiopenmp mpiopenmpg
