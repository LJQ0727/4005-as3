n_body := 1000
n_iterations := 5000

n_thds := 4

checkpoint_folder := ./checkpoints/sequential_100_20221109192429/

# change to -O2 in prod
C_FLAG := -g3

default: run_pthread

run_seq: seq
	./seq $(n_body) $(n_iterations)
run_seqg: seqg
	./seqg $(n_body) $(n_iterations)
run_pthread: pthread
	./pthread $(n_body) $(n_iterations) $(n_thds)
run_pthreadg: pthreadg
	./pthreadg $(n_body) $(n_iterations) $(n_thds)
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
	g++ ./src/openmp.cpp -o openmp -fopenmp -O2 -std=c++11
openmpg:
	g++ ./src/openmp.cpp -o openmpg -fopenmp -I/usr/include -L/usr/local/lib -L/usr/lib -lglut -lGLU -lGL -lm -O2 -DGUI -std=c++11
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
	make video
clean:
	rm -f seq mpi pthread seqg mpig pthreadg cuda cudag openmp openmpg video
