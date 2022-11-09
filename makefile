n_body := 100
n_iterations := 100

checkpoint_folder := ./checkpoints/sequential_1000_20221107093003


default: run_seq

run_seq: seq
	./seq $(n_body) $(n_iterations)
run_video: video
	./video $(checkpoint_folder)

seq:
	g++ -g3 ./src/sequential.cpp -o seq -O2 -std=c++11
mpi:	
	mpic++ ./src/mpi.cpp -o mpi -std=c++11
pthread:
	g++ ./src/pthread.cpp -o pthread -lpthread -O2 -std=c++11
cuda:
	nvcc ./src/cuda.cu -o cuda -O2 --std=c++11
seqg:
	g++ ./src/sequential.cpp -o seqg -I/usr/include -L/usr/local/lib -L/usr/lib -lglut -lGLU -lGL -lm -DGUI -O2 -std=c++11
mpig:
	mpic++ ./src/mpi.cpp -o mpig -I/usr/include -L/usr/local/lib -L/usr/lib -lglut -lGLU -lGL -lm -DGUI -std=c++11
pthreadg:
	g++ ./src/pthread.cpp -o pthreadg -I/usr/include -L/usr/local/lib -L/usr/lib -lglut -lGLU -lGL -lm -lpthread -DGUI -O2 -std=c++11
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
