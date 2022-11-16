#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

#include <assert.h>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/checkpoint.h"


int block_size = 1024;


int n_body;
int n_iteration;

__device__ __managed__ int bound_x_d = 4000;
__device__ __managed__ int bound_y_d = 4000;
__device__ __managed__ int max_mass_d = 400;
__device__ __managed__ double error_d = 1e-5f;
__device__ __managed__ double dt_d = 0.0001f;
__device__ __managed__ double gravity_const_d = 1000000.0f;
__device__ __managed__ double radius2_d = 0.01f;


__global__ void update_position(double *x, double *y, double *vx, double *vy, int n) {
    // update position 
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) {
        x[i] += vx[i] * dt_d;
        y[i] += vy[i] * dt_d;
    }
}

__global__ void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int n) {
    //TODO: calculate force and acceleration, update velocity
    int i = blockDim.x * blockIdx.x + threadIdx.x;  // update the ith element in the array
    if (i < n) {
        for (int j = 0; j < n; j++)
        {
            if (i == j) continue;

            // for each pair of bodies
            // calculate distance
            double distance_x = x[i] - x[j];
            double distance_y = y[i] - y[j];
            double distance = sqrt(distance_x * distance_x + distance_y * distance_y);

            
            // calculate force
            double force = gravity_const_d * m[i] * m[j] / (distance * distance + error_d); // the force scalar
            if (force > 35000) force = 35000;
            // calculate acceleration from j to i
            double acceleration_i = force / m[i];
            double acceleration_x_i = -acceleration_i * distance_x / distance;  // acceleration on i on x axis
            double acceleration_y_i = -acceleration_i * distance_y / distance;

            if (distance*distance < radius2_d)
            {
                // if the distance is too small, we will reverse the velocity of the two bodies
                vx[i] = -vx[i];
                vy[i] = -vy[i];
                break;
            }

            // update velocity
            vx[i] += acceleration_x_i * dt_d;
            vy[i] += acceleration_y_i * dt_d;
        }
        
    }
}

__global__ void check_bounds(double *x, double *y, double *vx, double *vy, int n_body) {
    //check if the body will go out of bounds. If so, it will bounce back
    for (int i = 0; i < n_body; i++)
    {
        if (x[i] <= 0 || x[i] >= bound_x_d)
        {
            vx[i] = -vx[i];
        }
        if (y[i] <= 0 || y[i] >= bound_y_d)
        {
            vy[i] = -vy[i];
        }
    }

} 



void generate_data(double *m, double *x,double *y,double *vx,double *vy, int n) {
    // TODO: Generate proper initial position and mass for better visualization
    srand((unsigned)time(NULL));
    for (int i = 0; i < n; i++) {
        m[i] = rand() % max_mass + 1.0f;
        x[i] = 2000.0f + rand() % (bound_x / 4);
        y[i] = 2000.0f + rand() % (bound_y / 4);
        vx[i] = 0.0f;
        vy[i] = 0.0f;
    }
}


void master() {
    double* m = new double[n_body];
    double* x = new double[n_body];
    double* y = new double[n_body];
    double* vx = new double[n_body];
    double* vy = new double[n_body];

    generate_data(m, x, y, vx, vy, n_body);

    Logger l = Logger("cuda", n_body, bound_x, bound_y);

    double *device_m;
    double *device_x;
    double *device_y;
    double *device_vx;
    double *device_vy;


    cudaMalloc(&device_m, n_body * sizeof(double));
    cudaMalloc(&device_x, n_body * sizeof(double));
    cudaMalloc(&device_y, n_body * sizeof(double));
    cudaMalloc(&device_vx, n_body * sizeof(double));
    cudaMalloc(&device_vy, n_body * sizeof(double));

    cudaMemcpy(device_m, m, n_body * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(device_x, x, n_body * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(device_y, y, n_body * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(device_vx, vx, n_body * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(device_vy, vy, n_body * sizeof(double), cudaMemcpyHostToDevice);

    int n_block = n_body / block_size + 1; 

    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();


        update_velocity<<<n_block, block_size>>>(device_m, device_x, device_y, device_vx, device_vy, n_body);
        update_position<<<n_block, block_size>>>(device_x, device_y, device_vx, device_vy, n_body);
        check_bounds<<<1,1>>>(device_x, device_y, device_vx, device_vy, n_body);


        cudaMemcpy(x, device_x, n_body * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(y, device_y, n_body * sizeof(double), cudaMemcpyDeviceToHost);

        l.save_frame(x, y);

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = t2 - t1;
        printf("Iteration %d, elapsed time: %.3f\n", i, time_span);

        #ifdef GUI
        glClear(GL_COLOR_BUFFER_BIT);
        glColor3f(1.0f, 0.0f, 0.0f);
        glPointSize(2.0f);
        glBegin(GL_POINTS);
        double xi;
        double yi;
        for (int i = 0; i < n_body; i++){
            xi = x[i];
            yi = y[i];
            glVertex2f(xi, yi);
        }
        glEnd();
        glFlush();
        glutSwapBuffers();
        #else

        #endif

    }

    cudaFree(device_m);
    cudaFree(device_x);
    cudaFree(device_y);
    cudaFree(device_vx);
    cudaFree(device_vy);

    cudaFree(device_m);
    cudaFree(device_x);
    cudaFree(device_y);
    cudaFree(device_vx);
    cudaFree(device_vy);

    delete[] m;
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;
    
}


int main(int argc, char *argv[]){
    
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);

    #ifdef GUI
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(500, 500);
    glutCreateWindow("N Body Simulation CUDA Implementation");
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    gluOrtho2D(0, bound_x, 0, bound_y);
    #endif

    master();

    printf("Student ID: 120090727\n"); // replace it with your student id
    printf("Name: Li Jiaqi\n"); // replace it with your name
    printf("Assignment 2: N Body Simulation CUDA Implementation\n");

    return 0;

}


