#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/checkpoint.h"


int n_body;
int n_iteration;

int n_omp_threads;


void generate_data(double *m, double *x,double *y,double *vx,double *vy, int n) {
    // TODO: Generate proper initial position and mass for better visualization
    for (int i = 0; i < n; i++) {
        m[i] = rand() % max_mass + 1.0f;
        x[i] = rand() % bound_x;
        y[i] = rand() % bound_y;
        vx[i] = 0.0f;
        vy[i] = 0.0f;
    }
}

void check_bounds(double *x, double *y, double *vx, double *vy, int n_body) {
    //check if the body will go out of bounds. If so, it will bounce back
    #pragma omp parallel for
    for (int i = 0; i < n_body; i++)
    {
        if (x[i] <= 0 || x[i] >= bound_x)
        {
            vx[i] = -vx[i];
        }
        if (y[i] <= 0 || y[i] >= bound_y)
        {
            vy[i] = -vy[i];
        }
    }

} 


void update_position(double *x, double *y, double *vx, double *vy, int idx) {
    //TODO: update position 

    // update position for the body at idx 
    double delta_x = vx[idx] * dt;    // the estimated distance that the body will move in the x direction
    double delta_y = vy[idx] * dt;    // the estimated distance that the body will move in the y direction

    // update the position
    x[idx] += delta_x;
    y[idx] += delta_y;


}

void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int idx, int n_body) {
    // calculate force and acceleration, update velocity
    for (int j = 0; j < n_body; j++)
    {
        if (j == idx) continue;

        // for each pair of bodies
        // calculate distance
        double distance_x = x[idx] - x[j];
        double distance_y = y[idx] - y[j];
        double distance = sqrt(distance_x * distance_x + distance_y * distance_y);

                    
        // calculate force
        double force = gravity_const * m[idx] * m[j] / (distance * distance + error); // the force scalar
        // calculate acceleration from j to i
        double acceleration_i = force / m[idx];
        double acceleration_x_i = -acceleration_i * distance_x / distance;  // acceleration on i on x axis
        double acceleration_y_i = -acceleration_i * distance_y / distance;

        if (distance < radius2)
        {
            // if the distance is too small, we will reverse the velocity
            // we give them a little push to prevent from colliding again
            vx[idx] = -vx[idx];
            vy[idx] = -vy[idx];
            break;
        }

        // update velocity
        vx[idx] += acceleration_x_i * dt;
        vy[idx] += acceleration_y_i * dt;

    }
    
        
}


void master() {
    double* m = new double[n_body];
    double* x = new double[n_body];
    double* y = new double[n_body];
    double* vx = new double[n_body];
    double* vy = new double[n_body];

    generate_data(m, x, y, vx, vy, n_body);

    Logger l = Logger("openmp", n_body, bound_x, bound_y);

    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        
        //TODO: choose better threads configuration
        omp_set_num_threads(n_omp_threads);
        #pragma omp parallel for
        for (int i = 0; i < n_body; i++) {
            update_velocity(m, x, y, vx, vy, i, n_body);
        }

        omp_set_num_threads(8);
        #pragma omp parallel for
        for (int i = 0; i < n_body; i++) {
            update_position(x, y, vx, vy, i);
        }

        check_bounds(x, y, vx, vy, n_body);

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

    delete[] m;
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;
    
}


int main(int argc, char *argv[]){
    
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);
    n_omp_threads = atoi(argv[3]);

    #ifdef GUI
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(500, 500);
    glutCreateWindow("N Body Simulation Sequential Implementation");
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    gluOrtho2D(0, bound_x, 0, bound_y);
    #endif
    master();

    printf("Student ID: 120090727\n"); // replace it with your student id
    printf("Name: Li Jiaqi\n"); // replace it with your name
    printf("Assignment 2: N Body Simulation OpenMP Implementation\n");
    
    return 0;

}


