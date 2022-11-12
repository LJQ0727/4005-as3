#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <pthread.h>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/checkpoint.h"

int n_thd; // number of threads

int n_body;
int n_iteration;

pthread_barrier_t barrier;  // barrier for synchronization



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


void update_position(double *x, double *y, double *vx, double *vy, int n, int start, int count) {
    // start: the start index that we will update
    // count: how many will update in this thread

    // update position for each body
    for (int i = start; i < n; i++)
    {
        double delta_x = vx[i] * dt;    // the estimated distance that the body will move in the x direction
        double delta_y = vy[i] * dt;    // the estimated distance that the body will move in the y direction

        // update the position
        x[i] += delta_x;
        y[i] += delta_y;

    }

    pthread_barrier_wait(&barrier);     // synchronize threads, make sure the positions will not change in the next part
    
    // check if the body will go out of bounds. If so, it will bounce back
    for (int i = 0; i < n; i++)
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

    pthread_barrier_wait(&barrier);     // sychronize threads, make sure that before next velocity update this is completed

}

void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int n, int start, int count) {
    // start: the start index that we will update
    // count: how many will update in this thread

    // calculate force and acceleration, update velocity
    for (int i = start; i < start+count; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j) continue;

            // for each pair of bodies
            // calculate distance
            double distance_x = x[i] - x[j];
            double distance_y = y[i] - y[j];
            double distance = sqrt(distance_x * distance_x + distance_y * distance_y);

                        
            // calculate force
            double force = gravity_const * m[i] * m[j] / (distance * distance + error); // the force scalar
            // calculate acceleration from j to i
            double acceleration_i = force / m[i];
            double acceleration_x_i = -acceleration_i * distance_x / distance;  // acceleration on i on x axis
            double acceleration_y_i = -acceleration_i * distance_y / distance;

            if (distance < radius2)
            {
                // if the distance is too small, we will reverse the velocity of the two bodies
                vx[i] = -vx[i];
                vy[i] = -vy[i];
                break;
            }

            // update velocity
            vx[i] += acceleration_x_i * dt;
            vy[i] += acceleration_y_i * dt;

        }
    }
        
}


typedef struct {
    //TODO: specify your arguments for threads
    int rank;
    int world_size;
    double *m;
    double *x;
    double *y;
    double *vx;
    double *vy;
    int n;  // n elements in the array
    //TODO END
} Args;


void* worker(void* args) {
    //TODO: procedure in each threads

    Args* my_arg = (Args*) args;
    // int a = my_arg->a;
    // int b = my_arg->b;
    int rank = my_arg->rank;
    int world_size = my_arg->world_size;
    int n = my_arg->n;

    int num_my_element = n / world_size;
    int start = rank * num_my_element;
    // for the last thread, there can be different elements than others
    if (rank == world_size-1) {
        num_my_element = n - (world_size-1) * num_my_element;
    }

    double *m = my_arg->m;
    double *x = my_arg->x;
    double *y = my_arg->y;
    double *vx = my_arg->vx;
    double *vy = my_arg->vy;


    for (int i = 0; i < n_iteration; i++)
    {
        update_velocity(m,x,y,vx,vy,n,start,num_my_element);    // update velocity of this portion of body
        update_position(x,y,vx,vy,n,start,num_my_element);      // update position of this postion of body
    }
    

    // TODO END


    pthread_exit(NULL); // let the master thread join
}


void master(){
    double* m = new double[n_body];
    double* x = new double[n_body];
    double* y = new double[n_body];
    double* vx = new double[n_body];
    double* vy = new double[n_body];

    generate_data(m, x, y, vx, vy, n_body);

    Logger l = Logger("pthread", n_body, bound_x, bound_y);

    // initialize threads and args
    pthread_t thds[n_thd];
    Args args[n_thd];
    pthread_barrier_init(&barrier, NULL, n_thd+1);
    for (int i = 0; i < n_thd; i++)
    {
        args[i].rank = i;
        args[i].world_size = n_thd;
        args[i].m = m;
        args[i].x = x;
        args[i].y = y;
        args[i].vx = vx;
        args[i].vy = vy;
        args[i].n = n_body;
        pthread_create(&thds[i], NULL, worker, (void*) &args[i]);
    }


    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        //TODO: assign jobs

        pthread_barrier_wait(&barrier);     // synchronize all threads
        pthread_barrier_wait(&barrier);     // synchronize all threads
        //TODO End

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


int main(int argc, char *argv[]) {
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);
    n_thd = atoi(argv[3]);

    #ifdef GUI
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Pthread");
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0, bound_x, 0, bound_y);
    #endif
    master();

	return 0;
}

