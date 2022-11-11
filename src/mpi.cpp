#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <mpi.h>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/checkpoint.h"


int n_body;
int n_iteration;


int my_rank;
int world_size;


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


void update_position(double *x, double *y, double *vx, double *vy, int n) {
    //TODO: update position 

    // update position for each body
    for (int i = 0; i < n; i++)
    {
        double delta_x = vx[i] * dt;    // the estimated distance that the body will move in the x direction
        double delta_y = vy[i] * dt;    // the estimated distance that the body will move in the y direction

        // update the position
        x[i] += delta_x;
        y[i] += delta_y;

    }
    
    //check if the body will go out of bounds. If so, it will bounce back
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

}

void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int n) {
    // calculate force and acceleration, update velocity
    for (int i = 0; i < n; i++)
    {
        for (int j = i+1; j < n; j++)
        {
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
            // calculate acceleration from i to j
            double acceleration_j = force / m[j];
            double acceleration_x_j = acceleration_j * distance_x / distance;  // acceleration on j on x axis
            double acceleration_y_j = acceleration_j * distance_y / distance;

            if (distance < radius2)
            {
                // if the distance is too small, we will reverse the velocity of the two bodies
                // we give them a little push to prevent from colliding again
                vx[i] = -vx[i];
                vy[i] = -vy[i];
                vx[j] = -vx[j];
                vy[j] = -vy[j];
                break;
            }

            // update velocity
            vx[i] += acceleration_x_i * dt;
            vy[i] += acceleration_y_i * dt;
            vx[j] += acceleration_x_j * dt;
            vy[j] += acceleration_y_j * dt;

        }
    }
        
}


void slave(){
    // TODO: MPI routine
    double* local_m;
    double* local_x;
    double* local_y;
    double* local_vx;
    double* local_vy;
    // TODO End
}



void master() {
    double* total_m = new double[n_body];
    double* total_x = new double[n_body];
    double* total_y = new double[n_body];
    double* total_vx = new double[n_body];
    double* total_vy = new double[n_body];

    int num_elements = n_body;
    // use `world_size` and `my_rank` to determine the number of elements for each process
    int num_my_element = num_elements / world_size;

    // create local receive buffer
    double total_m_buf[num_my_element];
    double total_x_buf[num_my_element];
    double total_y_buf[num_my_element];
    double total_vx_buf[num_my_element];
    double total_vy_buf[num_my_element];

    generate_data(total_m, total_x, total_y, total_vx, total_vy, n_body);

    Logger l = Logger("mpi", n_body, bound_x, bound_y);

    // use MPI_scatterv to scatter the data
    int counts_send[world_size];
    for (int i = 0; i < world_size-1; i++) {
        counts_send[i] = num_my_element;
    }
    counts_send[world_size-1] = num_elements - (world_size-1) * num_my_element;

    int displacements[world_size];
    for (int i = 0; i < world_size; i++) {
        displacements[i] = i * num_my_element;
    }
    MPI_Scatterv(total_m, counts_send, displacements, MPI_DOUBLE, total_m_buf, num_my_element, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // scatter the rest of the data
    MPI_Scatterv(total_x, counts_send, displacements, MPI_DOUBLE, total_x_buf, num_my_element, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(total_y, counts_send, displacements, MPI_DOUBLE, total_y_buf, num_my_element, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(total_vx, counts_send, displacements, MPI_DOUBLE, total_vx_buf, num_my_element, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(total_vy, counts_send, displacements, MPI_DOUBLE, total_vy_buf, num_my_element, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        // TODO: MPI routine
        
        // TODO End
        l.save_frame(total_x, total_y);

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
            xi = total_x[i];
            yi = total_y[i];
            glVertex2f(xi, yi);
        }
        glEnd();
        glFlush();
        glutSwapBuffers();
        #else

        #endif
    }

    delete total_m;
    delete total_x;
    delete total_y;
    delete total_vx;
    delete total_vy;

}




int main(int argc, char *argv[]) {
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	if (my_rank == 0) {
		#ifdef GUI
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize(500, 500); 
		glutInitWindowPosition(0, 0);
		glutCreateWindow("N Body Simulation MPI Implementation");
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glMatrixMode(GL_PROJECTION);
		gluOrtho2D(0, bound_x, 0, bound_y);
		#endif
        master();
	} else {
        slave();
    }

	if (my_rank == 0){
		printf("Student ID: 120090727\n"); // replace it with your student id
		printf("Name: Li Jiaqi\n"); // replace it with your name
		printf("Assignment 2: N Body Simulation MPI Implementation\n");
	}

	MPI_Finalize();

	return 0;
}

