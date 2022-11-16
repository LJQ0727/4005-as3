#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <mpi.h>
#include <omp.h>

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

int n_omp_threads;


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
        if (force > 35000) force = 35000;
        // calculate acceleration from j to i
        double acceleration_i = force / m[idx];
        double acceleration_x_i = -acceleration_i * distance_x / distance;  // acceleration on i on x axis
        double acceleration_y_i = -acceleration_i * distance_y / distance;

        if (distance*distance < radius2)
        {
            // if the distance is too small, we will reverse the velocity
            vx[idx] = -vx[idx];
            vy[idx] = -vy[idx];
            break;
        }

        // update velocity
        vx[idx] += acceleration_x_i * dt;
        vy[idx] += acceleration_y_i * dt;

    }
    
        
}

void check_bounds(double *x, double *y, double *vx, double *vy, int n_body, int idx, int num_my_elements) {
    //check if the body will go out of bounds. If so, it will bounce back
    omp_set_num_threads(n_omp_threads);
    #pragma omp parallel for
    for (int i = idx; i < idx + num_my_elements; i++)
    {
        if (x[i] <= 2000 || x[i] >= 3000)
        {
            vx[i] = -vx[i];
        }
        if (y[i] <= 2000 || y[i] >= 3000)
        {
            vy[i] = -vy[i];
        }
    }

} 


void slave(){
    int num_my_element = n_body / (world_size - 1);
    int num_elements = num_my_element * world_size;

    double* total_m = new double[num_elements];
    double* total_x = new double[num_elements];
    double* total_vy = new double[num_elements];
    double* total_y = new double[num_elements];
    double* total_vx = new double[num_elements];
    

    // receive the data and store in the total array
    MPI_Bcast(total_m, num_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(total_x, num_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(total_y, num_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(total_vx, num_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(total_vy, num_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    for (int i = 0; i < n_iteration; i++)
    {
        omp_set_num_threads(n_omp_threads);
        #pragma omp parallel for
        for (int j = my_rank * num_my_element; j < my_rank * num_my_element + num_my_element; j++)
        {
            update_velocity(total_m, total_x, total_y, total_vx, total_vy, j, n_body);
            update_position(total_x, total_y, total_vx, total_vy, j);
        }
        check_bounds(total_x, total_y, total_vx, total_vy, num_elements, my_rank * num_my_element, num_my_element);
    

        MPI_Allgather(MPI_IN_PLACE, num_my_element, MPI_DOUBLE, total_x, num_my_element, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgather(MPI_IN_PLACE, num_my_element, MPI_DOUBLE, total_y, num_my_element, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgather(MPI_IN_PLACE, num_my_element, MPI_DOUBLE, total_vx, num_my_element, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgather(MPI_IN_PLACE, num_my_element, MPI_DOUBLE, total_vy, num_my_element, MPI_DOUBLE, MPI_COMM_WORLD);        
    }
    delete[] total_m;
    delete[] total_x;
    delete[] total_y;
    delete[] total_vx;
    delete[] total_vy;
}



void master() {
    int num_my_element;
    if (world_size == 1) {
        num_my_element = n_body;
    } else {
        num_my_element = n_body / (world_size - 1);
    }
    int num_elements = num_my_element * world_size;

    double* total_m = new double[num_elements];
    double* total_x = new double[num_elements];
    double* total_vy = new double[num_elements];
    double* total_y = new double[num_elements];
    double* total_vx = new double[num_elements];

    for (int i = 0; i < num_elements; i++)
    {
        total_m[i] = 0;
    }
    
    // initialize data
    generate_data(total_m, total_x, total_y, total_vx, total_vy, n_body);

    Logger l = Logger("mpiopenmp", n_body, bound_x, bound_y);

    // copy the data to each other node
    MPI_Bcast(total_m, num_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(total_x, num_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(total_y, num_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(total_vx, num_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(total_vy, num_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        omp_set_num_threads(n_omp_threads);
        #pragma omp parallel for
        for (int j = my_rank * num_my_element; j < my_rank * num_my_element + num_my_element; j++)
        {
            update_velocity(total_m, total_x, total_y, total_vx, total_vy, j, n_body);
            update_position(total_x, total_y, total_vx, total_vy, j);
        }
        check_bounds(total_x, total_y, total_vx, total_vy, num_elements, my_rank * num_my_element, num_my_element);
    
        MPI_Allgather(MPI_IN_PLACE, num_my_element, MPI_DOUBLE, total_x, num_my_element, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgather(MPI_IN_PLACE, num_my_element, MPI_DOUBLE, total_y, num_my_element, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgather(MPI_IN_PLACE, num_my_element, MPI_DOUBLE, total_vx, num_my_element, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgather(MPI_IN_PLACE, num_my_element, MPI_DOUBLE, total_vy, num_my_element, MPI_DOUBLE, MPI_COMM_WORLD);        
        
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

    delete[] total_m;
    delete[] total_x;
    delete[] total_y;
    delete[] total_vx;
    delete[] total_vy;

}




int main(int argc, char *argv[]) {
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);
    n_omp_threads = atoi(argv[3]);

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	if (my_rank == 0) {
		#ifdef GUI
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize(800, 800); 
		glutInitWindowPosition(0, 0);
		glutCreateWindow("N Body Simulation hybrid MPI OpenMP Implementation");
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
		printf("Assignment 2: N Body Simulation MPI+OpenMP Implementation\n");
	}

	MPI_Finalize();

	return 0;
}

