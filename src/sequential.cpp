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

    // make a copy so that we can refer to the original position when determining collision
    double x_ori[n];
    double y_ori[n];
    for (int i = 0; i < n; i++)
    {
        x_ori[i] = x[i];
        y_ori[i] = y[i];
    }

    // update position for each body
    for (int i = 0; i < n; i++)
    {
        double delta_x = vx[i] * dt;    // the estimated distance that the body will move in the x direction
        double delta_y = vy[i] * dt;    // the estimated distance that the body will move in the y direction

        // check if the body will go out of bounds. If so, it will be located at the boundary with a bounced velocity
        if (x[i] + delta_x <= 0) {
            x[i] = 0;
            vx[i] = -vx[i];
        }
        else if (x[i] + delta_x >= bound_x) {
            x[i] = bound_x;
            vx[i] = -vx[i];
        }
        if (y[i] + delta_y <= 0) {
            y[i] = 0;
            vy[i] = -vy[i];
        }
        else if (y[i] + delta_y >= bound_y) {
            y[i] = bound_y;
            vy[i] = -vy[i];
        }

        // update the position
        x[i] += delta_x;
        y[i] += delta_y;


        // check if the updated position will cause collision, if so, it will be in the updated position with a bounced velocity
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                continue;
            }
            double distance = pow(x[i] + delta_x - x[j], 2) + pow(y[i] + delta_y - y[j], 2);
            if (distance < radius2)
            {
                vx[i] = -vx[i];
                vy[i] = -vy[i];
                break;
            }
        }

    }
    

}

void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int n) {
    // calculate force and acceleration, update velocity
    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < n; j++)
        {
            // for each pair of bodies
            // calculate distance
            double distance_x = x[i] - x[j];
            double distance_y = y[i] - y[j];
            double distance = sqrt(distance_x * distance_x + distance_y * distance_y);

                        
            // calculate force
            double force = gravity_const * m[i] * m[j] / (distance * distance + error); // the force scalar
            // calculate acceleration from j to i
            double acceleration = force / m[i];
            double acceleration_x_jtoi = -acceleration * distance_x / distance;  // acceleration on i on x axis
            double acceleration_y_jtoi = -acceleration * distance_y / distance;
            // calculate acceleration from i to j
            double acceleration_x_itoj = -acceleration_x_jtoi;
            double acceleration_y_itoj = -acceleration_y_jtoi;

            // update velocity
            vx[i] -= acceleration_x_jtoi * dt;
            vy[i] -= acceleration_y_jtoi * dt;
            vx[j] -= acceleration_x_itoj * dt;
            vy[j] -= acceleration_y_itoj * dt;
        }
    }
        
}
    


void master() {
    double* m = new double[n_body]; // mass of each body
    double* x = new double[n_body]; // x location of each body
    double* y = new double[n_body]; // y location of each body
    double* vx = new double[n_body];    // initial x velocity of each body
    double* vy = new double[n_body];    // initial y velocity of each body

    generate_data(m, x, y, vx, vy, n_body);

    Logger l = Logger("sequential", n_body, bound_x, bound_y);

    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        update_velocity(m, x, y, vx, vy, n_body);
        update_position(x, y, vx, vy, n_body);

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

    delete m;
    delete x;
    delete y;
    delete vx;
    delete vy;
    
}


int main(int argc, char *argv[]){
    
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);

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
    printf("Assignment 2: N Body Simulation Sequential Implementation\n");
    
    return 0;

}


