//
//  main.c
//  RollerCoaster
//
//  Created by Ryan MacLeod on 2012-11-09.
//  Copyright (c) 2012 Ryan MacLeod. All rights reserved.
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <GLUT/glut.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define RAD2DEG 180.0/M_PI
#define DEG2RAD M_PI/180.0

// Define macro for the number of segments per curve.
#define NUMBER_SEGMENTS 100

// The radius of the cylinger that makes up the skybox.
#define SKYBOX_RADIUS 200

// The maximum number of data points the roller coaster can have.
#define COASTER_POINTS 100

// Display Callbacks for the Three Screens
static void myTimer(int value);
static void	myDisplay(void);
static void	myReshape(int w, int h);
static void keyPress(int key, int x, int y);
static void keyRelease(int key, int x, int y);
static void myKey(unsigned char key, int x, int y);
// Draw functions for the different objects in the scene.
static void drawCurve(void);
static void drawGround(void);
static void drawSkybox(void);
static void drawConnectors();
static void drawColumns();

// Initialize functions for roller coaster.
static void	init(void);

static double xMax, yMax;
static double phi = 0;
static double ctrlpoints[COASTER_POINTS][3] = {
    {15.0, 10.0, 20.0}, {0.0, 15.0, 15.0}, {-10.0, 15.0, 10.0}, {-25.0, 10.0, 5.0}, {-30.0, -5.0, 0.0},
    {-30.0, -20.0, 0.0}, {-40.0, -30.0, 5.0}, {-55.0, -30.0, 10.0}, {-65.0, -20.0, 10.0}, {-70.0, -5.0, 5.0},
    {-70.0, 5.0, 30.0}, {-70.0, 15.0, 10.0}, {-65.0, 30.0, 0.0}, {-55.0, 40.0, -15.0}, {-45.0, 40.0, -20.0},
    {-35.0, 40.0, -25.0}, {-25.0, 40.0, -20.0}, {-15.0, 40.0, -15.0}, {-5.0, 40.0, -10.0}, {5.0, 40.0, -5.0},
    {15.0, 40.0, -2.5}, {25.0, 35.0, 0.0}, {30.0, 25.0, 2.5}, {40.0, 20.0, 5.0}, {50.0, 25.0, 7.5},
    {55.0, 35.0, 10.0}, {50.0, 45.0, 12.5}, {40.0, 50.0, 15.0}, {30.0, 45.0, 17.5}, {25.0, 35.0, 20.0},
    {30.0, 25.0, 22.5}, {40.0, 20.0, 25.0}, {50.0, 25.0, 27.5}, {55.0, 35.0, 30.0}, {50.0, 45.0, 32.5},
    {40.0, 50.0, 35.0}, {30.0, 45.0, 37.5}, {25.0, 35.0, 40.0}, {30.0, 25.0, 42.5}, {40.0, 20.0, 45.0},
    {50.0, 25.0, 47.5}, {55.0, 35.0, 50.0}, {50.0, 45.0, 52.5}, {40.0, 50.0, 55.0}, {30.0, 45.0, 57.5},
    {25.0, 35.0, 60.0}, {30.0, 25.0, 60.0}, {40.0, 20.0, 60.0}, {50.0, 20.0, 60.0}, {60.0, 20.0, 60.0},
    {70.0, 20.0, 0.0}, {80.0, 12.5, 0.0}, {80.0, 2.5, 0.0}, {70.0, -5.0, 0.0}, {60.0, -5.0, 0.0},
    {50.0, -5.0, 10.0}, {40.0, -10.0, 20.0}, {30.0, -15.0, 20.0}, {20.0, -25.0, 10.0}, {20.0, -35.0, 0.0},
    {30.0, -40.0, -10.0}, {40.0, -40.0, -15.0}, {50.0, -35.0, -20.0}, {50.0, -25, -15.0}, {40.0, -15.0, -10.0},
    {30.0, -10.0, 0.0}, {20.0, -5.0, 10.0}, {15.0, 10.0, 20.0}, {0.0, 15.0, 15.0}, {-10.0, 15.0, 10.0}};
static double** leftRail;
static double** rightRail;
static double** centerRail;
static double** columnTopRight;
static double** columnTopLeft;
static double** qValues;
static double** dqValues;
static double** ddqValues;
static double** uValues;
static double** vValues;
static double** nValues;
// The number of control points in this data set.
static int numberPoints = 70;
// Used to choose whether to use the rotating camera or the on coaster camera.
static int cameraPosition = 1;
// The current position of the roller coaster camera.
static int currentPosition = 0;
// The current velocity of the roller coaster camera.
static int velocity = 0;

//Display list that will hold the points and the drawing.
static GLuint scene;

/* Holds the texture information for the ground which I didn't end up using because I couldn't
 get them looking very nice. The colours from the ground interferred with the sky so I left them 
 out.
static int textureWidth, textureHeight;
static GLuint texture1, texture2;
static GLubyte  *texture;*/ 

// Helper functions to determine cross-products...
static void calculateVectors(void);
static double *crossProduct(double *A, double *B);
static double *negationProduct(double *A);
static double *normalizeProduct(double *A);
static double *q2D(double u, double *P0, double *P1, double *P2, double *P3);
static double *q(double u, double *P0, double *P1, double *P2, double *P3);
static double *dq(double u, double *P0, double *P1, double *P2, double *P3);
static double *ddq(double u, double *P0, double *P1, double *P2, double *P3);

// Set the position of the right, left, center rails, as well as the top of the two support columns.
static double *setLeftRail(double *P0, double *P1);
static double *setRightRail(double *P0, double *P1);
static double *setCenterRail(double *P0);
double *setColumnTopRight(double *LEFT, double *CENTER);
double *setColumnTopLeft(double *RIGHT, double *CENTER);

/* Run through the main section of code for the program, this will simulate a roller coaster and change
 cameras based on user input.*/
int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1000, 600);
    glutInitWindowPosition(0, 0);
    glutCreateWindow("Roller Coaster");
    
    glutSpecialFunc(keyPress);
    glutSpecialUpFunc(keyRelease);
    glutKeyboardFunc(myKey);
    
    // Set up the data points to use with the coaster and the display lists.
    init();
    
    glutDisplayFunc(myDisplay);
    glutReshapeFunc(myReshape);
    glutTimerFunc(33, myTimer, 0);
    
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    
    
    
    glutMainLoop();
    
    return 0;
}

/* Initialize anything necessary to set up the scene for the roller coaster simulation. */
void init(void){
    glShadeModel(GL_SMOOTH);
    glEnable(GL_DEPTH_TEST);
    
    // Read in the control points from a file, first lets test without that feature.
    leftRail = (double **)malloc(COASTER_POINTS*NUMBER_SEGMENTS * sizeof(double *));
    rightRail = (double **)malloc(COASTER_POINTS*NUMBER_SEGMENTS * sizeof(double *));
    centerRail = (double **)malloc(COASTER_POINTS*NUMBER_SEGMENTS * sizeof(double *));
    columnTopRight = (double **)malloc(COASTER_POINTS*NUMBER_SEGMENTS * sizeof(double *));
    columnTopLeft = (double **)malloc(COASTER_POINTS*NUMBER_SEGMENTS * sizeof(double *));
    qValues = (double **)malloc(COASTER_POINTS*NUMBER_SEGMENTS * sizeof(double *));
    dqValues = (double **)malloc(COASTER_POINTS*NUMBER_SEGMENTS * sizeof(double *));
    ddqValues = (double **)malloc(COASTER_POINTS*NUMBER_SEGMENTS * sizeof(double *));
    uValues = (double **)malloc(COASTER_POINTS*NUMBER_SEGMENTS * sizeof(double *));
    vValues = (double **)malloc(COASTER_POINTS*NUMBER_SEGMENTS * sizeof(double *));
    nValues = (double **)malloc(COASTER_POINTS*NUMBER_SEGMENTS * sizeof(double *));
    calculateVectors();
    
    // Generate a display list that will hold the scene.
    scene = glGenLists(1);
    glNewList(scene, GL_COMPILE);
        // Draw the ground and colour it green.
        drawGround();
    
        // Draw the sky and colour it blue.
        drawSkybox();
        // Draw the coaster.
        drawCurve();
    
        // Draw the connection pieces for the rails.
        drawConnectors();
    
        // Draw the columnst that support the rails.
        drawColumns();
    glEndList();
}

// Calculate the vertices of the curve as well as the normal, up, and tangent vectors.
void calculateVectors(){
    double up[3] = {0.0, 0.0, 1.0};
    int index = 0;
    
    for (int j = 0; j < numberPoints-3; j++){
        for (int i = 0; i <= NUMBER_SEGMENTS; i++){
            qValues[(j*NUMBER_SEGMENTS + i)] = q((double)i/NUMBER_SEGMENTS, ctrlpoints[j], ctrlpoints[j+1], ctrlpoints[j+2], ctrlpoints[j+3]);
            dqValues[(j*NUMBER_SEGMENTS + i)] = dq((double)i/NUMBER_SEGMENTS, ctrlpoints[j], ctrlpoints[j+1], ctrlpoints[j+2], ctrlpoints[j+3]);
            ddqValues[(j*NUMBER_SEGMENTS + i)] = ddq((double)i/NUMBER_SEGMENTS, ctrlpoints[j], ctrlpoints[j+1], ctrlpoints[j+2], ctrlpoints[j+3]);
            nValues[(j*NUMBER_SEGMENTS + i)] = normalizeProduct(negationProduct(dqValues[(j*NUMBER_SEGMENTS + i)]));
            uValues[(j*NUMBER_SEGMENTS + i)] = normalizeProduct(crossProduct(up, nValues[(j*NUMBER_SEGMENTS + i)]));
            vValues[(j*NUMBER_SEGMENTS + i)] = normalizeProduct(crossProduct(nValues[(j*NUMBER_SEGMENTS + i)], uValues[(j*NUMBER_SEGMENTS + i)]));
        }
    }
    for (int j = 0; j < numberPoints-3; j++){
        for (int i = 0; i <= NUMBER_SEGMENTS; i++){
            index = j*NUMBER_SEGMENTS + i;
            // Set the three rails.
            leftRail[index] = setLeftRail(qValues[index], uValues[index]);
            rightRail[index] = setRightRail(qValues[index], uValues[index]);
            centerRail[index] = setCenterRail(qValues[index]);
            columnTopRight[index] = setColumnTopRight(leftRail[index], centerRail[index]);
            columnTopLeft[index] = setColumnTopRight(rightRail[index], centerRail[index]);
        }
    }
    
    
}

// Solve for a point along a curve segment.
double *q(double u, double *P0, double *P1, double *P2, double *P3){
    double *C = malloc(3*sizeof(double));
    if(C == NULL){
		fprintf(stderr, "out of memory\n");
		exit(EXIT_FAILURE);
    }
    
    // Initialize the values for u, u squared, and u cubed. Limits the number of multiplications.
    double c = u;
    double b = c * c;
    double a = b * c;
    
    // Calculate the r values for the equation based on the u value.
    double r0 = (1.0/6.0)*a;
    double r1 = ((-1.0/2.0)*a) + ((1.0/2.0)*b) + ((1.0/2.0)*c) + (1.0/6.0);
    double r2 = ((1.0/2.0)*a) - b + (2.0/3.0);
    double r3 = ((-1.0/6.0)*a) + ((1.0/2.0)*b) + ((-1.0/2.0)*c) + (1.0/6.0);
    
    C[0] = P0[0]*r3 + P1[0]*r2 + P2[0]*r1 + P3[0]*r0;
    C[1] = P0[1]*r3 + P1[1]*r2 + P2[1]*r1 + P3[1]*r0;
    C[2] = P0[2]*r3 + P1[2]*r2 + P2[2]*r1 + P3[2]*r0;
    
    return C;
}
// Solve for a derivative of a point along the curve segment.
double *dq(double u, double *P0, double *P1, double *P2, double *P3){
    double *C = malloc(3*sizeof(double));
    if(C == NULL){
		fprintf(stderr, "out of memory\n");
		exit(EXIT_FAILURE);
    }
    

    // Initialize the values for u, u squared, and u cubed. Limits the number of multiplications.
    double c = u;
    double b = c * c;
    
    // Calculate the r values for the equation based on the u value.
    double r0 = (1.0/2.0)*b;
    double r1 = ((-3.0/2.0)*b) + (c) + ((1.0/2.0));
    double r2 = ((3.0/2.0)*b) - 2*c;
    double r3 = ((-1.0/2.0)*b) + (c) + ((-1.0/2.0));
    
    C[0] = P0[0]*r3 + P1[0]*r2 + P2[0]*r1 + P3[0]*r0;
    C[1] = P0[1]*r3 + P1[1]*r2 + P2[1]*r1 + P3[1]*r0;
    C[2] = P0[2]*r3 + P1[2]*r2 + P2[2]*r1 + P3[2]*r0;
    
    return C;
}
// Solve for the second derivative of a point along the curve.
double *ddq(double u, double *P0, double *P1, double *P2, double *P3){
    double *C = malloc(3*sizeof(double));
    if(C == NULL){
		fprintf(stderr, "out of memory\n");
		exit(EXIT_FAILURE);
    }
    
    // Initialize the values for u, u squared, and u cubed. Limits the numberf multiplications.
    double c = u;
    
    // Calculate the r values for the equation based on the u value.
    double r0 = c;
    double r1 = (-3.0*c) + 1;
    double r2 = ((3.0)*c) - 2;
    double r3 = (-c) + 1;
    
    C[0] = P0[0]*r3 + P1[0]*r2 + P2[0]*r1 + P3[0]*r0;
    C[1] = P0[1]*r3 + P1[1]*r2 + P2[1]*r1 + P3[1]*r0;
    C[2] = P0[2]*r3 + P1[2]*r2 + P2[2]*r1 + P3[2]*r0;
    
    
    return C;
}

// Negate the values of a vector and return the new vector.
double *negationProduct(double *A){
    double *B = malloc(3*sizeof(double));
    if(B == NULL){
		fprintf(stderr, "out of memory\n");
		exit(EXIT_FAILURE);
    }
    
    B[0] = -1.0*A[0];
    B[1] = -1.0*A[1];
    B[2] = -1.0*A[2];
    return B;
}
// Normalize a vector and return the normalized vector.
double *normalizeProduct(double *A){
    double *B = malloc(3*sizeof(double));
    if(B == NULL){
		fprintf(stderr, "out of memory\n");
		exit(EXIT_FAILURE);
    }
    
    double divisor = sqrt(pow(A[0],2) + pow(A[1],2) + pow(A[2],2));
    if(divisor > 0){
        B[0] = A[0]/divisor;
        B[1] = A[1]/divisor;
        B[2] = A[2]/divisor;
    }
    else{
        B[0] = A[0];
        B[1] = A[1];
        B[2] = A[2];
    }
    return B;
}

// Find the cross products of two vectors.
double *crossProduct(double *A, double *B){
    double *C = malloc(3*sizeof(double));
    if(C == NULL){
		fprintf(stderr, "out of memory\n");
		exit(EXIT_FAILURE);
    }
    
    C[0] = A[1]*B[2] - A[2]*B[1];
    C[1] = A[2]*B[0] - A[0]*B[2];
    C[2] = A[0]*B[1] - A[1]*B[0];
    return C;
}


// Find the correct position of the left rail based on the direction of the coaster.
double *setLeftRail(double *P1, double *P2){
    double *A = malloc(3*sizeof(double));
    if(A == NULL){
		fprintf(stderr, "out of memory\n");
		exit(EXIT_FAILURE);
    }
    
    A[0] = P1[0] - P2[0];
    A[1] = P1[1] - P2[1];
    A[2] = P1[2];
    
    return A;
}

// Find the correct position of the right rail based on the direction of the coaster.
double *setRightRail(double *P1, double *P2){
    double *A = malloc(3*sizeof(double));
    if(A == NULL){
		fprintf(stderr, "out of memory\n");
		exit(EXIT_FAILURE);
    }
    
    A[0] = P1[0] + P2[0];
    A[1] = P1[1] + P2[1];
    A[2] = P1[2];
    return A;
}

// Set the center rail one point below the curve.
double *setCenterRail(double *POINT){
    double *A = malloc(3*sizeof(double));
    if(A == NULL){
		fprintf(stderr, "out of memory\n");
		exit(EXIT_FAILURE);
    }
    
    A[0] = POINT[0];
    A[1] = POINT[1];
    A[2] = POINT[2] - 1;
    
    return A;
}

// Get the postion for the right support column.
double *setColumnTopRight(double *LEFT, double *CENTER){
    double *A = malloc(3*sizeof(double));
    if(A == NULL){
		fprintf(stderr, "out of memory\n");
		exit(EXIT_FAILURE);
    }
    
    A[0] = CENTER[0] + 4*(CENTER[0] - LEFT[0]);
    A[1] = CENTER[1] + 4*(CENTER[1] - LEFT[1]);
    A[2] = CENTER[2] + 2*(CENTER[2] - LEFT[2]);
    
    return A;
}

// Get the position for the left support column.
double *setColumnTopLeft(double *RIGHT, double *CENTER){
    double *A = malloc(3*sizeof(double));
    if(A == NULL){
		fprintf(stderr, "out of memory\n");
		exit(EXIT_FAILURE);
    }
    
    A[0] = CENTER[0] + 4*(CENTER[0] - RIGHT[0]);
    A[1] = CENTER[1] + 4*(CENTER[1] - RIGHT[1]);
    A[2] = CENTER[2] + 2*(CENTER[2] - RIGHT[2]);
    
    return A;
}

/* Draw large planes that will act as the ground and sky to give the roller coaster a feeling of place
 in the 3d world. */
void drawGround(){
    //glBindTexture(GL_TEXTURE_2D, texture1);
    // Colour the ground green.
    glColor3f(0.137255, 0.556863, 0.137255);
    glBegin(GL_POLYGON);
        glVertex3f(200, 200, -40);
        glVertex3f(200, -200, -40);
        glVertex3f(-200, -200, -40);
        glVertex3f(-200, 200, -40);
    glEnd();
}

void drawSkybox(){
    // Draw the walls of the skybox.
    glBegin(GL_QUAD_STRIP);
        for(int i = 0; i <= 360; i++){
            // Navy blue.
            glColor3f( 0.137255, 0.137255, 0.556863);
            glVertex3f(SKYBOX_RADIUS*cos(i), SKYBOX_RADIUS*sin(i), 100);
            // Light blue.
            glColor3f(0.74902, 0.847059, 0.847059);
            glVertex3f(SKYBOX_RADIUS*cos(i), SKYBOX_RADIUS*sin(i), -100);
        }
    glEnd();
    // Draw the top of the skybox.
    glBegin(GL_TRIANGLE_FAN);
        // Navy blue.
        glColor3f( 0.137255, 0.137255, 0.556863);
        glVertex3f(0, 0, 100);
        for(int i = 0; i <= 360; i++){
            glVertex3f(SKYBOX_RADIUS*cos(i), SKYBOX_RADIUS*sin(i), 100);
        }
    glEnd();
}

// Draw the connectors that attach the rails to each other.
void drawConnectors(){
    for (int j = 0; j < (numberPoints*NUMBER_SEGMENTS)-300; j += NUMBER_SEGMENTS){
        glBegin(GL_QUAD_STRIP);
            for(int i = 0; i <= 360; i++){
                glColor3f(1.0, 0.0, 1.0);
                glVertex3f(leftRail[j][0] + cos(i)*0.05, leftRail[j][1] + cos(i)*0.05, (leftRail[j][2] + 0.05*sin(i)));
                glColor3f(1.0, 1.0, 0.0);
                glVertex3f(centerRail[j][0] + cos(i)*0.05, centerRail[j][1] + cos(i)*0.05, (centerRail[j][2] + 0.05*sin(i)));
                glColor3f(1.0, 0.0, 0.0);
                glVertex3f(columnTopRight[j][0] + cos(i)*0.05, columnTopRight[j][1] + cos(i)*0.05, (columnTopRight[j][2] + 0.05*sin(i)));
            }
        glEnd();
        glBegin(GL_QUAD_STRIP);
            for(int i = 0; i <= 360; i++){
                glColor3f(1.0, 0.0, 1.0);
                glVertex3f(rightRail[j][0] + cos(i)*0.05, rightRail[j][1] + cos(i)*0.05, (rightRail[j][2] + 0.05*sin(i)));
                glColor3f(1.0, 1.0, 0.0);
                glVertex3f(centerRail[j][0] + cos(i)*0.05, centerRail[j][1] + cos(i)*0.05, (centerRail[j][2] + 0.05*sin(i)));
                glColor3f(1.0, 0.0, 0.0);
                glVertex3f(columnTopLeft[j][0] + cos(i)*0.05, columnTopLeft[j][1] + cos(i)*0.05, (columnTopLeft[j][2] + 0.05*sin(i)));
        }
        glEnd();
    }
}

// Draw the columns to the ground.
void drawColumns(){
    for (int j = 0; j < (numberPoints*NUMBER_SEGMENTS)-300; j += NUMBER_SEGMENTS){
        glBegin(GL_QUAD_STRIP);
            for(int i = 0; i <= 360; i++){
                // Hunter green.
                glColor3f(0.13, 0.37, 0.31);
                glVertex3f(columnTopRight[j][0] + 0.1*cos(i), columnTopRight[j][1] + 0.1*sin(i), columnTopRight[j][2]);
                // Silver
                glColor3f(0.90, 0.91, 0.98);
                glVertex3f(columnTopRight[j][0] + 0.1*cos(i), columnTopRight[j][1] + 0.1*sin(i), -40);
            }
        glEnd();
        glBegin(GL_QUAD_STRIP);
            for(int i = 0; i <= 360; i++){
                // Orange.
                glColor3f(1.0, 0.5, 0.0);
                glVertex3f(columnTopLeft[j][0] + 0.1*cos(i), columnTopLeft[j][1] + 0.1*sin(i), columnTopLeft[j][2]);
                // Silver.
                glColor3f(0.90, 0.91, 0.98);
                glVertex3f(columnTopLeft[j][0] + 0.1*cos(i), columnTopLeft[j][1] + 0.1*sin(i), -40);
        }
        glEnd();
    }
}

/* Draw uniform B-spline curves for a specific section of the roller coaster.*/
void drawCurve()
{
    //Left rail purple
    glBegin(GL_QUAD_STRIP);
    for (int j = 0; j < (numberPoints*NUMBER_SEGMENTS)-300; j++){
        for(int i = 0; i <= 360; i++){
           glColor3f(1.0, 0.0, 1.0);
           glVertex3f(leftRail[j][0] + uValues[j][0]*cos(i)*0.2, leftRail[j][1] + uValues[j][1]*cos(i)*0.2, (leftRail[j][2] + 0.2*sin(i)));
            glColor3f(1.0, 0.5, 1.0);
           glVertex3f(leftRail[j+1][0] + uValues[j+1][0]*cos(i)*0.2, leftRail[j+1][1] + uValues[j+1][1]*cos(i)*0.2, (leftRail[j+1][2] + 0.2*sin(i)));
        }
    }
    glEnd();
    //Right Rail Red
    glBegin(GL_QUAD_STRIP);
    for (int j = 0; j < (numberPoints*NUMBER_SEGMENTS)-300; j++){
        for(int i = 0; i <= 360; i++){
            glColor3f(1.0, 0.0, 0.0);
            glVertex3f(rightRail[j][0] + uValues[j][0]*cos(i)*0.2, rightRail[j][1] + uValues[j][1]*cos(i)*0.2, (rightRail[j][2] + 0.2*sin(i)));
            glColor3f(1.0, 0.5, 0.5);
            glVertex3f(rightRail[j+1][0] + uValues[j+1][0]*cos(i)*0.2, rightRail[j+1][1] + uValues[j+1][1]*cos(i)*0.2, (rightRail[j+1][2] + 0.2*sin(i)));
        }
    }
    glEnd();
    //Center Rail Yellow
    glBegin(GL_QUAD_STRIP);
    for (int j = 0; j < (numberPoints*NUMBER_SEGMENTS)-300; j++){
        for(int i = 0; i <= 360; i++){
            glColor3f(1.0, 1.0, 0.0);
            glVertex3f(centerRail[j][0] + uValues[j][0]*cos(i)*0.2, centerRail[j][1] + uValues[j][1]*cos(i)*0.2, (centerRail[j][2] + 0.2*sin(i)));
            glColor3f(1.0, 1.0, 0.5);
            glVertex3f(centerRail[j+1][0] + uValues[j+1][0]*cos(i)*0.2, centerRail[j+1][1] + uValues[j+1][1]*cos(i)*0.2, (centerRail[j+1][2] + 0.2*sin(i)));
        }
    }
    glEnd();
}

/* The display calllback function for the program. This will be called for each frame of the simulation. */
void myDisplay(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    // Look at the origin and rotate around the roller coaster while maintaining a distance of 50.
    glLoadIdentity();
    if(cameraPosition == 1){
        gluLookAt(cos(phi)*100.0, sin(phi)*100.0, 60.0,
                  0.0, 0.0, 0.0,
                  0.0, 0.0, 1.0);
    }
    // Use the up vector as the v vector.
    else{
        gluLookAt(qValues[currentPosition][0], qValues[currentPosition][1], qValues[currentPosition][2]+2,
                  qValues[currentPosition][0] + dqValues[currentPosition][0],  qValues[currentPosition+1][1]+ dqValues[currentPosition][1],  qValues[currentPosition+1][2] + dqValues[currentPosition][2] +2,
                  vValues[currentPosition][0], vValues[currentPosition][1], vValues[currentPosition][2]);
    }

    glCallList(scene);
    
    glutSwapBuffers();
}

/* The timer call back for the game over last as long as the time wait macro
 * is specified for. It does not have anything to update other than the timer.
 */
void myTimer(int value){
    /*
     * Timer callback for screen between the levels.
     */
    phi = (phi + 0.005);
    
    // Update the velocity of the cart based on the height.
    velocity = (int)((1.0/2.0)*((sqrt(2.0*((10000.0/10.0)-(9.81*(qValues[currentPosition][2]+40)))))));//qValues[currentPosition][2])+40)));
    currentPosition = (currentPosition + velocity) % (NUMBER_SEGMENTS*numberPoints-(3*NUMBER_SEGMENTS));
    
    glutPostRedisplay();

    glutTimerFunc(66, myTimer, 0);
    
}
/* The reshape function for this program sets up a new window based on the specified
 aspect ratio. It also sets up the new viewing coordinates for rendering.*/
void myReshape(int w, int h)
{
    /*
     *	reshape callback function; the upper and lower boundaries of the
     *	window are at 100.0 and 0.0, respectively; the aspect ratio is
     *  determined by the aspect ratio of the viewport
     */
    
    xMax = 100.0*w/h;
    yMax = 100.0;
    
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //glOrtho(-100, 100, -100, 100, 10, 600.0);
    gluPerspective(70, (GLdouble) 2, 0.1, 300.0);
    
    glMatrixMode(GL_MODELVIEW);
}

/* We can use the arrow keys to move around the roller coaster scene. This will help with develop
 of the shape of the roller coaster and also with viewing the coaster upon completion. This will
 check for a press. */
void
keyPress(int key, int x, int y)
{
    /*
     *	this function is called when a special key is pressed; we are
     *	interested in the cursor keys only
     */
    
    switch (key)
    {
        // Key commands removed and rotation speed increased on timer callback.
    }
}

// Takes keyboard input and acts accordingly based on the key pressed.
void myKey(unsigned char key, int x, int y)
{
    switch (key)
    {
        // Exit the program and free any allocated memory by pressing q.
        case 'q':
            free(leftRail);
            free(rightRail);
            free(centerRail);
            free(columnTopRight);
            free(columnTopLeft);
            free(qValues);
            free(dqValues);
            free(ddqValues);
            free(uValues);
            free(vValues);
            free(nValues);
            exit(0);
        // Switch camera positions by pressing r.
        case 'r':
            cameraPosition = -1 * cameraPosition;
            break;
    }
}

/* We can use the arrow keys to move around the roller coaster scene. This will help with develop
 of the shape of the roller coaster and also with viewing the coaster upon completion. This will
 check for a release. */
void keyRelease(int key, int x, int y)
{
    /*
     *	this function is called when a special key is released; we are
     *	interested in the cursor keys only
     */
    
}
