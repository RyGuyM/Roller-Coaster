Roller-Coaster
==============

Roller Coaster written in C using OpenGL. This program is used to demonstrate C2 continuous motion of a camera along a roller coaster track. To perform this task uniform b-splines were used to created a track for the camera to ride along. To run this program we need to run the following in the command line:


    $ gcc -std=c99 -o Asteroids Asteroids.c -framework OPENGL -framework GLUT 

    $ ./Asteroids

To change the view from being in the world perspective where the track is rotating to riding the coaster, use the following commands in the program:

    Press 'r' - Toggle the position of the camera.
    
    Press 'q' - Quit the program.
