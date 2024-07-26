# Irregulator
Code for making irregular and deformed cuboids of user-determined length, for input into scattering models (e.g. DDA), for example in Hamil et al. (2024) for modelling realistic KCl particles.

The code is very simple to use; to set the length of the cuboid, change these parameters on line 171 of "Irregulator_random.c":

    /* SET DIMENSIONS OF CUBOID HERE */

    cuboid_dim_a = 20;
    cuboid_dim_b = 15;
    cuboid_dim_c = 15;

To set the maximum depth of the depressions of the cuboid, change these parameters on line 177 of "Irregulator_random.c":

    /* SET MAXIMUM DEPTH OF DEPRESSIONS HERE */

    max_depth = 3;

To compile the code, open a terminal in the same folder/directory as "Irregulator_random.c" and "STAG_irregulator.c", and type:

    gcc -o Irregulator Irregulator_random.c

To run the code:

    ./Irregulator

A random cuboid will be created and a 3D rotatable image similar to the below should load (click and drag to rotate and view from other angles).

<img width="289" alt="image" src="https://github.com/user-attachments/assets/4b862d6a-8ca4-4be7-9b9f-86d774458680">

Increasing the dimensions above will increase the number of dipoles (the 'resolution' of the shape):

<img width="877" alt="image" src="https://github.com/user-attachments/assets/fada8937-560c-44a5-ada7-f8ab60fe0373">

A shape file containing the list of x,y,z dipole positions is also created: "shape.dat". This can be directly input into DDA scattering models to analyse the absorption and scattering properties using codes such as DDSCAT (Draine and Flatau, 1994), Adda (Yurkin, 2011) or CORAL (Lodge, 2024).

General algorithm:

    1) Create an irregular cuboid, from user-defined dimensions and user-defined maximum depth (of irregularities)
    2) Export this to view in STAG (image at this step disabled by default within STAG, but can be enabled by user)
    3) Adjust cube using random 'depressions' model to each surface of the cuboid (1 depression per surface, at a random position)
    4) Define 2 triangular planes at random points on opposite vertices
    5) Remove dipoles above or below these planes (depending on which vertex)
    6) Display irregular cuboid on screen using STAG.py (Simulated Three-dimensional Aerosol Geometries)

You are free to use and modify the code, but please attitribute credit to the author (Matt Lodge) if it is used.

Although every care has been taken to benchmark the code and test it under a wide range of conditions, this code is provided without a warantee of any kind. The author cannot be held liable or accept any resposibility for any claim, loss or damages that arise in connection with the use of this code.
