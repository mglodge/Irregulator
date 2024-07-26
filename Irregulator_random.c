/*

28/03/24 - Code to create irregular cuboids and display them on the screen. Algorithm is as follows:

    1) Create an irregular cuboid, from given dimensions
    2) Export this to view in STAG
    3) Adjust cube using random 'depressions' model to each surface of the cuboid (1 depression per surface, at a random position)
    4) Define 2 triangular planes at random points on opposite vertices
    5) Remove dipoles above or below these planes (depending on which vertex)
    6) Display irregular cuboid on screen using STAG.py (Simulated Three-dimensional Aerosol Geometries)
    
*/
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <time.h>


double find_dipole_height(int v, int w, int dep_v, int dep_w, int v_max, int w_max, int depth){
    double get_dipole_height, test_L, test_R, fraction_to_edge;

        // function to determine dipole height for each (v,w) position on each surface

        if((v==dep_v)&&(w==dep_w)){
            if(depth==1){
                get_dipole_height = -1.0*depth; //if the current dipole is at the depression point, just set it equal to the max depth
            }
            else{
                get_dipole_height = -1.0*(depth-1); //slight aesthetic adjustment: raise the depression point up by +1, so that it sits with the same height as the dipoles immediatedly adjacent to it. The condition above prevents this occuring on the only case we wouldn't want it to -- when +1 would mean being back on the surface (height=0)
            }
        }
        else if( (v==0) || (v==v_max) || (w==0) || (w==w_max) ){
            get_dipole_height = 0; //if the current dipole is at the edge of the grid, set the depression height == 0
        }
        else{

            //determine how close this dipole is to the edge of the rectangular surface, as a fraction of the distance between the depression and the edge of the rectangular surface (at the same gradient in (v,w) space)
            
            if(v>=dep_v){ //if we are to the right of the depression
                if(w>=dep_w){ //if we are above the depression
                    
                    //UPPER-RIGHT
                    test_L= (float)(v-dep_v)/(float)(v_max-dep_v);
                    test_R= (float)(w-dep_w)/(float)(w_max-dep_w);

                    if(test_L>=test_R){
                        //we hit wall v first
                        fraction_to_edge = test_L;
                    }
                    else{
                        //we hit wall w first
                        fraction_to_edge = test_R;
                    }
                }
                else{ //if we are below the depression
                    
                    //BOTTOM-RIGHT
                    test_L= (float)(v-dep_v)/(float)(v_max-dep_v);
                    test_R= (float)(dep_w-w)/(float)(dep_w);

                    if(test_L>=test_R){
                        //we hit wall v first
                        fraction_to_edge = test_L;
                    }
                    else{
                        //we hit wall w first
                        fraction_to_edge = test_R;
                    }
                }
            }
            else{ //if we are to the left of the depression
                if(w>=dep_w){ //if we are above the depression
                    
                    //UPPER-LEFT
                    test_L= (float)(dep_v-v)/(float)(dep_v);
                    test_R= (float)(w-dep_w)/(float)(w_max-dep_w);

                    if(test_L>=test_R){
                        //we hit wall v first
                        fraction_to_edge = test_L;
                    }
                    else{
                        //we hit wall w first
                        fraction_to_edge = test_R;
                    }
                }
                else{ //if we are below the depression
                    
                    //BOTTOM-LEFT
                    test_L= (float)(dep_v-v)/(float)(dep_v);
                    test_R= (float)(dep_w-w)/(float)(dep_w);

                    if(test_L>=test_R){
                        //we hit wall v first
                        fraction_to_edge = test_L;
                    }
                    else{
                        //we hit wall w first
                        fraction_to_edge = test_R;
                    }
                }
            }

            //work out the height of the dipole at this (v,w) position for the adjusted surface
            get_dipole_height = -1*floor((1.0-fraction_to_edge) * depth);
        }
    return get_dipole_height;
}
void get_triangular_plane(int vertex_A[3], int vertex_B[3], int vertex_C[3], double *normal, double *k_coeff) { //the notation double* means 'an array', and the notation *normal means we alloate a pointer to return multiple variables (normal[0,1,2] and k_coeff) from this function
    double B_minus_A[3], C_minus_A[3];

    // function to define the coordinates of the triangular plane that joins vertices A, B and C

    /* find two vectors that are in this plane */

    B_minus_A[0]= vertex_B[0] - vertex_A[0]; 
    B_minus_A[1]= vertex_B[1] - vertex_A[1]; 
    B_minus_A[2]= vertex_B[2] - vertex_A[2]; 

    C_minus_A[0]= vertex_C[0] - vertex_A[0]; 
    C_minus_A[1]= vertex_C[1] - vertex_A[1]; 
    C_minus_A[2]= vertex_C[2] - vertex_A[2]; 

    /* find the normal, by calculating the cross-product of 'B-A' X 'C-A' */
    
    normal[0]= B_minus_A[1]*C_minus_A[2] - B_minus_A[2]*C_minus_A[1];
    normal[1]= B_minus_A[2]*C_minus_A[0] - B_minus_A[0]*C_minus_A[2];
    normal[2]= B_minus_A[0]*C_minus_A[1] - B_minus_A[1]*C_minus_A[0];

    //calculate coefficient k by substituting (x,y,z) values of vector A into the equation:   (n_x * x)  +  (n_y * y)  +  (n_z * z)  =  -k
    *k_coeff = -(normal[0]*vertex_A[0]) - (normal[1]*vertex_A[1]) - (normal[2]*vertex_A[2]);
}
FILE* DDSCAT_outfile;
FILE* original_grid_outfile;
FILE* new_grid_outfile;
int diagnostics, cuboid_dim_a, cuboid_dim_b, cuboid_dim_c, max_lattice_dim, v, w, v_max, w_max, into_page_max, i, j, k, x, y , z, original_N, original_lattice_dim, new_lattice_dim, dipole_count, JA, IX, IY, IZ, ICOMPX, ICOMPY, ICOMPZ, STAG_lattice_dim, edgecase, dep_v, dep_w, depth, dipole_height, max_depth, vertex_A[3], vertex_B[3], vertex_C[3];
double STAG_odd_even_offset, min[3], max[3], STAG_offset[3];
double normal[3], k_coeff, triangular_plane;
char buf[1000];
int*** original_grid;
int*** new_grid;
double** dipole_info;
double** STAG_dipole_positions;



int main()
{


    printf("\n\n ---------------------------------------------------------------------------------------------------------------------");
    printf("\n ---------------------------------------------------------------------------------------------------------------------");
    printf("\n\n                                                  WELCOME TO IRREGULATOR!                            \n\n");


    printf("\n                  ________________                  ________               ");
    printf("\n                 |       |        |                |        |              ");
    printf("\n                 |       |        |                |        |              ");
    printf("\n                 |_______|________|                |________|________      ");
    printf("\n                 |       |        |      --->      |        |        |     ");
    printf("\n                 |       |        |                |        |        |     ");
    printf("\n                 |_______|________|                |________|________|     ");


    printf("\n\n\n ---------------------------------------------------------------------------------------------------------------------");
    printf("\n ---------------------------------------------------------------------------------------------------------------------\n\n");


    /* SET DIMENSIONS OF CUBOID HERE */

    cuboid_dim_a = 20;
    cuboid_dim_b = 15;
    cuboid_dim_c = 15;
    
    /* SET MAXIMUM DEPTH OF DEPRESSIONS HERE */

    max_depth = 3;





    /* begin code */

    original_N = cuboid_dim_a * cuboid_dim_b * cuboid_dim_c;

    diagnostics=0; //set = 1 to print diagnostic statements

    
    /* initialise matrix to store dipoles, now that we know how many there are in the cuboid */
    dipole_info=(double**)malloc(original_N*sizeof(double*));
    for (i=0; i<original_N; i++) {
        dipole_info[i]=(double*)malloc((3)*sizeof(double));   //matrix will record each of the values for "IX IY IZ ICOMPX ICOMPY ICOMPZ" in array elements [0]->[5]
    }

    /* create grid of x,y,z positions for cuboid */

    k=0;
    for(x=0; x<cuboid_dim_a; x++){
        for(y=0; y<cuboid_dim_b; y++){
            for(z=0; z<cuboid_dim_c; z++){
                dipole_info[k][0] = x;
                dipole_info[k][1] = y;
                dipole_info[k][2] = z;
                k++;
            }
        }
    }
    
    printf("\n\n %d dipole cuboid (%d x %d x %d) successfully created.", k, cuboid_dim_a, cuboid_dim_b, cuboid_dim_c);
    
        /* save original resolution output to STAG_spherify data file (to visualise it in 3D) */

    printf(" Exporting original data to S.T.A.G...");

    // find the largest dimension

    if(cuboid_dim_a>=cuboid_dim_b){
        if(cuboid_dim_a>=cuboid_dim_c){
            max_lattice_dim = cuboid_dim_a; // a is the largest dimension
        }
        else{
            max_lattice_dim = cuboid_dim_c; // c is the largest dimension
        }
    }
    else{
        if(cuboid_dim_b>=cuboid_dim_c){
            max_lattice_dim = cuboid_dim_b; // b is the largest dimension
        }
        else{
            max_lattice_dim = cuboid_dim_c; // c is the largest dimension
        }
    }

    original_grid_outfile=fopen("original_cuboid.txt","w"); //open file for saving dipole positions

    for(x=0; x<cuboid_dim_a; x++){
        for(y=0; y<cuboid_dim_b; y++){
            for(z=0; z<cuboid_dim_c; z++){
                fprintf(original_grid_outfile,"%d, %d, %d\n", x,y,z); //save x-y-z coords of all original cuboid dipoles
            }
        }
    }

    /* The final row is extra data needed for the python visualisation S.T.A.G program, NOT a dipole! */
    fprintf(original_grid_outfile,"%d, %d, %d\n", max_lattice_dim, original_N, max_lattice_dim); //the final row contains the number of dipoles, the grid size, and a random number just to keep the shape of three columns for python to read */
    fclose(original_grid_outfile);

    printf(" Export complete.\n");

    /* create a binary grid (1s or 0s) to store positions of dipoles that we want to keep */

    original_grid = (int***)malloc(cuboid_dim_a*sizeof(int**));
    for (i=0;i<cuboid_dim_a;i++) {
        original_grid[i] = (int**)malloc(cuboid_dim_b*sizeof(int*));
        for (j=0;j<cuboid_dim_b;j++) {
            original_grid[i][j] = (int*)malloc(cuboid_dim_c*sizeof(int));
        }
    }

    // initialise grid
    for(x=0;x<cuboid_dim_a;x++){
        for(y=0;y<cuboid_dim_b;y++){
            for(z=0;z<cuboid_dim_c;z++){
                original_grid[x][y][z]=1; //set values at all positions to 1 initially (we keep all dipoles in the cuboid, unless they are deleted by the algorithm below)
            }
        }
    }


        /* CREATE RANDOM DEPRESSIONS ON SURFACES */

    srand(time(0)); // initialise random number generator with seed based on computer time

    /* front surface (y=0, x-z plane) */

    //choose coordinates to place the depression point -- uses (v,w) space, matching the orientation for any given surface plane
    dep_v = (rand() % (cuboid_dim_a-2)) + 1; // creates list of random numbers between 1 -> max-2 (so that depression points are not on boundaries)
    dep_w = (rand() % (cuboid_dim_c-2)) + 1; // creates list of random numbers between 1 -> max-2 (so that depression points are not on boundaries)
    depth = rand() % (max_depth+1); // creates list of random numbers between "0 -> max_depth" important keep as a positive number at this stage due to ceil() function -- we will make it negative later

    printf("\n  For the front surface, the depression is at (x,z) = (%d, %d) with a depth of %d (remember that we will actually only go as deep as 'depth-1')", dep_v, dep_w, depth);

    v_max = cuboid_dim_a-1; // x - set the boundaries of the surface rectangle from whichever plane we are looking at each time
    w_max = cuboid_dim_c-1; // z - remember to -1 from all three coords because we begin at element [0]
    into_page_max = cuboid_dim_b-1; // y -- find the maximum depth of the shape into the page, used to remove dipoles from 'max' surfaces rather than 'origin' surfaces

    //go through each dipole in this rectangular (v,w) surface one-by-one
    for(w=0; w<=w_max; w++){
        for(v=0; v<=v_max; v++){

            dipole_height = find_dipole_height(v, w, dep_v, dep_w, v_max, w_max, depth);

            //delete cubes above this point by setting their values == 0 in the binary grid
            for(i=0; i<(-1*dipole_height); i++){
                original_grid[v][i][w] = 0;   // delete cubes between y=0 and y=dipole_height (simply given by the integers 'i') at this dipole (v,w)
            }

        }
    }


    /* back surface (y=max, x-z plane) */   

    //choose coordinates to place the depression point -- uses (v,w) space, matching the orientation for any given surface plane
    dep_v = (rand() % (cuboid_dim_a-2)) + 1; // creates list of random numbers between 1 -> max-2 (so that depression points are not on boundaries)
    dep_w = (rand() % (cuboid_dim_c-2)) + 1; // creates list of random numbers between 1 -> max-2 (so that depression points are not on boundaries)
    depth = rand() % (max_depth+1); // creates list of random numbers between "0 -> max_depth" important keep as a positive number at this stage due to ceil() function -- we will make it negative later

    printf("\n  For the back surface, the depression is at (x,z) = (%d, %d) with a depth of %d (remember that we will actually only go as deep as 'depth-1')", dep_v, dep_w, depth);

    v_max = cuboid_dim_a-1; // x - set the boundaries of the surface rectangle from whichever plane we are looking at each time
    w_max = cuboid_dim_c-1; // z - remember to -1 because we begin at element [0]
    into_page_max = cuboid_dim_b-1; // y - find the maximum depth of the shape into the page, used to remove dipoles from 'max' surfaces rather than 'origin' surfaces

    //go through each dipole in this rectangular (v,w) surface one-by-one
    for(w=0; w<=w_max; w++){
        for(v=0; v<=v_max; v++){

            dipole_height = find_dipole_height(v, w, dep_v, dep_w, v_max, w_max, depth);

            //delete cubes above this point by setting their values == 0 in the binary grid
            for(i=0; i<(-1*dipole_height); i++){
                original_grid[v][into_page_max-i][w] = 0;   // delete cubes between y=y_max and y=dipole_height (simply given by the integers 'into_page_max-i') at this dipole (v,w)
            }
            
        }
    }


    /* left surface (x=0, y-z plane) */

    //choose coordinates to place the depression point -- uses (v,w) space, matching the orientation for any given surface plane
    dep_v = (rand() % (cuboid_dim_b-2)) + 1; // creates list of random numbers between 1 -> max-2 (so that depression points are not on boundaries)
    dep_w = (rand() % (cuboid_dim_c-2)) + 1; // creates list of random numbers between 1 -> max-2 (so that depression points are not on boundaries)
    depth = rand() % (max_depth+1); // creates list of random numbers between "0 -> max_depth" important keep as a positive number at this stage due to ceil() function -- we will make it negative later

    printf("\n  For the left surface, the depression is at (y,z) = (%d, %d) with a depth of %d (remember that we will actually only go as deep as 'depth-1')", dep_v, dep_w, depth);

    v_max = cuboid_dim_b-1; // y - set the boundaries of the surface rectangle from whichever plane we are looking at each time
    w_max = cuboid_dim_c-1; // z - remember to -1 because we begin at element [0]
    into_page_max = cuboid_dim_a-1; // x - find the maximum depth of the shape into the page, used to remove dipoles from 'max' surfaces rather than 'origin' surfaces

    //go through each dipole in this rectangular (v,w) surface one-by-one
    for(w=0; w<=w_max; w++){
        for(v=0; v<=v_max; v++){

            dipole_height = find_dipole_height(v, w, dep_v, dep_w, v_max, w_max, depth);

            //delete cubes above this point by setting their values == 0 in the binary grid
            for(i=0; i<(-1*dipole_height); i++){
                original_grid[i][v][w] = 0;   // delete cubes between x=0 and x=dipole_height (simply given by the integers 'i') at this dipole (v,w)
            }

        }
    }


    /* right surface (x=max, y-z plane) */

    //choose coordinates to place the depression point -- uses (v,w) space, matching the orientation for any given surface plane
    dep_v = (rand() % (cuboid_dim_b-2)) + 1; // creates list of random numbers between 1 -> max-2 (so that depression points are not on boundaries)
    dep_w = (rand() % (cuboid_dim_c-2)) + 1; // creates list of random numbers between 1 -> max-2 (so that depression points are not on boundaries)
    depth = rand() % (max_depth+1); // creates list of random numbers between "0 -> max_depth" important keep as a positive number at this stage due to ceil() function -- we will make it negative later

    printf("\n  For the right surface, the depression is at (y,z) = (%d, %d) with a depth of %d (remember that we will actually only go as deep as 'depth-1')", dep_v, dep_w, depth);

    v_max = cuboid_dim_b-1; // y - set the boundaries of the surface rectangle from whichever plane we are looking at each time
    w_max = cuboid_dim_c-1; // z - remember to -1 because we begin at element [0]
    into_page_max = cuboid_dim_a-1; // x -find the maximum depth of the shape into the page, used to remove dipoles from 'max' surfaces rather than 'origin' surfaces

    //go through each dipole in this rectangular (v,w) surface one-by-one
    for(w=0; w<=w_max; w++){
        for(v=0; v<=v_max; v++){

            dipole_height = find_dipole_height(v, w, dep_v, dep_w, v_max, w_max, depth);

            //delete cubes above this point by setting their values == 0 in the binary grid
            for(i=0; i<(-1*dipole_height); i++){
                original_grid[into_page_max-i][v][w] = 0;   // delete cubes between x=x_max and x=dipole_height (simply given by the integers 'into_page_max-i') at this dipole (v,w)
            }

        }
    }
    

    /* bottom surface (z=0, x-y plane) */

    //choose coordinates to place the depression point -- uses (v,w) space, matching the orientation for any given surface plane
    dep_v = (rand() % (cuboid_dim_a-2)) + 1; // creates list of random numbers between 1 -> max-2 (so that depression points are not on boundaries)
    dep_w = (rand() % (cuboid_dim_b-2)) + 1; // creates list of random numbers between 1 -> max-2 (so that depression points are not on boundaries)
    depth = rand() % (max_depth+1); // creates list of random numbers between "0 -> max_depth" important keep as a positive number at this stage due to ceil() function -- we will make it negative later

    printf("\n  For the bottom surface, the depression is at (x,y) = (%d, %d) with a depth of %d (remember that we will actually only go as deep as 'depth-1')", dep_v, dep_w, depth);

    v_max = cuboid_dim_a-1; // x - set the boundaries of the surface rectangle from whichever plane we are looking at each time
    w_max = cuboid_dim_b-1; // y - remember to -1 from all three coords because we begin at element [0]
    into_page_max = cuboid_dim_c-1; // z -- find the maximum depth of the shape into the page, used to remove dipoles from 'max' surfaces rather than 'origin' surfaces

    //go through each dipole in this rectangular (v,w) surface one-by-one
    for(w=0; w<=w_max; w++){
        for(v=0; v<=v_max; v++){

            dipole_height = find_dipole_height(v, w, dep_v, dep_w, v_max, w_max, depth);

            //delete cubes above this point by setting their values == 0 in the binary grid
            for(i=0; i<(-1*dipole_height); i++){
                original_grid[v][w][i] = 0;   // delete cubes between z=0 and z=dipole_height (simply given by the integers 'i') at this dipole (v,w)
            }

        }
    }


    /* top surface (z=max, x-y plane) */

    //choose coordinates to place the depression point -- uses (v,w) space, matching the orientation for any given surface plane
    dep_v = (rand() % (cuboid_dim_a-2)) + 1; // creates list of random numbers between 1 -> max-2 (so that depression points are not on boundaries)
    dep_w = (rand() % (cuboid_dim_b-2)) + 1; // creates list of random numbers between 1 -> max-2 (so that depression points are not on boundaries)
    depth = rand() % (max_depth+1); // creates list of random numbers between "0 -> max_depth" important keep as a positive number at this stage due to ceil() function -- we will make it negative later

    printf("\n  For the top surface, the depression is at (y,z) = (%d, %d) with a depth of %d (remember that we will actually only go as deep as 'depth-1')", dep_v, dep_w, depth);

    v_max = cuboid_dim_a-1; // x - set the boundaries of the surface rectangle from whichever plane we are looking at each time
    w_max = cuboid_dim_b-1; // y - remember to -1 from all three coords because we begin at element [0]
    into_page_max = cuboid_dim_c-1; // z -- find the maximum depth of the shape into the page, used to remove dipoles from 'max' surfaces rather than 'origin' surfaces

    //go through each dipole in this rectangular (v,w) surface one-by-one
    for(w=0; w<=w_max; w++){
        for(v=0; v<=v_max; v++){

            dipole_height = find_dipole_height(v, w, dep_v, dep_w, v_max, w_max, depth);

            //delete cubes above this point by setting their values == 0 in the binary grid
            for(i=0; i<(-1*dipole_height); i++){
                original_grid[v][w][into_page_max-i] = 0;   // delete cubes between z=0 and z=dipole_height (simply given by the integers 'into_page_max-i') at this dipole (v,w)
            }

        }
    }

            /* REMOVE TRIANGULAR SECTIONS */

    /* REMOVE TRIANGLE NEAR ORIGIN (0,0,0) */

    /* choose random vertices for the triangles */

    vertex_A[0] = (rand() % (cuboid_dim_a-6))+3; // selects a random number between 3 -> 'dim_a - 3'
    vertex_A[1] = 0;
    vertex_A[2] = 0;
    
    vertex_B[0] = 0;
    vertex_B[1] = (rand() % (cuboid_dim_a-6))+3; // selects a random number between 3 -> 'dim_a - 3'
    vertex_B[2] = 0;

    vertex_C[0] = 0;
    vertex_C[1] = 0;
    vertex_C[2] = (rand() % (cuboid_dim_a-6))+3; // selects a random number between 3 -> 'dim_a - 3'

    printf("\n\n ORIGIN TRIANGLE:");

    printf("\n\n\t Vertex A: \n\t\t\t %d \n\t\t\t %d \n\t\t\t %d", vertex_A[0], vertex_A[1], vertex_A[2]);
    printf("\n\n\t Vertex B: \n\t\t\t %d \n\t\t\t %d \n\t\t\t %d", vertex_B[0], vertex_B[1], vertex_B[2]);
    printf("\n\n\t Vertex C: \n\t\t\t %d \n\t\t\t %d \n\t\t\t %d", vertex_C[0], vertex_C[1], vertex_C[2]);

    /* find the triangular plane joining these three vertices */
    get_triangular_plane(vertex_A, vertex_B, vertex_C, &normal, &k_coeff);

    printf("\n\n\t Normal: \n\t\t\t %f \n\t\t\t %f \n\t\t\t %f \n\t\t\t %f (k coefficient)", normal[0], normal[1], normal[2], k_coeff);

    /* go through every dipole, and if they are BELOW this plane, remove them:   (n_x * x)  +  (n_y * y)  +  (n_z * z)  +  k = 0    */

    for(x=0; x<cuboid_dim_a; x++){
        for(y=0; y<cuboid_dim_b; y++){
            for(z=0; z<cuboid_dim_c; z++){
                
                triangular_plane = normal[0]*x + normal[1]*y + normal[2]*z + k_coeff;

                if(triangular_plane<0){
                    original_grid[x][y][z] = 0;   // delete cubes that are BELOW the triangular plane
                }

            }
        }
    }

     /* REMOVE TRIANGLE ON OPPOSITE SIDE FROM NEAR MAXIMUMS (x_max, y_max, z_max) */

    /* choose random vertices for the triangles */

    vertex_A[0] = (rand() % (cuboid_dim_a-6))+3; // selects a random number between 3 -> 'dim_a - 3'
    vertex_A[1] = cuboid_dim_b-1;
    vertex_A[2] = cuboid_dim_c-1; //set equal to the maximum values
    
    vertex_B[0] = cuboid_dim_a-1;
    vertex_B[1] = (rand() % (cuboid_dim_a-6))+3; // selects a random number between 3 -> 'dim_a - 3'
    vertex_B[2] = cuboid_dim_c-1;

    vertex_C[0] = cuboid_dim_a-1;
    vertex_C[1] = cuboid_dim_b-1;
    vertex_C[2] = (rand() % (cuboid_dim_a-6))+3; // selects a random number between 3 -> 'dim_a - 3'

    printf("\n\n MAXIMA TRIANGLE:");
    
    printf("\n\n\t Vertex A: \n\t\t\t %d \n\t\t\t %d \n\t\t\t %d", vertex_A[0], vertex_A[1], vertex_A[2]);
    printf("\n\n\t Vertex B: \n\t\t\t %d \n\t\t\t %d \n\t\t\t %d", vertex_B[0], vertex_B[1], vertex_B[2]);
    printf("\n\n\t Vertex C: \n\t\t\t %d \n\t\t\t %d \n\t\t\t %d", vertex_C[0], vertex_C[1], vertex_C[2]);

    /* find the triangular plane joining these three vertices */
    get_triangular_plane(vertex_A, vertex_B, vertex_C, &normal, &k_coeff);

    printf("\n\n\t Normal: \n\t\t\t %f \n\t\t\t %f \n\t\t\t %f \n\t\t\t %f (k coefficient)", normal[0], normal[1], normal[2], k_coeff);

    /* go through every dipole, and if they are ABOVE this plane, remove them:   (n_x * x)  +  (n_y * y)  +  (n_z * z)  +  k = 0    */

    for(x=0; x<cuboid_dim_a; x++){
        for(y=0; y<cuboid_dim_b; y++){
            for(z=0; z<cuboid_dim_c; z++){
                
                triangular_plane = normal[0]*x + normal[1]*y + normal[2]*z + k_coeff;

                if(triangular_plane>0){
                    original_grid[x][y][z] = 0;   // delete cubes that are ABOVE the triangular plane
                }

            }
        }
    }




    /* go through resulting grid and count the number of dipoles still remaining */
    dipole_count=0;
    for(x=0; x<cuboid_dim_a; x++){
        for(y=0; y<cuboid_dim_b; y++){
            for(z=0; z<cuboid_dim_c; z++){
                if(original_grid[x][y][z]==1){
                    dipole_count++;
                } 
            }
        }
    }

    printf("\n\n\n OVERALL RESULT: %d dipoles in irregular cube (out of an original %d)...", dipole_count, original_N);

    /* save irregularised cuboid to STAG_spherify data file (to visualise it in 3D) */

    printf("\n\n Exporting irregularised cuboid to S.T.A.G...");

    new_grid_outfile=fopen("irregular_cuboid.txt","w"); //open file for saving dipole positions

    for(x=0; x<cuboid_dim_a; x++){
        for(y=0; y<cuboid_dim_b; y++){
            for(z=0; z<cuboid_dim_c; z++){
                if(original_grid[x][y][z]==1){
                    fprintf(new_grid_outfile,"%d, %d, %d\n", x,y,z); //save x-y-z coords of all cuboid dipoles
                }
            }
        }
    }

    /* The final row is extra data needed for the python visualisation S.T.A.G program, NOT a dipole! */
    fprintf(new_grid_outfile,"%d, %d, %d\n", max_lattice_dim, dipole_count, max_lattice_dim); //the final row contains the number of dipoles, the grid size, and a random number just to keep the shape of three columns for python to read */
    fclose(new_grid_outfile);

    printf(" Export complete.\n");



    /* SAVE DATA IN DDSCAT FORMAT */

    /* save output to DDSCAT type file */
    printf(" Exporting irregularised cuboid to shape.dat in DDSCAT format...");
    DDSCAT_outfile=fopen("shape.dat","w"); //open file for saving dipole positions

    fprintf(DDSCAT_outfile," Irregular %d x %d x %d Cuboid, created with irregulator.c \n", cuboid_dim_a, cuboid_dim_b, cuboid_dim_c);
    fprintf(DDSCAT_outfile," %6d   = NAT \n", dipole_count);
    fprintf(DDSCAT_outfile,"   1.0000    .0000    .0000 = A_1 vector \n");
    fprintf(DDSCAT_outfile,"    .0000   1.0000    .0000 = A_2 vector \n");
    fprintf(DDSCAT_outfile,"  1.00000   1.0000   1.0000 = lattice spacings (d_x,d_y,d_z)/d \n");
    fprintf(DDSCAT_outfile," 0 0 0 = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d for dipole 0 0 0 \n");
    fprintf(DDSCAT_outfile,"     JA  IX  IY  IZ ICOMP(x,y,z) \n");

    i=0;
    for(x=0; x<cuboid_dim_a; x++){
        for(y=0; y<cuboid_dim_b; y++){
            for(z=0; z<cuboid_dim_c; z++){
                if(original_grid[x][y][z]==1){
                    fprintf(DDSCAT_outfile,"     %7d  %4d  %4d  %4d     1     1     1 \n", i+1, x, y, z);; //save x-y-z coords of all cuboid dipoles
                    i++;
                }
            }
        }
    }

    fclose(DDSCAT_outfile);
    
    printf(" Export complete.\n");

    // run xSTAG_spherify on each file - views both the original and higher resolution images in 3D

    //system("xSTAG_irregulator.bat"); //WINDOWS VERSION: opens a batch file with a command to run STAG_spherify as a python script
    system("python STAG_irregulator.py"); // MAC version -- open file in python */

    /* free memory for arrays */

    for(i=0;i<original_lattice_dim;i++){
        for(j=0;j<original_lattice_dim;j++){
            free((void*)original_grid[i][j]);
        }
        free((void*)original_grid[i]);
    }
    free((void*)original_grid);


    for (i=0; i<original_N; i++) {
        free((void*)dipole_info[i]);
    }
    free((void*)dipole_info);

    return 0;
}
