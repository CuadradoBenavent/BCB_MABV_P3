#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Give names of input and output files
//const char* filename = "inp.txt";
const char* filename = "tests/test-6.txt";
const char* filename_out = "trajectory-6.xyz";

//Function to allocate 2D arrays
double** malloc_2d(size_t m, size_t n) {
 
  // Allocate an array of double pointers with size m
  double** a=malloc(m*sizeof(double*));
  if (a == NULL) {
    return NULL;
  }

  // Allocate a contiguos block of memory for the 2D array elements
  a[0] = malloc (n*m*sizeof(double));
  if (a[0] == NULL) {
    free(a);
    return NULL;
  }

  // Set the pointers in the array of double prointers
  // to point to the correct locations
  for (size_t i=1 ; i<m ; i++) {
    a[i] = a[i-1]+n;
  } 

  return a;
}

// Function to free 2D arrays
void free_2d(double** a) {  // void to define function that does not return a value
  free(a[0]); // free when not used anymore
  a[0] = NULL; // set freed pointer to NULL
  free(a);
}


//Function to read the number of atoms from xyz file
size_t read_Natoms(FILE* input_file) {
  size_t Natoms;  // declare variable of type size_t
  if (fscanf(input_file, "%zu", &Natoms) != 1) { // %zu for size_t values
    fprintf(stderr, "Failed to read the number of atoms\n");  // prints error if no first line found with 1 element
    exit(1);  // terminates program if error exists
  }
  return Natoms;  // reads only first line and pointer at end
}


// Function to read the molecule data (starts reading at second line)
void read_molecule(FILE* input_file, size_t Natoms, double** coord, double* mass) {
  for (size_t i = 0; i < Natoms; i++) {  // iterating over all lines (except first)
    if (fscanf(input_file, "%lf %lf %lf %lf", &coord[i][0], &coord[i][1], &coord[i][2], &mass[i]) != 4) {  // checks for errors if line has more or less that 4 columns
      fprintf(stderr, "Failed to read data for atom %zu\n", i);
      exit(1);
    }
  }
}

 
//Define a function to compute the interatomic distances
void compute_distances(size_t Natoms, double** coord, double** distance) {
  for (size_t i = 0; i < Natoms; i++) {
    for (size_t j = 0; j < Natoms; j++) {
      //Compute distances between atoms i and j
      double dx = coord[i][0] - coord[j][0];
      double dy = coord[i][1] - coord[j][1];
      double dz = coord[i][2] - coord[j][2];
      distance[i][j] = sqrt(dx * dx + dy * dy + dz * dz);
    }
  }
}


//Define a function to calculate the Lennard-Jones potential
double V(double epsilon, double sigma, size_t Natoms, double** distance) {
  double total_potential_energy = 0.0;
  for (size_t i = 0; i < Natoms; i++) {
    for (size_t j = i + 1; j < Natoms; j++) {
      double r = distance[i][j];
      if (r > 0) {
        double sigma_r = sigma/r;
	double sigma_r6 = pow(sigma_r, 6);
	double sigma_r12 = pow(sigma_r, 12);
	double VLJ = 4.0 * epsilon * (sigma_r12 - sigma_r6);
	total_potential_energy += VLJ;
      }
    }
  }

  return total_potential_energy;
}

// Define a function to calculate kinetic energy
double T(size_t Natoms, double** velocity,double* mass){
 double total_kinetic_energy = 0.0;
 for (size_t i = 0; i < Natoms; i++){
   double velocity_squared = 0.0;
   for (size_t j = 0; j < 3; j++) {
	   velocity_squared += pow(velocity[i][j], 2);  // vix^2 + viy^2 + viz^2 
   }
   double T = mass[i] * velocity_squared;
   total_kinetic_energy += T;
 }
 total_kinetic_energy *= 0.5;   
 return total_kinetic_energy;

}


// Define a function to calculate total energy
double E(double V, double T){
  double total_E = V + T;
  return total_E;

}

// Define a function to compute the acceleration
void compute_acc(size_t Natoms, double** coord, double* mass, double** distance, double** acceleration, double epsilon, double sigma){
	for (size_t i = 0; i < Natoms; i++){
	  double axi = 0.0;
	  double ayi = 0.0;
	  double azi = 0.0;
	  for (size_t j = 0; j < Natoms; j++) { // Include all other atoms
            if (i != j) { // Exclude self-interaction
              double r = distance[i][j];
	      if (r>0){
                double sigma_r = sigma/r;
      	        double sigma_r6 = pow(sigma_r, 6);
      	        double sigma_r12 = pow(sigma_r, 12);
      	        double U = 24.0 * epsilon / r * (sigma_r6 - 2.0 * sigma_r12);
                axi -= (1 / mass[i] * U * ( coord[i][0] - coord[j][0]) / r); 
                ayi -= (1 / mass[i] * U * ( coord[i][1] - coord[j][1]) / r);
                azi -= (1 / mass[i] * U * ( coord[i][2] - coord[j][2]) / r);
	      }
	    }
	  }
	  acceleration[i][0] = axi;
	  acceleration[i][1] = ayi; 
	  acceleration[i][2] = azi; 
	}
}



void write_xyz(FILE* output_file, size_t Natoms, double** coord, double potential_energy, double kinetic_energy, double total_energy, size_t iteration) {
  fprintf(output_file, "Iteration: %zu\n", iteration);
  fprintf(output_file, "%zu\n", Natoms);
  fprintf(output_file, "Potential Energy: %lf, Kinetic Energy: %lf, Total Energy: %lf\n", potential_energy, kinetic_energy, total_energy);
  for (size_t i = 0; i < Natoms; i++) {
    fprintf(output_file, "Ar %lf %lf %lf\n", coord[i][0], coord[i][1], coord[i][2]);
  }
}
	     

/////////////////////////////////////////////////////////////////////////////////////////////

int main() {

  FILE* input_file = fopen(filename, "r");  // open inp.txt for reading
  if (input_file == NULL) {
    printf("Failed to open file");
    return 1;  // exit with error if no file found
  }
  printf("Reading molecular data from: %s\n", filename);

  //Read the number of atoms
  size_t Natoms = read_Natoms(input_file); 
  //printf("Number of atoms: %zu\n", Natoms);


  //Allocate memory for coords and masses
  double** coord = malloc_2d(Natoms, 3);  // 3 coordinates per atom
  double* mass = malloc(Natoms * sizeof(double));  // 1 mass per atom
  if (coord == NULL || mass == NULL) {
    fprintf(stderr, "Memory allocation failed\n");
    fclose(input_file);
    return 1;
  }
  
  //Read molecule data
  read_molecule(input_file, Natoms, coord, mass);
  //for (size_t i = 0; i < Natoms; i++) {
  //  printf("Atom %zu: x = %.1lf, y = %.1lf, z = %.1lf, mass = %.3lf\n", i + 1, coord[i][0], coord[i][1], coord[i][2], mass[i]);
  //}

  fclose(input_file);


  //Calculate interatomic distances
  double** distance = malloc_2d(Natoms, Natoms);
  if (distance == NULL) {
    fprintf(stderr, "Memory allocation for distance array has failed\n");
    free_2d(coord);
    free(mass);
    return 1;
  }


  //Compute and print the distances
  compute_distances(Natoms, coord, distance);
  //printf("\nInternuclear distances:\n");
  //for (size_t i = 0; i < Natoms; i++) {
  //  for (size_t j = 0; j < Natoms; j++) {
  //    printf("Interatomic distance[%zu][%zu] = %lf\n", i, j, distance[i][j]);  // %lf to print double
  //  }
  //}

  //Compute and print LJ potential energy
  double epsilon = 0.0661;
  double sigma = 0.3345;
  double total_V = V(epsilon, sigma, Natoms, distance);
  //printf("\nTotal potential energy: %lf J/mol\n", total_V);

  //Compute and print kinetic energy
  double** velocity = malloc_2d(Natoms, 3);
  for (size_t i = 0; i < Natoms; i++) {
    for (size_t j = 0; j < 3; j++){
      velocity[i][j] = 0.0;
    }
  }  // initialize velocities to 0
  double total_T = T(Natoms, velocity, mass);
  //printf("Total kinetic energy: %lf J/mol\n", total_T);
  

  //Compute and print total energy
  double total_E = E(total_V, total_T);
  //printf("Total energy: %lf J/mol\n", total_E);



  //Compute acceleration
  double** acceleration = malloc_2d(Natoms, 3);
  if (acceleration == NULL) {
    fprintf(stderr, "Memory allocation for acceleration array has failed\n");
    free_2d(coord);
    free_2d(distance);
    free(mass);
    free_2d(velocity);
    return 1;
  }
  compute_acc(Natoms, coord, mass, distance, acceleration, epsilon, sigma);

   //Define output file
   FILE* output_file = fopen(filename_out, "w");
   if (output_file == NULL) {
     fprintf(stderr, "Failed to open trajectory file for writing\n");
     free_2d(coord);
     free_2d(distance);
     free(mass);
     free_2d(velocity);
     free_2d(acceleration);
     return 1;
   }

  
  //Verlet algorithm
  double dt = 0.2;
  size_t nsteps = 1000;
  //printf("\n\n*****Beginning molecular dynamics simulation based on the Verlet algorithm*****\n\n");
  size_t M = 10; //To print trajectory every M numbers of steps 
  for (size_t n = 0; n < nsteps; n++){
    for (size_t i = 0; i < Natoms; i++){
      for (size_t j = 0; j < 3; j++){
        coord[i][j] += velocity[i][j] * dt + acceleration[i][j] * pow(dt, 2) / 2.0;  // coord update
        velocity[i][j] += 0.5 * acceleration[i][j] * dt;
      }
    }
    compute_distances(Natoms, coord, distance); //a^n
    compute_acc(Natoms, coord, mass, distance, acceleration, epsilon, sigma);
    for (size_t i = 0; i < Natoms; i++){
      for (size_t j = 0; j < 3; j++){
        velocity[i][j] += 0.5 * acceleration[i][j] * dt;  //a^n+1
      }
    }
    double total_V = V(epsilon, sigma, Natoms, distance);
    double total_T = T(Natoms, velocity, mass);
    double total_E = E(total_V, total_T);


 //Write trajectory every M steps
    if (n % M == 0) {
      write_xyz(output_file, Natoms, coord, total_V, total_T, total_E, n/M);
    }
   
    //printf("\nIteration: %zu", n+1);  // to print starting at 1
    //printf("\nEnergies: \n   Potential %lf J/mol", total_V);
    //printf("\n   Kinetic %lf J/mol", total_T);
    //printf("\n   Total %lf J/mol", total_E);
    //printf("\nCoordinates:\n");
    //for (size_t i = 0; i < Natoms; i++) {
    //  for (size_t j = 0; j < 3; j++) {
    //    printf("%lf ", coord[i][j]);  // %lf to print double
    //  } 
    //printf("\n"); 
    //}
  }

  fclose(output_file);

  printf("An output file has been created: %s \n", filename_out);

  // Free allocated memory
  free_2d(coord);
  free_2d(distance);
  free(mass);
  free_2d(velocity);
  free_2d(acceleration);

  return 0;
}
