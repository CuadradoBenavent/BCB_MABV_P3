#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double** malloc_2d(size_t m, size_t n) {  // ** pointer to a pointer
  
  // Allocate an array of double pointers with size m
  double** a = malloc(m*sizeof(double*)); // request m double pointers
  if (a == NULL) {  // to check that the allocation ok
    return NULL;
  }
  
  // Allocate a contiguous block of memory for the 2D array elements
  a[0] = malloc(n*m*sizeof(double));  // to store the first element of the 2D array 
  if (a[0] == NULL) {
    free(a);
    return NULL;
  }

  // Set the pointers in the array of double pointers
  // to point to the correct locations
  for (size_t i=1 ; i<m ; i++) {  // initialize array
    a[i] = a[i-1]+n;
  }

  return a;
}

void free_2d(double** a) {
  free(a[0]);  // free when not used anymore
  a[0] = NULL;  // set freed pointer to NULL
  free(a);
}


size_t read_Natoms(FILE* input_file){
    size_t Natoms;
    fscanf(input_file, "%zu", &Natoms);  //%zu for size_t values 
    return Natoms;
}


void read_molecule(FILE* input_file, size_t Natoms, double** coord, double* mass){ // void to define function that does not return a value 

	for (size_t i = 0; i < Natoms; i++) {
        	fscanf(input_file, "%lf %lf %lf %lf", &coord[i][0], &coord[i][1], &coord[i][2], &mass[i]);
	}
}
	

int main(){

    FILE* input_file = fopen("inp.txt", "r"); // Open input.txt for reading
    if (input_file == NULL) {
        printf("Error opening input file for reading");
        exit(-1);
    }

     size_t Natoms = read_Natoms(input_file);
     // Allocate memory for coordinates and masses
     double** coord = malloc_2d(Natoms, 3);  // 3 coordinates per atom
     double* mass = malloc(Natoms * sizeof(double));  // 1 mass per atom
     read_molecule(input_file, Natoms, coord, mass);

     // Print the read data for verification
    printf("Number fo atoms: %zu \n", Natoms);
    for (size_t i = 0; i < Natoms; i++) {
        printf("Atom %zu: x = %.1lf, y = %.1lf, z = %.1lf, mass = %.3lf\n",
               i + 1, coord[i][0], coord[i][1], coord[i][2], mass[i]);
    }

    // Free allocated memory
    free_2d(coord);
    free(mass);

    fclose(input_file);
    return 0;



}
