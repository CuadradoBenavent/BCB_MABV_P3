#include <stdio.h>
#include <stdlib.h>
#include <math.h>


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

void free_2d(double** a) {
  free(a[0]);
  a[0] = NULL;
  free(a);
}


//Function to read the number of atoms from xyz file
size_t read_Natoms(FILE* input_file) {
  size_t Natoms;
  if (fscanf(input_file, "%zu", &Natoms) != 1) {
    fprintf(stderr, "Failed to read the number of atoms\n");
    exit(1);
  }
  return Natoms;
}


// Function to read the molecule data
void read_molecule(FILE* input_file, size_t Natoms, double** coord, double* mass) {
  for (size_t i = 0; i < Natoms; i++) {
    if (fscanf(input_file, "%lf %lf %lf %lf", &coord[i][0], &coord[i][1], &coord[i][2], &mass[i]) != 4) {
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



int main() {
  const char* filename = "inp.txt";
  FILE* file = fopen(filename, "r");
  if (file == NULL) {
    perror("Failed to open file");
    return 1;
  }

  //Read the number of atoms
  size_t Natoms = read_Natoms(file); 
  printf("Number of atoms: %zu\n", Natoms);


  //Allocate memory for coords and masses
  double** coord = malloc_2d(Natoms, 3);
  double* mass = malloc(Natoms * sizeof(double));
  if (coord == NULL || mass == NULL) {
    fprintf(stderr, "Memory allocation failed\n");
    fclose(file);
    return 1;
  }

  //Read molecule data
  read_molecule(file, Natoms, coord, mass);
  fclose(file);


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
  printf("\nInternuclear distances:\n");
  for (size_t i = 0; i < Natoms; i++) {
    for (size_t j = 0; j < Natoms; j++) {
      printf("Interatomic distance[%zu][%zu] = %lf\n", i, j, distance[i][j]);
    }
  }

  double epsilon = 0.0661;
  double sigma = 0.3345;
  double total_potential = V(epsilon, sigma, Natoms, distance);
  printf("Total potential energy: %lf J/mol\n", total_potential);





  //Print data
  printf("Molecule data from %s:\n", filename);
  for (size_t i = 0; i < Natoms; i++) {
    printf("Atom %zu; X=%lf Y=%lf Z=%lf mass=%lf \n", i + 1, coord[i][0], coord[i][1], coord[i][2], mass[i]);
  }



  // Free allocated memory
  free_2d(coord);
  free_2d(distance);
  free(mass);

  return 0;
}
