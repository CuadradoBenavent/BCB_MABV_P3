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


int main() {
  const char* filename = "inp.txt";
  FILE* file = fopen(filename, "r");
  if (file == NULL) {
    perror("Failed to open file");
    return 1;
  }


  //Read the number of rows
  size_t rows;
  if (fscanf(file, "%zu", &rows) != 1) {
    fprintf(stderr, "Failed to read the number of rows\n");
    fclose(file);
    return 1;
  }

  //Number of columns (fixed for this example)
  size_t cols = 4;

  //Allocate the 2D array
  double** array = malloc_2d(rows,cols);
  if (array == NULL) {
    fprintf(stderr, "Memory allocated failed\n");
    fclose(file);
    return 1;
  }


  //Read the data into the array
  for (size_t i=0; i < rows; i++) {
    if (fscanf(file, "%lf %lf %lf %lf", &array[i][0], &array[i][1], &array[i][2], &array[i][3]) != 4) {
      fprintf(stderr, "Failed to read data for row %zu\n", i);
      free_2d(array);
      fclose(file);
      return 1;
    }
  }

  fclose(file);


  //Print the data 
  printf("Read data from %s:\n", filename);
  for (size_t i=0; i < rows; i++) {
    printf("Row %zu: %lf %lf %lf %lf\n", i, array[i][0], array[i][1], array[i][2], array[i][3]);
  }


  // Free allocated memory
  free_2d(array);

  return 0;
}
