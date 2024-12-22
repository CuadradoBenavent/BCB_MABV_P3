#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double** malloc_2d(size_t m, size_t n) {
  double** a = malloc(m*sizeof(double*)); // double pointers of size m
  if (a == NULL) {
    return NULL;
  }

  a[0] = malloc(n*m*sizeof(double));  //contiguous block of memory 
  if (a[0] == NULL) {
    free(a);
    return NULL;
  }

  for (size_t i=1 ; i<m ; i++) {  
    a[i] = a[i-1]+n;  // for pointers of array to point to correct location
  }
  return a;
  }

void free_2d(double** a) {
  free(a[0]);
  a[0] = NULL;
  free(a);
}



int main(){
    size_t m = 3, n = 4;
    // Allocate a 3x4 matrix
    double** matrix = malloc_2d(m, n);

    if (matrix == NULL) {
        printf("Memory allocation failed!\n");
        return 1;
    }

    // Set values to test the matrix
    // Filling the matrix with some values for testing
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            matrix[i][j] = (i * n) + j;  // Just an example of filling values
        }
    }

    // Print the matrix
    printf("3x4 Matrix:\n");
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            printf("%f ", matrix[i][j]);  // Print the values in the matrix
        }
        printf("\n");
    }

    // Free the allocated memory
    free_2d(matrix);

    return 0;



}
