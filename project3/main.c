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
  const char* filename = "input.txt";
  FILE* file = fopen(filename, "r");
  if (file == NULL) {
    perror("Failed to open file");
