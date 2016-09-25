#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NUM_POINTS 524288 

unsigned int X_axis[NUM_POINTS];
unsigned int Y_axis[NUM_POINTS];

unsigned int X0_axis[NUM_POINTS];
unsigned int Y0_axis[NUM_POINTS];

int logfunc(int n)
{
  int a = 0;
  int b = 1;
  while (b < n) {
    b = b * 2;
    a++;
  }
  return a;
}

int expfunc(int n)
{
  int a;
  int b = 1;
  for (a = 0; a < n; a++) {
    b = b * 2;
  }
  return b;
}

float squareRoot(int x) {
  if(x < 0) return -1;
  float n = x;
  while(abs(n*n-x)>0.0001) {
    n=(n+x/n)/2;
  }
  return n;
}

float quadrant_cost(int *x, int *y, int num_quadrants, int myquadrant) {
  int i, j;
  float cost;
  for (i = 0; i < NUM_POINTS / num_quadrants; i++) {
    for (j = 1; j < NUM_POINTS / num_quadrants - i; j++) {
      cost += squareRoot((x[myquadrant + i] - x[myquadrant + i + j]) * (x[myquadrant + i] - x[myquadrant + i + j]) + 
                        (y[myquadrant + i] - y[myquadrant + i + j]) * (y[myquadrant + i] - y[myquadrant + i + j]));
    }
  }
  return cost;
}

/*
 * merge() 
 * Merge two sorted arrays, A1 with a integers and 
 * A2 with b integers, into a sorted array A.
 */
void merge(int *A1, int *B1, int a, int *A2, int *B2, int b, int *A, int *B) 
{
  int i,j,k;
  i = 0;
  j = 0;
  k = 0;
  while (i < a && j < b) {
    if (A1[i] <= A2[j]) {
      /* copy A[i] to C[k] and move the pointer i and k forward */
      A[k] = A1[i];
      B[k] = B1[i];
      i++;
      k++;
    }
    else {
      /* copy B[j] to C[k] and move the pointer j and k forward */
      A[k] = A2[j];
      B[k] = B2[j];
      j++;
      k++;
    }
  }
  /* move the remaining elements in A into C */
  while (i < a) {
    A[k] = A1[i];
    B[k] = B1[i];
    i++;
    k++;
  }
  /* move the remaining elements in B into C */
  while (j < b) {
    A[k] = A2[j];
    B[k] = B2[j];
    j++;
    k++;
  }
}

/*
 * merge_sort()
 * Sort array A with n integers, using merge-sort algorithm.
 */
void merge_sort1(int *A, int *B, int n)
{
  int i;
  int *A1, *A2, *B1, *B2;
  int n1, n2;

  if (n < 2) return;  /* the array is sorted when n=1.*/

  /* divide A into two array A1 and A2 */
  n1 = n / 2; /* the number of elements in A1 and B1 */
  n2 = n - n1;  /* the number of elements in A2 and B2 */

  A1 = (int*)malloc(sizeof(int) * n1);
  A2 = (int*)malloc(sizeof(int) * n2);
  B1 = (int*)malloc(sizeof(int) * n1);
  B2 = (int*)malloc(sizeof(int) * n2);

  /* move the first n/2 elements to A1 */
  for (i = 0; i < n1; i++) {
    A1[i] = A[i];
    B1[i] = B[i];
  }

  /* move the rest to A2 */
  for (i = 0; i < n2; i++) {
    A2[i] = A[i+n1];
    B2[i] = B[i+n1];
  }

  /* recursive call */
  merge_sort1(A1, B1, n1);
  merge_sort1(A2, B2, n2);

  /* conquer */
  merge(A1, B1, n1, A2, B2, n2, A, B);

  free(A1);
  free(A2);
  free(B1);
  free(B2);
}

void merge_sort2(int *A, int *B, int n, int unit)
{
  int i;
  int *A1, *A2, *B1, *B2;
  int n1, n2;

  if (n < unit + 1) return;

  /* divide A into two array A1 and A2 */
  n1 = n / 2; /* the number of elements in A1 and B1 */
  n2 = n - n1;  /* the number of elements in A2 and B2 */

  A1 = (int*)malloc(sizeof(int) * n1);
  A2 = (int*)malloc(sizeof(int) * n2);
  B1 = (int*)malloc(sizeof(int) * n1);
  B2 = (int*)malloc(sizeof(int) * n2);

  /* move the first n/2 elements to A1 */
  for (i = 0; i < n1; i++) {
    A1[i] = A[i];
    B1[i] = B[i];
  }

  /* move the rest to A2 */
  for (i = 0; i < n2; i++) {
    A2[i] = A[i+n1];
    B2[i] = B[i+n1];
  }

  /* recursive call */
  merge_sort2(A1, B1, n1, unit);
  merge_sort2(A2, B2, n2, unit);

  /* conquer */
  merge(A1, B1, n1, A2, B2, n2, A, B);

  free(A1);
  free(A2);
  free(B1);
  free(B2);
}

void find_quadrants(int num_quadrants, int myid, int numprocs)
{
  /* YOU NEED TO FILL IN HERE */
  int i, j, k;
  int dir;
  int *coordinates, *coordinates1;

  float local_cost;
  float total_cost;
  int myquadrant;

  coordinates = (int*)malloc(sizeof(int) * (4 * num_quadrants));
  coordinates1 = (int*)malloc(sizeof(int) * (4 * num_quadrants));
  unsigned int localX[NUM_POINTS/numprocs];
  unsigned int localY[NUM_POINTS/numprocs];

  coordinates[0] = X_axis[0];
  coordinates[1] = X_axis[0];
  coordinates[2] = Y_axis[0];
  coordinates[3] = Y_axis[0];

  for (i = 0; i < NUM_POINTS; i++) {
    if (X_axis[i] > coordinates[1]) {
      coordinates[1] = X_axis[i];
    }
    else if (X_axis[i] < coordinates[0]) {
      coordinates[0] = X_axis[i];
    }

    if (Y_axis[i] > coordinates[3]) {
      coordinates[3] = Y_axis[i];
    }
    else if (Y_axis[i] < coordinates[2]) {
      coordinates[2] = Y_axis[i];
    }
  }

  for (i = 0; i < logfunc(num_quadrants); i++) {
    for (j = 0; j < expfunc(i); j++) {
      for (k = 0; k < NUM_POINTS / ((expfunc(i)) * numprocs); k++) {
        localX[k] = X_axis[k + j * (NUM_POINTS / (expfunc(i))) + myid * (NUM_POINTS / ((expfunc(i)) * numprocs))];
        localY[k] = Y_axis[k + j * (NUM_POINTS / (expfunc(i))) + myid * (NUM_POINTS / ((expfunc(i)) * numprocs))];
      }

      if ((coordinates[4 * j + 1] - coordinates[4 * j]) >= 
        (coordinates[4 * j + 3] - coordinates[4 * j + 2]))
        dir = 1;
      else dir = 0;

      if (dir == 1) merge_sort1(localX, localY, NUM_POINTS / ((expfunc(i)) * numprocs)); 
      else merge_sort1(localY, localX, NUM_POINTS / ((expfunc(i)) * numprocs));

      MPI_Gather(&localX, NUM_POINTS / ((expfunc(i)) * numprocs), MPI_INT, X0_axis, NUM_POINTS / ((expfunc(i)) * numprocs), MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Gather(&localY, NUM_POINTS / ((expfunc(i)) * numprocs), MPI_INT, Y0_axis, NUM_POINTS / ((expfunc(i)) * numprocs), MPI_INT, 0, MPI_COMM_WORLD);

      if (myid == 0) {
        if (dir == 1) {
          merge_sort2(X0_axis, Y0_axis, NUM_POINTS / (expfunc(i)), NUM_POINTS / ((expfunc(i)) * numprocs));
          coordinates1[4 * 2 * j + 2] = coordinates[4 * j + 2];
          coordinates1[4 * 2 * j + 3] = coordinates[4 * j + 3];
          coordinates1[4 * (2 * j + 1) + 2] = coordinates[4 * j + 2];
          coordinates1[4 * (2 * j + 1) + 3] = coordinates[4 * j + 3];
          coordinates1[4 * 2 * j] = coordinates[4 * j];
          coordinates1[4 * 2 * j + 1] = X0_axis[NUM_POINTS / ((expfunc(i)) * 2)];
          coordinates1[4 * (2 * j + 1)] = X0_axis[NUM_POINTS / ((expfunc(i)) * 2)];
          coordinates1[4 * (2 * j + 1) + 1] = coordinates[4 * j + 1]; 
        }
        else {
          merge_sort2(Y0_axis, X0_axis, NUM_POINTS / (expfunc(i)), NUM_POINTS / ((expfunc(i)) * numprocs)); 
          coordinates1[4 * 2 * j] = coordinates[4 * j];
          coordinates1[4 * 2 * j + 1] = coordinates[4 * j + 1];
          coordinates1[4 * (2 * j + 1)] = coordinates[4 * j];
          coordinates1[4 * (2 * j + 1) + 1] = coordinates[4 * j + 1];
          coordinates1[4 * 2 * j + 2] = coordinates[4 * j + 2];
          coordinates1[4 * 2 * j + 3] = Y0_axis[NUM_POINTS / ((expfunc(i)) * 2)];
          coordinates1[4 * (2 * j + 1) + 2] = Y0_axis[NUM_POINTS / ((expfunc(i)) * 2)];
          coordinates1[4 * (2 * j + 1) + 3] = coordinates[4 * j + 3]; 
        }

        for (k = 0; k < NUM_POINTS / (expfunc(i)); k++) {
          X_axis[k + j * (NUM_POINTS / (expfunc(i)))] = X0_axis[k];
          Y_axis[k + j * (NUM_POINTS / (expfunc(i)))] = Y0_axis[k];
        }
      }
    }

    if (myid == 0) {
      for (k = 0; k < 4 * expfunc(i + 1); k ++) {
        coordinates[k] = coordinates1[k];
      }
    }
    MPI_Bcast(coordinates, 4 * num_quadrants, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&X_axis, NUM_POINTS, MPI_INT, 0, MPI_COMM_WORLD);  
    MPI_Bcast(&Y_axis, NUM_POINTS, MPI_INT, 0, MPI_COMM_WORLD);
  }

  if (myid == 0){
  total_cost = 0.0;

  //if (numprocs <= num_quadrants) {
    for (i = 0; i < num_quadrants; i++) {
      myquadrant = i * (NUM_POINTS / num_quadrants);
      total_cost += quadrant_cost(X_axis, Y_axis, num_quadrants, myquadrant);
    }
}
  //}
  //else {
    //if (myid < num_quadrants) {
      //myquadrant = myid * (NUM_POINTS / num_quadrants);
      //local_cost = quadrant_cost(X_axis, Y_axis, num_quadrants, myquadrant);
    //}
  //}

  //MPI_Reduce(&local_cost, &total_cost, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (myid == 0) {
    for (i = 0; i < num_quadrants; i++) {
      printf("The coordinates of quadrant %d are: x_axis[%d, %d], y_axis[%d, %d] \n", i,
       coordinates[4 * i], coordinates[4 * i + 1], coordinates[4 * i + 2], coordinates[4 * i + 3]);
    }
    printf("\nThe total partition cost of the quadrant distribution is %f\n", total_cost);
  }

}

int main(int argc, char *argv[])
{
  int num_quadrants;
  int myid, numprocs;
  int namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  double startwtime, endwtime;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Get_processor_name(processor_name,&namelen);

  if (argc != 2)
  {
    fprintf(stderr, "Usage: recursive_bisection <#of quadrants>\n");
    MPI_Finalize();
    exit(0);
  }

  fprintf (stderr,"Process %d on %s\n", myid, processor_name);

  num_quadrants = atoi(argv[1]);

  if (myid == 0)
    fprintf (stdout, "Extracting %d quadrants with %d processors \n", num_quadrants, numprocs);

  if (myid == 0)
  {
    int i;

    srand(10000);

    for (i = 0; i < NUM_POINTS; i++)
      X_axis[i] = (unsigned int)rand();

    for (i = 0; i < NUM_POINTS; i++)
      Y_axis[i] = (unsigned int)rand();

  }

  MPI_Bcast(&X_axis, NUM_POINTS, MPI_INT, 0, MPI_COMM_WORLD);  
  MPI_Bcast(&Y_axis, NUM_POINTS, MPI_INT, 0, MPI_COMM_WORLD);  

  if (myid == 0) startwtime = MPI_Wtime();
  find_quadrants(num_quadrants, myid, numprocs);
  if (myid == 0) {
    endwtime = MPI_Wtime();
    printf("\nWall clock time = %f s.\n", endwtime - startwtime);
  }
  MPI_Finalize();
  return 0;
}
