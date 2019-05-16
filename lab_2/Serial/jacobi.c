#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/timeb.h>

void gen_laplace_mat(int n, int **diag, int **l_upper_band, int **l_lower_band, int **u_upper_band, int **u_lower_band);
void calc_newx(int, int, int, int, int *, int *, int *, int *, int *, double *, double *, double *);
void init_x(int, double *);
void init_true_x(int, double *);
void init_b(int, int, int, int, int *, int *, int *, int *, int *, double *, double *);
double calc_resid(int, int, int, int, int *, int *, int *, int *, int *, double *, double *);

int main(int argc, char *argv[])
{
  int n, i, diag_size, mid_size, short_size, execs = 0;
  int *diag, *l_upper, *l_lower, *u_upper, *u_lower;
  double *x, *newx, *b, *true_x;
  double residual = 10000;
  struct timeb start, stop;
  int time_diff;
  printf("argc: %u %d \n", argc, argv[1]);
  if (argc != 2)
  {
    printf("You loose args\n");
    return 1;
  }

  n = atoi(argv[1]);

  ftime(&start);

  gen_laplace_mat(n, &diag, &l_upper, &l_lower, &u_upper, &u_lower);

  diag_size = n * n;
  mid_size = n * n - 1;
  short_size = n * n - n;

  x = malloc(diag_size * sizeof(double));
  newx = malloc(diag_size * sizeof(double));
  true_x = malloc(diag_size * sizeof(double));
  b = malloc(diag_size * sizeof(double));

  init_x(diag_size, x);
  init_true_x(diag_size, true_x);
  init_b(n, diag_size, mid_size, short_size, diag, l_upper, l_lower, u_upper, u_lower, true_x, b);

  while (residual > 0.001 && execs < 100000)
  {
    calc_newx(n, diag_size, mid_size, short_size, diag, l_upper, l_lower, u_upper, u_lower, x, b, newx);
    for (i = 0; i < diag_size; i++)
    {
      x[i] = newx[i];
    }
    residual = calc_resid(n, diag_size, mid_size, short_size, diag, l_upper, l_lower, u_upper, u_lower, x, b);
    execs++;
  }

  // end timing here
  ftime(&stop);
  time_diff = (int)(1000.0 * (stop.time - start.time) + (stop.millitm - start.millitm));
  printf("Finished in %u milliseconds.\n", time_diff);
  printf("In %d iterations, got residual of %f\n", execs, residual);

  return 0;
}

void calc_newx(int n, int diag_size, int mid_size, int short_size,
               int *diag, int *l_upper, int *l_lower, int *u_upper, int *u_lower,
               double *x, double *b, double *newx)
{
  int i;

  for (i = 0; i < diag_size; i++) {
    newx[i] = ((i < short_size ? u_upper[i] * x[i + n] : 0) +
               (i < mid_size ? u_lower[i] * x[i + 1] : 0) +
               (i > 0 ? l_upper[i - 1] * x[i - 1] : 0) +
               (i >= n ? l_lower[i - n] * x[i - n] : 0)) /
              (diag[i]);
    newx[i] = b[i] / (diag[i]) - newx[i];
  }
}

void init_x(int size, double *x) {
  int i;
  for (i = 0; i < size; i++)
  {
    x[i] = (i % 2) == 0 ? 2 : -2;
  }
  return;
}

void init_true_x(int size, double *true_x) {
  int i;
  for (i = 0; i < size; i++)
  {
    true_x[i] = 1;
  }
  return;
}

void init_b(int n, int diag_size, int mid_size, int short_size,
            int *diag, int *l_upper, int *l_lower, int *u_upper, int *u_lower,
            double *true_x, double *b)
{
  int i;
  for (i = 0; i < diag_size; i++)
  {
    // get b by multiplying Ax, or (U + L + D)x
    b[i] = ((i < short_size ? u_upper[i] * true_x[i + n] : 0) +
            (i < mid_size ? u_lower[i] * true_x[i + 1] : 0) +
            (i > 0 ? l_upper[i - 1] * true_x[i - 1] : 0) +
            (i >= n ? l_lower[i - n] * true_x[i - n] : 0) +
            (true_x[i] * diag[i]));
  }
}

double calc_resid(
    int n, int diag_size, int mid_size, int short_size,
    int *diag, int *l_upper, int *l_lower, int *u_upper, int *u_lower,
    double *x, double *b)
{
  int i;
  double sum = 0, r;
  for (i = 0; i < diag_size; i++)
  {
    // calculate Ax - b for each row
    r = ((i < short_size ? u_upper[i] * x[i + n] : 0) +
         (i < mid_size ? u_lower[i] * x[i + 1] : 0) +
         (i > 0 ? l_upper[i - 1] * x[i - 1] : 0) +
         (i >= n ? l_lower[i - n] * x[i - n] : 0) +
         (x[i] * diag[i])) -
        b[i];
    sum += r * r;
  }
  return sum * sum;
}

void stencil(int, int, int *, int *);

// Generates a 2D laplacian matrix that is
// n^2 x n^2
// O(N) algorithm
void gen_laplace_mat(
    int n,
    int **diag,
    int **l_upper_band,
    int **l_lower_band,
    int **u_upper_band,
    int **u_lower_band) {

  int diag_len = n * n;
  int mid_band_len = diag_len - 1;
  int outer_band_len = diag_len - n;
  int i, low_fill, up_fill;

  *diag = malloc(diag_len * sizeof(int));
  *l_upper_band = malloc(mid_band_len * sizeof(int));
  *u_lower_band = malloc(mid_band_len * sizeof(int));
  *l_lower_band = malloc(outer_band_len * sizeof(int));
  *u_upper_band = malloc(outer_band_len * sizeof(int));

  for (i = 0; i < diag_len; i++) {
    (*diag)[i] = -4;

    stencil(n, i + 1, &low_fill, &up_fill);
    if (i > 0)
      (*l_upper_band)[i - 1] = low_fill;
    if (i < mid_band_len)
      (*u_lower_band)[i] = up_fill;

    if (i < outer_band_len) {
      (*l_lower_band)[i] = 1;
      (*u_upper_band)[i] = 1;
    }
  }
  return;
}

void stencil(int n, int center_pos, int *lower_fill, int *upper_fill) {
  *lower_fill = center_pos % n == 1 ? 0 : 1;
  *upper_fill = center_pos % n == 0 ? 0 : 1;

  return;
}
