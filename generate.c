#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float uniform_rand() {
  return ((float)rand()/(float)(RAND_MAX + 1.0f));
}

int main(int argc, char** argv) {
  if (argc != 4) {
    printf("Usage: generator [atoms per side] [number of lipids] [data filename]\n");
    exit(1);
  }

  int atoms_per_side = atoi(argv[1]);
  int molecules = atoi(argv[2]);

  float grid = 1.19523f;
  float box = grid*(float)atoms_per_side;
  float bond_length = 0.75f;

  int grid_atoms = atoms_per_side*atoms_per_side;
  int natoms = grid_atoms + 2*molecules;

  float *x = malloc(natoms*sizeof(float));
  float *y = malloc(natoms*sizeof(float));
  int *h = malloc(natoms*sizeof(float));
  int *t = malloc(natoms*sizeof(int));
  int *m = malloc(natoms*2*sizeof(int));

  for (int i = 0; i < atoms_per_side; i++) {
    for (int j = 0; j < atoms_per_side; j++) {
      x[i*atoms_per_side + j] = grid*j;
      y[i*atoms_per_side + j] = grid*i;
    }
  }

  for (int i = 0; i < grid_atoms; i++) {
    h[i] = i;
    t[i] = 1;
  }

  for (int i = 0; i < molecules; i++) {
    int r = (int) (uniform_rand()*(grid_atoms - i));
    r += i;
    int tmp = h[r];
    h[r] = h[i];
    h[i] = tmp;
  }

  for (int i = 0; i < molecules; i++) {
    t[h[i]] = 2;
    int lipid_i = grid_atoms + 2*i;
    t[lipid_i] = 3;
    t[lipid_i + 1] = 4;

    float lipid_angle = uniform_rand()*2.0*3.14159;
    float lax = cos(lipid_angle);
    float lay = sin(lipid_angle);

    x[lipid_i] = x[h[i]] + bond_length*lax;
    y[lipid_i] = y[h[i]] + bond_length*lay;
    x[lipid_i + 1] = x[h[i]] + 2.0*bond_length*lax;
    y[lipid_i + 1] = y[h[i]] + 2.0*bond_length*lay;
  }

  for (int i = 0; i < natoms; i++) {
    if (x[i] < 0.0) x[i] += box;
    if (x[i] >= box) x[i] -= box;
    if (y[i] < 0.0) y[i] += box;
    if (y[i] >= box) y[i] -= box;
  }

  FILE *data_file = fopen(argv[3], "w");

  if (data_file == NULL) {
    printf("Problem opening data file %s\n", argv[3]);
    exit(1);
  }

  fprintf(data_file, "\n\n");
  fprintf(data_file, "%12d atoms\n", natoms);
  fprintf(data_file, "%12d bonds\n", 2*molecules);
  fprintf(data_file, "\n\n\n\n\n\n\n\n\n\n");
  fprintf(data_file, "%4.5f %4.5f xlo xhi\n", 0.0, box);
  fprintf(data_file, "%4.5f %4.5f ylo yhi\n", 0.0, box);
  fprintf(data_file, "%4.5f %4.5f zlo zhi\n\n", -0.1, 0.1);

  fprintf(data_file, " Atoms\n\n");

  for (int i = 0; i < natoms; i++) {
    fprintf(data_file, "%7d %7d %7d %4.3f %4.3f %4.3f\n", i + 1, 0, t[i], x[i], y[i], 0.0);
  }

  fprintf(data_file, "\n Bonds\n\n");

  for (int i = 0; i < molecules; i++) {
    int lipid_i = grid_atoms + 2*i;
    fprintf(data_file, "%7d %7d %7d %7d\n", 2*i + 1, 1, h[i] + 1, lipid_i + 1);
    fprintf(data_file, "%7d %7d %7d %7d\n", 2*i + 2, 1, lipid_i + 1, lipid_i + 2);
  }

  fflush(data_file);
  fclose(data_file);

  free(x);
  free(y);
  free(h);
  free(t);
  free(m);
}
