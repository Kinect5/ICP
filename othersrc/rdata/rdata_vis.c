//
//   Copyright (C) 2009 ViGIR Lab (http://vigir.missouri.edu)
//   Written by Guilherme N. DeSouza <desouzag@missouri.edu>
//
//
/* 
 * June 2001, Johnny Park
 * 
 * June 2003,
 *
*/
/*
 * Range data programs
 */

#include "rdata_vis.h"


/* 
 * Reads range data
 */
int ReadRdataVis(char *name, rdata_vis *rd)
{
  register int  i, j, k;
  FILE          *fp;
  char          filename[256];

  // read. xyz file
  sprintf(filename, "%s.xyz", name);
  if ( (fp = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "Cannot read %s\n", filename);
    return 0;
  }
  if (fread(&(rd->num_pt), sizeof(int), 1, fp) != 1) {
    fprintf(stderr, "IO error in %s\n", filename);
    return 0;
  }
  printf("\tNumber of points: %d\n", rd->num_pt); fflush(stdout);
  rd->xyz = (float *) malloc (3 * (rd->num_pt) * sizeof(float));
  for (j=0; j<rd->num_pt; j++) {
    if (fread(&(rd->xyz[3*j]), sizeof(float), 3, fp) != 3) {
      fprintf(stderr, "IO error in %s\n", filename);
      return 0;
    }
  }
  fclose(fp);

  // read .nor file
  sprintf(filename, "%s.nor", name);
  if ( (fp = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "Cannot read %s\n", filename);
    return 0;
  }
  if (fread(&(rd->num_pt), sizeof(int), 1, fp) != 1) {
    fprintf(stderr, "IO error in %s\n", filename);
    return 0;
  }
  printf("\tNumber of points: %d\n", rd->num_pt); fflush(stdout);
  rd->nor = (float *) malloc (3 * (rd->num_pt) * sizeof(float));
  for (j=0; j<rd->num_pt; j++) {
    if (fread(&(rd->nor[3*j]), sizeof(float), 3, fp) != 3) {
      fprintf(stderr, "IO error in %s\n", filename);
      return 0;
    }
  }
  fclose(fp);

  // read .tri file
  sprintf(filename, "%s.tri", name);
  if ( (fp = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "Cannot read %s\n", filename);
    return 0;
  }
  if (fread(&(rd->num_tri), sizeof(int), 1, fp) != 1) {
    fprintf(stderr, "IO error in %s\n", filename);
    return 0;
  }
  printf("\tNumber of triangles: %d\n", rd->num_tri); fflush(stdout);
  rd->tri = (int *) malloc (3 * (rd->num_tri) * sizeof(int));
  for (j=0; j<rd->num_tri; j++) {
    if (fread(&(rd->tri[3*j]), sizeof(int), 3, fp) != 3) {
      fprintf(stderr, "IO error in %s\n", filename);
      return 0;
    }
  }
  fclose(fp);

  // read .rgb file if available
  sprintf(filename, "%s.rgb", name);
  if ( (fp = fopen(filename, "rb")) != NULL) {
    printf("\t.rgb available\n"); fflush(stdout);
    if (fread(&(rd->num_pt), sizeof(int), 1, fp) != 1) {
      fprintf(stderr, "IO error in %s\n", filename);
      return 0;
    }
    rd->rgb = (float *) malloc (3 * (rd->num_pt) * sizeof(float));
    for (j=0; j<rd->num_pt; j++) {
      if (fread(&(rd->rgb[3*j]), sizeof(float), 3, fp) != 3) {
        fprintf(stderr, "IO error in %s\n", filename);
        return 0;
      }
    }
    rd->flag_rgb = 1; 
    fclose(fp);
  }
  else {
    printf("\t.rgb NOT available\n"); fflush(stdout);
    rd->flag_rgb = 0;
  }
   
  // read .Rt file
  sprintf(filename, "%s_new.Rt", name);
  if ( (fp = fopen(filename, "r")) == NULL) {
    sprintf(filename, "%s.Rt", name);
    if ( (fp = fopen(filename, "r")) == NULL) {
      fprintf(stderr, "Cannot read %s\n", filename);
      return 0;
    }
  }
  for (j=0; j<3; j++) for (k=0; k<3; k++) {
    if (fscanf(fp, "%lf", &(rd->M[4*k+j])) != 1) {
      fprintf(stderr, "IO error in %s\n", filename);
      return 0;
    }
  }
  for (j=0; j<3; j++) {
    if (fscanf(fp, "%lf", &(rd->M[12+j])) != 1) {
      fprintf(stderr, "IO error in %s\n", filename);
      return 0;
    }
  }
  rd->M[15] = 1.0;
  fclose(fp);

  // assign name
  strcpy(rd->filename, name);
  i = strlen(name);
  k = 0;
  for (j=i-1; j>=0; j--) {
    if (name[j] == '/') break;
    k++;
  }
  strncpy(rd->name, &name[j+1], k);

  return 1;


}
