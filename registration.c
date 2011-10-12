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

#include "registration.h"


extern int        window_main_anc,
                  window_main_mov;
extern rdata_vis  *rd_anc;            // array of "anchor" range data 
extern rdata_vis  *rd_mov;            // array of "moving" range data
extern int        *corres_rd_anc,
                  *corres_rd_mov,
                  *corres_pt_anc,
                  *corres_pt_mov;
extern int        num_rdata_anc,
                  num_rdata_mov,
                  num_corres_anc,
                  num_corres_mov;
extern double     new_M[16];

/*
 * Peform a course registration using selected corresponding points by user
 */
void CourseRegistration(void)
{
  register int  i;
  double        Sxx, Sxy, Sxz,
                Syx, Syy, Syz,
                Szx, Szy, Szz; 
  double        max_eval;
  point_xyz     *p, *q;
  point_xyz     mean_corres_p,
                mean_corres_q;  
  int           num_corres; 
  matrix        *Q, *R;
  vector        *max_evec, *t;

  // initialize values
  Sxx = Sxy = Sxz = Syx = Syy = Syz = Szx = Szy = Szz = 0.0;
  mean_corres_p.x = mean_corres_p.y = mean_corres_p.z = 0.0;
  mean_corres_q.x = mean_corres_q.y = mean_corres_q.z = 0.0;

  // take the smaller one
  num_corres = (num_corres_anc<num_corres_mov)?num_corres_anc:num_corres_mov;

  // allocate memory
  Q = AllocateMatrix(4, 4);
  R = AllocateMatrix(3, 3);
  max_evec = AllocateVector(4);
  t = AllocateVector(3);
  p = (point_xyz *) malloc (num_corres * sizeof(point_xyz));
  q = (point_xyz *) malloc (num_corres * sizeof(point_xyz));
  
  // transform
  for (i=0; i<num_corres; i++) {
    p[i].x = (rd_mov[corres_rd_mov[i]].xyz[3*corres_pt_mov[i]] *
              rd_mov[corres_rd_mov[i]].M[0]) +
             (rd_mov[corres_rd_mov[i]].xyz[3*corres_pt_mov[i]+1] *
              rd_mov[corres_rd_mov[i]].M[4]) +
             (rd_mov[corres_rd_mov[i]].xyz[3*corres_pt_mov[i]+2] *
              rd_mov[corres_rd_mov[i]].M[8]) +
              rd_mov[corres_rd_mov[i]].M[12];
    p[i].y = (rd_mov[corres_rd_mov[i]].xyz[3*corres_pt_mov[i]] *
              rd_mov[corres_rd_mov[i]].M[1]) +
             (rd_mov[corres_rd_mov[i]].xyz[3*corres_pt_mov[i]+1] *
              rd_mov[corres_rd_mov[i]].M[5]) +
             (rd_mov[corres_rd_mov[i]].xyz[3*corres_pt_mov[i]+2] *
              rd_mov[corres_rd_mov[i]].M[9]) +
              rd_mov[corres_rd_mov[i]].M[13];
    p[i].z = (rd_mov[corres_rd_mov[i]].xyz[3*corres_pt_mov[i]] *
              rd_mov[corres_rd_mov[i]].M[2]) +
             (rd_mov[corres_rd_mov[i]].xyz[3*corres_pt_mov[i]+1] *
              rd_mov[corres_rd_mov[i]].M[6]) +
             (rd_mov[corres_rd_mov[i]].xyz[3*corres_pt_mov[i]+2] *
              rd_mov[corres_rd_mov[i]].M[10]) +
              rd_mov[corres_rd_mov[i]].M[14];

    q[i].x = (rd_anc[corres_rd_anc[i]].xyz[3*corres_pt_anc[i]] *
              rd_anc[corres_rd_anc[i]].M[0]) +
             (rd_anc[corres_rd_anc[i]].xyz[3*corres_pt_anc[i]+1] *
              rd_anc[corres_rd_anc[i]].M[4]) +
             (rd_anc[corres_rd_anc[i]].xyz[3*corres_pt_anc[i]+2] *
              rd_anc[corres_rd_anc[i]].M[8]) +
              rd_anc[corres_rd_anc[i]].M[12];
    q[i].y = (rd_anc[corres_rd_anc[i]].xyz[3*corres_pt_anc[i]] *
              rd_anc[corres_rd_anc[i]].M[1]) +
             (rd_anc[corres_rd_anc[i]].xyz[3*corres_pt_anc[i]+1] *
              rd_anc[corres_rd_anc[i]].M[5]) +
             (rd_anc[corres_rd_anc[i]].xyz[3*corres_pt_anc[i]+2] *
              rd_anc[corres_rd_anc[i]].M[9]) +
              rd_anc[corres_rd_anc[i]].M[13];
    q[i].z = (rd_anc[corres_rd_anc[i]].xyz[3*corres_pt_anc[i]] *
              rd_anc[corres_rd_anc[i]].M[2]) +
             (rd_anc[corres_rd_anc[i]].xyz[3*corres_pt_anc[i]+1] *
              rd_anc[corres_rd_anc[i]].M[6]) +
             (rd_anc[corres_rd_anc[i]].xyz[3*corres_pt_anc[i]+2] *
              rd_anc[corres_rd_anc[i]].M[10]) +
              rd_anc[corres_rd_anc[i]].M[14];

    //printf("%d: %f %f %f\n", i, q[i].x, q[i].y, q[i].z);
  }
 
  
  for (i=0; i<num_corres; i++) {
    mean_corres_p.x += p[i].x;
    mean_corres_p.y += p[i].y;
    mean_corres_p.z += p[i].z;
      
    mean_corres_q.x += q[i].x;
    mean_corres_q.y += q[i].y;
    mean_corres_q.z += q[i].z;
  }
  mean_corres_p.x /= (double)num_corres;
  mean_corres_p.y /= (double)num_corres;
  mean_corres_p.z /= (double)num_corres;
  mean_corres_q.x /= (double)num_corres;
  mean_corres_q.y /= (double)num_corres;
  mean_corres_q.z /= (double)num_corres;

     
  for (i=0; i<num_corres; i++) {
    Sxx += p[i].x * q[i].x;
    Sxy += p[i].x * q[i].y;
    Sxz += p[i].x * q[i].z;
    Syx += p[i].y * q[i].x;
    Syy += p[i].y * q[i].y;
    Syz += p[i].y * q[i].z;
    Szx += p[i].z * q[i].x;
    Szy += p[i].z * q[i].y;
    Szz += p[i].z * q[i].z;
  }

  Sxx = Sxx / (double)num_corres - (mean_corres_p.x * mean_corres_q.x);
  Sxy = Sxy / (double)num_corres - (mean_corres_p.x * mean_corres_q.y);
  Sxz = Sxz / (double)num_corres - (mean_corres_p.x * mean_corres_q.z);
  Syx = Syx / (double)num_corres - (mean_corres_p.y * mean_corres_q.x);
  Syy = Syy / (double)num_corres - (mean_corres_p.y * mean_corres_q.y);
  Syz = Syz / (double)num_corres - (mean_corres_p.y * mean_corres_q.z);
  Szx = Szx / (double)num_corres - (mean_corres_p.z * mean_corres_q.x);
  Szy = Szy / (double)num_corres - (mean_corres_p.z * mean_corres_q.y);
  Szz = Szz / (double)num_corres - (mean_corres_p.z * mean_corres_q.z);

  // construct N
  Q->entry[0][0] = Sxx + Syy + Szz;
  Q->entry[1][0] = Q->entry[0][1] = Syz - Szy;
  Q->entry[2][0] = Q->entry[0][2] = Szx - Sxz;
  Q->entry[3][0] = Q->entry[0][3] = Sxy - Syx;
  Q->entry[1][1] = Sxx - Syy - Szz;
  Q->entry[1][2] = Q->entry[2][1] = Sxy + Syx;
  Q->entry[1][3] = Q->entry[3][1] = Szx + Sxz;
  Q->entry[2][2] = -Sxx + Syy - Szz;
  Q->entry[2][3] = Q->entry[3][2] = Syz + Szy;
  Q->entry[3][3] = -Sxx - Syy + Szz;

  // --- compute largest eigenvalues and eigenvectors of Q ---
  SymmetricLargestEigens(Q, max_evec, &max_eval); 
  // make sure max_evec[0] > 0
  if (max_evec->entry[0] < 0) {
    for (i=0; i<4; i++) max_evec->entry[i] *= -1.0;
  }
  // --- compute rotation matrix ---
  RotationQuaternion(max_evec, R);
 
  // --- compute translation vector ---
  t->entry[0] = mean_corres_q.x - 
                R->entry[0][0] * mean_corres_p.x -
                R->entry[0][1] * mean_corres_p.y -
                R->entry[0][2] * mean_corres_p.z;
  t->entry[1] = mean_corres_q.y - 
                R->entry[1][0] * mean_corres_p.x -
                R->entry[1][1] * mean_corres_p.y -
                R->entry[1][2] * mean_corres_p.z;
  t->entry[2] = mean_corres_q.z - 
                R->entry[2][0] * mean_corres_p.x -
                R->entry[2][1] * mean_corres_p.y -
                R->entry[2][2] * mean_corres_p.z;

  //PrintMatrix(R);
  //PrintVector(t);

  new_M[0] = R->entry[0][0]; new_M[4] = R->entry[0][1]; new_M[8] = R->entry[0][2];
  new_M[1] = R->entry[1][0]; new_M[5] = R->entry[1][1]; new_M[9] = R->entry[1][2];
  new_M[2] = R->entry[2][0]; new_M[6] = R->entry[2][1]; new_M[10] = R->entry[2][2];
  new_M[12] = t->entry[0];   new_M[13] = t->entry[1];   new_M[14] = t->entry[2];
  new_M[3] = new_M[7] = new_M[11] = 0; new_M[15] = 1;

  // free memory
  FreeMatrix(Q); FreeMatrix(R);
  FreeVector(max_evec); FreeVector(t);
  free(p); free(q);
}


/*
 * Perfrom a fine registration using the ICP algorithm
 */
void FineRegistration(void)
{
  int         i, j, count;
  int         num_p, num_q;
  point_xyz   *p, *q, cp;
  matrix      *R;
  vector      *t;
  double      Mf[16];

  // compute total number of points
  num_p = num_q = 0;
  for (i=0; i<num_rdata_anc; i++) num_q += rd_anc[i].num_pt;
  for (i=0; i<num_rdata_mov; i++) num_p += rd_mov[i].num_pt;

  // allocate memory
  R = AllocateMatrix(3, 3);
  t = AllocateVector(3);
  p = (point_xyz *) malloc (num_p * sizeof(point_xyz));
  q = (point_xyz *) malloc (num_q * sizeof(point_xyz));

  // copy xyz values
  count = 0;
  for (i=0; i<num_rdata_anc; i++) {
    
    for (j=0; j<rd_anc[i].num_pt; j++) {
      cp.x = rd_anc[i].xyz[3*j];
      cp.y = rd_anc[i].xyz[3*j+1];
      cp.z = rd_anc[i].xyz[3*j+2];
      // transform with its modeling transformation matrix
      q[count].x = cp.x * rd_anc[i].M[0] + cp.y * rd_anc[i].M[4] +
                   cp.z * rd_anc[i].M[8] + rd_anc[i].M[12];
      q[count].y = cp.x * rd_anc[i].M[1] + cp.y * rd_anc[i].M[5] +
                   cp.z * rd_anc[i].M[9] + rd_anc[i].M[13];
      q[count].z = cp.x * rd_anc[i].M[2] + cp.y * rd_anc[i].M[6] +
                   cp.z * rd_anc[i].M[10] + rd_anc[i].M[14]; 
      count++;
    }
  }
  count = 0;
  for (i=0; i<num_rdata_mov; i++) {
    // use modelview just for matrix multiplication
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMultMatrixd(new_M);
    glMultMatrixd(rd_mov[i].M);
    glGetDoublev(GL_MODELVIEW_MATRIX, Mf);
    glPopMatrix();

    for (j=0; j<rd_mov[i].num_pt; j++) {
      // transform 
      cp.x = rd_mov[i].xyz[3*j];
      cp.y = rd_mov[i].xyz[3*j+1];
      cp.z = rd_mov[i].xyz[3*j+2];
      
      p[count].x = cp.x*Mf[0] + cp.y*Mf[4] + cp.z*Mf[8] + Mf[12];
      p[count].y = cp.x*Mf[1] + cp.y*Mf[5] + cp.z*Mf[9] + Mf[13];
      p[count].z = cp.x*Mf[2] + cp.y*Mf[6] + cp.z*Mf[10] + Mf[14]; 
      count++;
    }
  }

  // perform ICP
  //ICPalgorithm(R, t, p, num_p, q, num_q, 5, 5.0, 100, 0.10, 0.0005);
  ICPalgorithm(R, t, p, num_p, q, num_q, 5, 5.0, 1000, 0.010, 0.000005);
  Mf[0] = R->entry[0][0]; Mf[4] = R->entry[0][1]; Mf[8] = R->entry[0][2];
  Mf[1] = R->entry[1][0]; Mf[5] = R->entry[1][1]; Mf[9] = R->entry[1][2];
  Mf[2] = R->entry[2][0]; Mf[6] = R->entry[2][1]; Mf[10] = R->entry[2][2];
  Mf[12] = t->entry[0];   Mf[13] = t->entry[1];   Mf[14] = t->entry[2];
  Mf[3] = Mf[7] = Mf[11] = 0; Mf[15] = 1;

  // update new_M
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glMultMatrixd(Mf);
  glMultMatrixd(new_M);
  glGetDoublev(GL_MODELVIEW_MATRIX, new_M);
  glPopMatrix();

  // free memory
  FreeMatrix(R); FreeVector(t);
  free(p); free(q);
}


/*
 * Save registration result
 */
void SaveRegistration(void)
{
  FILE *fp;
  char outfile[256];
  int i;

  glutSetWindow(window_main_mov);
  glMatrixMode(GL_MODELVIEW);
  for (i=0; i<num_rdata_mov; i++) {
    // update modeling matrix
    glPushMatrix();
    glLoadIdentity();
    glMultMatrixd(new_M);
    glMultMatrixd(rd_mov[i].M);
    glGetDoublev(GL_MODELVIEW_MATRIX, rd_mov[i].M);
    glPopMatrix();
    
    // write to .Rt file
    sprintf(outfile, "%s.Rt", rd_mov[i].filename);
    fp = fopen(outfile, "w");
    fprintf(fp, "%f %f %f\n", rd_mov[i].M[0], rd_mov[i].M[4], rd_mov[i].M[8]);
    fprintf(fp, "%f %f %f\n", rd_mov[i].M[1], rd_mov[i].M[5], rd_mov[i].M[9]);
    fprintf(fp, "%f %f %f\n", rd_mov[i].M[2], rd_mov[i].M[6], rd_mov[i].M[10]);
    fprintf(fp, "%f %f %f\n", rd_mov[i].M[12], rd_mov[i].M[13], 
                              rd_mov[i].M[14]);
    fclose(fp);
  }
  // set new_M to a identity matrix
  for (i=0; i<16; i++) new_M[i] = 0;
  new_M[0] = new_M[5] = new_M[10] = new_M[15] = 1;
}






