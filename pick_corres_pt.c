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

#include "pick_corres_pt.h"


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


void PickCorresAnc(double x1, double y1, double z1,
                   double x2, double y2, double z2)
{
  register int  i, j;
  int           corres_rd, corres_pt;
  double        lx, ly, lz;   // unit vector in the direction
                              // from (x2,y2,z2) to (x1,y1,z1)
  double        qx, qy, qz;
  double        px, py, pz;
  double        mag;
  double        d, min_dist;
  double        x, y, z;
  double        nx, ny, nz;

  // initialize values
  min_dist = 1.0;
  corres_rd = corres_pt = -1;

  lx = (double)(x1 - x2); 
  ly = (double)(y1 - y2);
  lz = (double)(z1 - z2);
  mag = sqrt(lx*lx + ly*ly + lz*lz);
  lx /= mag; ly /= mag; lz /= mag;

  for (i=0; i<num_rdata_anc; i++) {
    if (!rd_anc[i].flag_display) continue;
    for (j=0; j<rd_anc[i].num_pt; j++) {

      // transform by its modeling transformation matrix
      x =  rd_anc[i].M[0] * rd_anc[i].xyz[3*j] +
           rd_anc[i].M[4] * rd_anc[i].xyz[3*j+1] +
           rd_anc[i].M[8] * rd_anc[i].xyz[3*j+2] + rd_anc[i].M[12];
      y =  rd_anc[i].M[1] * rd_anc[i].xyz[3*j] +
           rd_anc[i].M[5] * rd_anc[i].xyz[3*j+1] +
           rd_anc[i].M[9] * rd_anc[i].xyz[3*j+2] + rd_anc[i].M[13];
      z =  rd_anc[i].M[2] * rd_anc[i].xyz[3*j] +
           rd_anc[i].M[6] * rd_anc[i].xyz[3*j+1] +
           rd_anc[i].M[10]* rd_anc[i].xyz[3*j+2] + rd_anc[i].M[14];
      nx = rd_anc[i].M[0] * rd_anc[i].nor[3*j] +
           rd_anc[i].M[4] * rd_anc[i].nor[3*j+1] +
           rd_anc[i].M[8] * rd_anc[i].nor[3*j+2];
      ny = rd_anc[i].M[1] * rd_anc[i].nor[3*j] +
           rd_anc[i].M[5] * rd_anc[i].nor[3*j+1] +
           rd_anc[i].M[9] * rd_anc[i].nor[3*j+2];
      nz = rd_anc[i].M[2] * rd_anc[i].nor[3*j] +
           rd_anc[i].M[6] * rd_anc[i].nor[3*j+1] +
           rd_anc[i].M[10]* rd_anc[i].nor[3*j+2];

      // if dot product is negative, skip it
      if ( (lx*nx + ly*ny + lz*nz) < 0) {
        continue;
      }
      // compute distance between each point to the line
      qx = x - x1;
      qy = y - y1;
      qz = z - z1;
      px = ly*qz - lz*qy;
      py = lz*qx - lx*qz;
      pz = lx*qy - ly*qx;
      d = sqrt(px*px + py*py + pz*pz); 
      if (d < min_dist) {
        min_dist = d;
        corres_rd = i; 
        corres_pt = j;
      }
    }
  } 
  if (min_dist < 1.0) {
    corres_rd_anc[num_corres_anc] = corres_rd;
    corres_pt_anc[num_corres_anc] = corres_pt;
    num_corres_anc++;
  }
  //printf("%f %f %f\n", lx, ly, lz);
  //printf("%f\n", min_dist);
}


void PickCorresMov(double x1, double y1, double z1,
                   double x2, double y2, double z2)
{
  register int  i, j;
  int           corres_rd, corres_pt;
  double        lx, ly, lz;   // unit vector in the direction
                              // from (x2,y2,z2) to (x1,y1,z1)
  double        qx, qy, qz;
  double        px, py, pz;
  double        mag;
  double        d, min_dist;
  double        x, y, z;
  double        nx, ny, nz;


  // initialize values
  min_dist = 1.0;
  corres_rd = corres_pt = -1;

  lx = (double)(x1 - x2); 
  ly = (double)(y1 - y2);
  lz = (double)(z1 - z2);
  mag = sqrt(lx*lx + ly*ly + lz*lz);
  lx /= mag; ly /= mag; lz /= mag;

  for (i=0; i<num_rdata_mov; i++) {
    if (!rd_mov[i].flag_display) continue;
    for (j=0; j<rd_mov[i].num_pt; j++) {
      // transform by its modeling transformation matrix
      x =  rd_mov[i].M[0] * rd_mov[i].xyz[3*j] +
           rd_mov[i].M[4] * rd_mov[i].xyz[3*j+1] +
           rd_mov[i].M[8] * rd_mov[i].xyz[3*j+2] + rd_mov[i].M[12];
      y =  rd_mov[i].M[1] * rd_mov[i].xyz[3*j] +
           rd_mov[i].M[5] * rd_mov[i].xyz[3*j+1] +
           rd_mov[i].M[9] * rd_mov[i].xyz[3*j+2] + rd_mov[i].M[13];
      z =  rd_mov[i].M[2] * rd_mov[i].xyz[3*j] +
           rd_mov[i].M[6] * rd_mov[i].xyz[3*j+1] +
           rd_mov[i].M[10]* rd_mov[i].xyz[3*j+2] + rd_mov[i].M[14];
      nx = rd_mov[i].M[0] * rd_mov[i].nor[3*j] +
           rd_mov[i].M[4] * rd_mov[i].nor[3*j+1] +
           rd_mov[i].M[8] * rd_mov[i].nor[3*j+2];
      ny = rd_mov[i].M[1] * rd_mov[i].nor[3*j] +
           rd_mov[i].M[5] * rd_mov[i].nor[3*j+1] +
           rd_mov[i].M[9] * rd_mov[i].nor[3*j+2];
      nz = rd_mov[i].M[2] * rd_mov[i].nor[3*j] +
           rd_mov[i].M[6] * rd_mov[i].nor[3*j+1] +
           rd_mov[i].M[10]* rd_mov[i].nor[3*j+2];

      // if dot product is negative, skip it
      if ( (lx*nx + ly*ny + lz*nz) < 0) {
        continue;
      }
      // compute distance between each point to the line
      qx = x - x1;
      qy = y - y1;
      qz = z - z1;
      px = ly*qz - lz*qy;
      py = lz*qx - lx*qz;
      pz = lx*qy - ly*qx;
      d = sqrt(px*px + py*py + pz*pz); 
      if (d < min_dist) {
        min_dist = d;
        corres_rd = i; 
        corres_pt = j;
      }
    }
  } 
  if (min_dist < 1.0) {
    corres_rd_mov[num_corres_mov] = corres_rd;
    corres_pt_mov[num_corres_mov] = corres_pt;
    num_corres_mov++;
  }
  //printf("%f %f %f\n", lx, ly, lz);
  //printf("%f\n", min_dist);
}


