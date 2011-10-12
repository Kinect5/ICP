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

#ifndef _ICP_H_
#define _ICP_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "matrix.h"
#include "point3D.h"
#include "time_keeper.h"


#define   WriteClosestPts     0
#define   WriteRotAndTrans    1
#define   WriteIterationInfo  0


/*
 * linked list pointing members of a bin
 */
typedef struct member {  
  struct member *next;    // pointer to next member
  int    thispt;          // array index of corresponding point
} member;

typedef struct bin {
  struct member *first;
  struct member *current;
} bin;



/*
 * Function Declaration
 */
void ICPalgorithm(matrix *fR, vector *ft,
                  point_xyz *data1, int num_data1, 
                  point_xyz *data2, int num_data2,
                  int subsample, double approx_error,
                  int max_num_iteration, double low_bound_mean_error,
                  double low_bound_delta_error);
void WeightedICP(matrix *fR, vector *ft,
                 point_xyz *data1, double *weight1, int num_data1, 
                 point_xyz *data2, double *weight2, int num_data2,
                 int subsample, double approx_error,
                 int max_num_iteration, double low_bound_mean_error,
                 double low_bound_delta_error);
void SimultaneousMaxMin(void);
void CreateBins(void);
void ClosestPointViaElias(void);
void ClosestPoint(void);
void ComputeRotationAndTranslation(void);
void ApplyRotationAndTranslation(void);


#endif

