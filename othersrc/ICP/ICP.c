//
//   Copyright (C) 2009 ViGIR Lab (http://vigir.missouri.edu)
//   Written by Guilherme N. DeSouza <desouzag@missouri.edu>
//
//
/*
 * ICP algorithm
 * 
 * June 2001, Johnny Park
 * 
 * June 2003,
 *   Added weighted ICP algorithm
 */
 
#include "ICP.h"

/*
 * Global Variables
 */
point_xyz *p,             /* pointer to array of data set p */
          *q,             /* pointer to array of data set q */
          *p_new;         /* pointer to array of data set p_new */
double *wp;               /* weight of data set p */
double *wq;               /* weight of data set q */
int weight_flag = 0;      /* 1: incorporate weight, 0: without weight */
int num_p,                /* number of points in data set p */ 
    num_q;                /* number of points in data set q */

bin ***qbin;              /* 3D array of bins */
int num_bin;              /* number of bins along each axis */

double max_x, min_x,      /* max and min values of x, y, z */
       max_y, min_y, 
       max_z, min_z;
double binsize_x,         /* size of bins along each axis */
       binsize_y, 
       binsize_z;

int *closest_pt;          /* closest_pt[i] = j :                     */
                          /*   p[i] and q[j] are the closest points  */
                          /*   if j = -1, then no point exists       */
int num_closest_pt;       /* number of corresponding closest point   */

double Dmax,              /* max distance allowed to be closest points at each ICP iteration */
       zeta,              /* absolute max distance allowed to be closest points */
       D,                 /* distance which represent a very good match */
       mean_error,        /* mean of error between closest points */ 
       std_error;         /* standard deviation of error between closest points */
int    sub_sample;        /* integer value for subsampling point set */

matrix *cp,               /* cp[i] (a point in p) */
       *cq_t,             /* transpose of cq[i] (a point in q that is */ 
                          /* closest to cp[i]) */
       *mean_cp,          /* mean of cp[] */
       *mean_cq_t,        /* mean of cq[]' */
       *cpcq_t,           /* cp[i] * cq[i]' */
       *sum_cpcq_t,       /* sum of cp[i] * cq[i]' */
       *mean,             /* mean_cpi * mean_cqi_t */
       *sigma,            /* cross covariance matrix */
       *sigma_t,          /* transpose of sigma */
       *A,                /* anti-symmetric matrix of sigma */
       *B,                /* sigma - sigma_t - tr(sigma)*I3 */
       *Q;                /* 4x4 matrix Q of (25) in Besl & McKay paper */
double max_eval;          /* largest eigenvalue of Q */
vector *max_evec;         /* eigenvector corresponding to max eigenvalue */

matrix *R;                /* rotation matrix */
vector *t;                /* translation vector */


FILE *fp_rot_trans,
     *fp_iteration;



/*
 * ICPalgorithm
 * 
 * Performs the ICP between two sets of point p and q
 * finds the rotation matrix R and the translation vector t such that 
 * q = Rp + t
 * 
 * The iteration terminates if at least one of the following 
 * three conditions is satisfied
 *  (1) Number of iteration exceeded max_num_iteration
 *  (2) Overall mean error is less than low_bound_mean_error
 *  (3) Delta error is less than low_bound_delta_error
 * 
 * Note that after termination of this function, data1 will be 
 * transformed points.
 */
void ICPalgorithm(matrix *fR, vector *ft,
                  point_xyz *data1, int num_data1, 
                  point_xyz *data2, int num_data2,
                  int subsample, double approx_error,
                  int max_num_iteration, double low_bound_mean_error,
                  double low_bound_delta_error)
{
  register int i, j;
  double previous_error,    /* error in the previous iteration */ 
         delta_error;       /* current error - previous error */
  int iteration_count;

  p_new = data1;           q = data2;
  num_p = num_data1;   num_q = num_data2;

  /* copy p */
  p = (point_xyz *) malloc (num_p * sizeof(point_xyz));
  for (i=0; i<num_p; i++) {
    p[i].x = p_new[i].x;
    p[i].y = p_new[i].y;
    p[i].z = p_new[i].z;
  }
  
  printf("num_p: %d, num_q: %d\n\n", num_p, num_q);  fflush(stdout);
  StartTime();

  zeta = approx_error;
  D = low_bound_mean_error;
  sub_sample = subsample;
  
  /* find max and min values of (x,y,z)*/ 
  SimultaneousMaxMin();

  /* create bins */
  CreateBins();
  
  /* start new file*/
  if (WriteRotAndTrans) {
    fp_rot_trans = fopen("rot_trans", "w+");
    fclose(fp_rot_trans);
    fp_rot_trans = fopen("rot_trans", "a");
  }
  if (WriteIterationInfo) {
    fp_iteration = fopen("iteration_info", "w+");
    fclose(fp_iteration);
    fp_iteration = fopen("iteration_info", "a");
  }

  /* allocate variables */
  closest_pt = (int *) malloc (num_p * sizeof(int));
  R = AllocateMatrix(3, 3);
  t = AllocateVector(3);
  cp = AllocateMatrix(3, 1); 
  cq_t = AllocateMatrix(1, 3);
  cpcq_t = AllocateMatrix(3, 3);
  sum_cpcq_t = AllocateMatrix(3, 3);
  mean = AllocateMatrix(3, 3);
  mean_cp = AllocateMatrix(3, 1);
  mean_cq_t = AllocateMatrix(1, 3);
  Q = AllocateMatrix(4, 4);
  sigma = AllocateMatrix(3, 3);
  sigma_t = AllocateMatrix(3, 3);
  A = AllocateMatrix(3, 3);
  B = AllocateMatrix(3, 3);
  max_evec = AllocateVector(4);

  /* initializie some variables */
  previous_error = 99999;
  mean_error = 99999;
  delta_error = 99999;
  iteration_count = 1;
  
  printf("Starting ICP...\n");  
  printf("\tmax_num_iteration: %d\n", max_num_iteration);
  printf("\tlow_bound_mean_error: %lf\n", low_bound_mean_error);
  printf("\tlow_bound_delta_error: %lf\n", low_bound_delta_error);
  fflush(stdout);
  
  /* --------------------------------------------------------------- *
   *                      Start of ICP iteration                     *
   * --------------------------------------------------------------- */                       
  while ( (iteration_count < max_num_iteration) &&
          (mean_error > low_bound_mean_error) &&
          (delta_error > low_bound_delta_error)    ) {
    printf("iteration %d\n", iteration_count);
    iteration_count++;
    
    /* for each point in p_new[], find closest point in q[] */
    ClosestPointViaElias();
    printf("\tnum_closest_pt:\t%d\n", num_closest_pt);
    printf("\terror:\t%f\n", mean_error);
    printf("\tstd_error:\t%f\n", std_error);
    printf("\tDmax:\t%f\n", Dmax);
    
    /* update error values */
    delta_error = fabs(mean_error - previous_error);
    previous_error = mean_error;
    printf("\tdelta_error: %f\n", delta_error); fflush(stdout);

    /* write iteration info */
    if (WriteIterationInfo) {
      fprintf(fp_iteration, "iteration %d\n", iteration_count-1);
      fprintf(fp_iteration, "\tnum_closest_pt:\t%d\n", num_closest_pt);
      fprintf(fp_iteration, "\terror:\t%f\n", mean_error);
      fprintf(fp_iteration, "\tstd_error:\t%f\n", std_error);
      fprintf(fp_iteration, "\tDmax:\t%f\n", Dmax);
      fprintf(fp_iteration, "\tdelta_error: %f\n\n", delta_error);
    }
    
    /* compute rotation matrix and traslation vector */
    ComputeRotationAndTranslation();

    /* update p_new by applying rotation and translation */
    ApplyRotationAndTranslation();

    /* output */
    if (WriteRotAndTrans) {
      for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
          fprintf(fp_rot_trans, "%f  ", R->entry[i][j]);
        }
        fprintf(fp_rot_trans, "\n");
      }
      fprintf(fp_rot_trans, "\n");
      for (i=0; i<3; i++) fprintf(fp_rot_trans, "%f  ", t->entry[i]);
      fprintf(fp_rot_trans, "\n\n");
    }
    
    PrintTime();    
  }
  /* --------------------------------------------------------------- *
   *                      End of ICP iteration                       *
   * --------------------------------------------------------------- */ 
  
  printf("\tR = \n"); PrintMatrix(R);
  printf("\tt = \n"); PrintVector(t);

  /* copy final R and t */
  CopyMatrix(fR, R);  CopyVector(ft, t);
  
  if (WriteRotAndTrans) fclose(fp_rot_trans);
  if (WriteIterationInfo) fclose(fp_iteration);

  // free memory
  free(p);
  free(closest_pt); 
  for (i=0; i<num_bin; i++) {
    for (j=0; j<num_bin; j++) {
      free(qbin[i][j]);
    }
  }
  for (i=0; i<num_bin; i++) {
    free(qbin[i]);
  }
  free(qbin); 
  FreeMatrix(R); FreeVector(t); FreeMatrix(cp); FreeMatrix(cq_t);
  FreeMatrix(cpcq_t); FreeMatrix(sum_cpcq_t); FreeMatrix(mean);
  FreeMatrix(mean_cp); FreeMatrix(mean_cq_t); FreeMatrix(Q); 
  FreeMatrix(sigma); FreeMatrix(sigma_t); FreeMatrix(A); FreeMatrix(B);
  FreeVector(max_evec);
}



void WeightedICP(matrix *fR, vector *ft,
                 point_xyz *data1, double *weight1, int num_data1, 
                 point_xyz *data2, double *weight2, int num_data2,
                 int subsample, double approx_error,
                 int max_num_iteration, double low_bound_mean_error,
                 double low_bound_delta_error)
{
  weight_flag = 1;

  wp = weight1;
  wq = weight2;

  ICPalgorithm(fR, ft, data1, num_data1, data2, num_data2, subsample,
               approx_error, max_num_iteration, low_bound_mean_error,
               low_bound_delta_error);
}



/*
 * Simultaneously finds max and min values of (x,y,z)
 * Requires 3(n/2) comparisons as opposed to naive 2n-2 comparisons
 */
void SimultaneousMaxMin(void)
{
  int i;

  printf("finding max and min of (x,y,z) values...\n"); fflush(stdout);
  
  /* initialize values */
  max_x = -9999999; max_y = -9999999; max_z = -9999999;
  min_x = 9999999;  min_y = 9999999;  min_z = 9999999;

  for (i=0; i<floor(num_q/2); i++) {
    if (q[2*i].x > q[2*i+1].x) {
      if (max_x < q[2*i].x)    max_x = q[2*i].x;
      if (min_x > q[2*i+1].x)  min_x = q[2*i+1].x;
    }
    else {
      if (max_x < q[2*i+1].x)  max_x = q[2*i+1].x;
      if (min_x > q[2*i].x)    min_x = q[2*i].x;
    }
    
    if (q[2*i].y > q[2*i+1].y) {
      if (max_y < q[2*i].y)    max_y = q[2*i].y;
      if (min_y > q[2*i+1].y)  min_y = q[2*i+1].y;
    }
    else {
      if (max_y < q[2*i+1].y)  max_y = q[2*i+1].y;
      if (min_y > q[2*i].y)    min_y = q[2*i].y;
    }

    if (q[2*i].z > q[2*i+1].z) {
      if (max_z < q[2*i].z)    max_z = q[2*i].z;
      if (min_z > q[2*i+1].z)  min_z = q[2*i+1].z;
    }
    else {
      if (max_z < q[2*i+1].z)  max_z = q[2*i+1].z;
      if (min_z > q[2*i].z)    min_z = q[2*i].z;
    }
  }

  /* check the last element in case num_p is an odd number*/ 
  if (max_x < q[num_q-1].x)  max_x = q[num_q-1].x;
  if (min_x > q[num_q-1].x)  min_x = q[num_q-1].x;
  if (max_y < q[num_q-1].y)  max_y = q[num_q-1].y;
  if (min_y > q[num_q-1].y)  min_y = q[num_q-1].y;
  if (max_z < q[num_q-1].z)  max_z = q[num_q-1].z;
  if (min_z > q[num_q-1].z)  min_z = q[num_q-1].z;

  printf("\tx: (%f, %f)\n\ty: (%f, %f)\n\tz: (%f, %f)\n\n",
         min_x, max_x, min_y, max_y, min_z, max_z);
  PrintTime();  fflush(stdout);
}



/*
 * Creates bins to store all the points
 * Number of bins and its size are determined by max & min of 
 * point set p, and num_p
 * According to "Acceleration of Binning Nearest Neightbour Method",
 * given 10^5 random points in a cubic, 50 bins along each axis gave
 * a optimum result.  So, we'll follow that.
 */
void CreateBins(void)
{
  int i, j;
  int bin_i, bin_j, bin_k;      /* index of bins */
  double min_halfbin_x,         /* min_halfbin_x = min_x - (binsize_x/2) */ 
        min_halfbin_y, 
        min_halfbin_z;          

  printf("Creating bins...\n"); fflush(stdout);
  
  /* calculate number of bins along each axis */
  //num_bin = (int)(num_q / 2000);
  //num_bin = (int)(num_q / 500);
  //if (num_bin < 5) num_bin = 5;
  num_bin = 50;

  /* calculate size of bins along each axis */   
  binsize_x = (double)(max_x - min_x) / (double)num_bin;
  binsize_y = (double)(max_y - min_y) / (double)num_bin;
  binsize_z = (double)(max_z - min_z) / (double)num_bin;

  min_halfbin_x = min_x - (binsize_x/2.0);
  min_halfbin_y = min_y - (binsize_y/2.0);
  min_halfbin_z = min_z - (binsize_z/2.0);

  printf("\tzeta: %f\n", zeta);
  printf("\tD: %f\n", D);

  /* initialize Dmax to zeta */
  Dmax = zeta;
                       
  num_bin++;
  printf("\tnum_bin: %d * %d * %d\n", num_bin, num_bin, num_bin); 
  fflush(stdout);
    
  /* allocate bins */
  qbin = (bin ***) malloc (sizeof(bin **) * num_bin);
  for (i=0; i<num_bin; i++) {
    qbin[i] = (bin **) malloc (sizeof(bin *) * num_bin);
    for (j=0; j<num_bin; j++) {
      qbin[i][j] = (bin *) calloc (num_bin, sizeof(bin));
    }
  }
  /* put all the points of p into a bin */
  for (i=0; i<num_q; i++) {
    bin_i = (int)((q[i].x - min_halfbin_x) / binsize_x);
    bin_j = (int)((q[i].y - min_halfbin_y) / binsize_y);
    bin_k = (int)((q[i].z - min_halfbin_z) / binsize_z);

    /* first member */
    if (qbin[bin_i][bin_j][bin_k].first == NULL) {
      qbin[bin_i][bin_j][bin_k].first = (member *) calloc (1, sizeof(member));
      qbin[bin_i][bin_j][bin_k].current = qbin[bin_i][bin_j][bin_k].first;
    }
    /* not a first member */
    else {
      qbin[bin_i][bin_j][bin_k].current->next = 
                                        (member *) calloc (1, sizeof(member));
      qbin[bin_i][bin_j][bin_k].current = 
                                      qbin[bin_i][bin_j][bin_k].current->next;
    }
    
    qbin[bin_i][bin_j][bin_k].current->thispt = i;
  }
  
  PrintTime();  fflush(stdout);
}



/*
 * find closest point using Elias method
 * returns mean square point matching error
 */
void ClosestPointViaElias(void)
{
  int i, j, k, m;
  int bin_i, bin_j, bin_k;  /* index of bins */
  int start_i, end_i,       /* start and end of neighbor index */
      start_j, end_j,
      start_k, end_k;
  double min_halfbin_x,     /* min_halfbin_x = min_x - (binsize_x/2) */
         min_halfbin_y, 
         min_halfbin_z;
  double max_halfbin_x,     /* max_halfbin_x = max_x + (binsize_z/2) */
         max_halfbin_y,
         max_halfbin_z;

  int current_q;
  int min_distance_pt;
  double min_distance;
  double init_min_distance;
  double current_distance;
  double dist_pq[num_p];    /* distance between closest points */
  FILE *fp_c_pts;
  int search_binsize;

  //printf("Finding closest points via Elias method...\n"); fflush(stdout);

  /* compute some variables */
  min_halfbin_x = min_x - (binsize_x/2.0);
  min_halfbin_y = min_y - (binsize_y/2.0);
  min_halfbin_z = min_z - (binsize_z/2.0);
  max_halfbin_x = max_x + (binsize_x/2.0);
  max_halfbin_y = max_y + (binsize_y/2.0);
  max_halfbin_z = max_z + (binsize_z/2.0);
  init_min_distance = sqrt(SQ(max_x - min_x) + 
                           SQ(max_y - min_y) + 
                           SQ(max_z - min_z)   ) + 1;

  /* initialize some values */
  num_closest_pt = 0;
  mean_error = 0;
  for (m=0; m<num_p; m++) 
    closest_pt[m] = -1;

  /* find closest point */                    
  for (m=0; m<num_p; m = m + (int) rint(((double) sub_sample) * (rand() / (RAND_MAX + 1.0))) ) {   //see manpage for rand()
  //for (m=0; m<num_p; m = m + (rand() % sub_sample) ) {  
  //for (m=0; m<num_p; m++) {
    /* set the min_distance to the largest value possible */
	  min_distance = init_min_distance;
    min_distance_pt = -1;
	  
    /* p_new[m] is out of bin boundary */
    if ( (-p_new[m].x + min_halfbin_x > Dmax) || 
         ( p_new[m].x - max_halfbin_x > Dmax) ||
         (-p_new[m].y + min_halfbin_y > Dmax) || 
         ( p_new[m].y - max_halfbin_y > Dmax) ||
         (-p_new[m].z + min_halfbin_z > Dmax) || 
         ( p_new[m].z - max_halfbin_z > Dmax)   ) {
      ;  // no closest point exists
    }
   
    
    /* p_new[m] is within bin boundary */
    else {
	    /* find corresponding bin for each point p_new[m] */
      bin_i = (int)((p_new[m].x - min_halfbin_x) / binsize_x);
      bin_j = (int)((p_new[m].y - min_halfbin_y) / binsize_y);
      bin_k = (int)((p_new[m].z - min_halfbin_z) / binsize_z);
  
      // find boundary of neighbor bins according to Dmax
      search_binsize = (int)ceil(Dmax / binsize_x);
      start_i = (bin_i-search_binsize < 0) ? 0 : bin_i-search_binsize;
      end_i = (bin_i+search_binsize > num_bin-1) ? 
                         num_bin-1 : bin_i+search_binsize;
      
      search_binsize = (int)ceil(Dmax / binsize_y);
      start_j = (bin_j-search_binsize < 0) ? 0 : bin_j-search_binsize;
      end_j = (bin_j+search_binsize > num_bin-1) ? 
                         num_bin-1 : bin_j+search_binsize;
      
      search_binsize = (int)ceil(Dmax / binsize_z);
      start_k = (bin_k-search_binsize < 0) ? 0 : bin_k-search_binsize;
      end_k = (bin_k+search_binsize > num_bin-1) ? 
                         num_bin-1 : bin_k+search_binsize;       
      
      for (i = start_i; i <= end_i; i++)
      for (j = start_j; j <= end_j; j++)
      for (k = start_k; k <= end_k; k++) {
	      /* at least one member is in this bin */
	      if (qbin[i][j][k].first != NULL) {
	        /* set a pointer to the first member of this bin */
	        qbin[i][j][k].current = qbin[i][j][k].first;
	        
	        /* find the member with the closest E-distance */
	        while (qbin[i][j][k].current != NULL) {
            current_q = qbin[i][j][k].current->thispt;

	          current_distance = E_distance(&p_new[m], &q[current_q]);
	          if (current_distance < min_distance) {
	            min_distance = current_distance;
	            min_distance_pt = current_q;
	          }
	          qbin[i][j][k].current = qbin[i][j][k].current->next;
	          //printf("\tcurrent_distance:%f\n", current_distance);
	          //printf("\tmin_distance:%f\n", min_distance);
	        }
	      }
	    }
	    /* is min_distance less than Dmax? */
	    if (min_distance < Dmax) {
	      /* store the array index of the found closest point */
	      closest_pt[m] = min_distance_pt;

	      /* update mean_error and distance_pq */
	      mean_error += min_distance;
	      dist_pq[num_closest_pt] = min_distance;

	      num_closest_pt++;
	    }	     
	  }
  }

  if (num_closest_pt == 0) {
    printf("No closest points exists!!! exiting...\n");
    exit(-1);
  }
  
  /* compute mean_error */
  mean_error /= (double)num_closest_pt;

  /* compute std_error */
  std_error = 0;
  for (i=0; i<num_closest_pt; i++) {
    std_error += (dist_pq[i] - mean_error) * (dist_pq[i] - mean_error);
  }
  std_error /= (double)num_closest_pt;
  std_error = sqrt(std_error);


/*
  // compute Dmax (Zhang paper) 
  if (mean_error < D)            // the registration is quite good 
    //Dmax = 3.0*std_error + mean_error;
    Dmax = 2.0*std_error + mean_error;
  else if (mean_error < 3.0*D)   // the registration is still good 
    //Dmax = 2.0*std_error + mean_error;
    Dmax = 1.5*std_error + mean_error;
  else if (mean_error < 6.0*D)   // the registration is not too good 
    Dmax = std_error + mean_error;
  else                           // the registration is really bad 
    Dmax = zeta;
 
*/
  Dmax = 3.0*std_error + mean_error;

  //printf("\tnum_closest_pt: %d\n", num_closest_pt);
  //printf("\tmean square error: %lf\n", error);
  //PrintTime();  fflush(stdout);

  /* output closest points */
  if (WriteClosestPts) {
    fp_c_pts = fopen("closest_pts", "w+");
    for (i=0; i<num_p; i++) {
      if ( closest_pt[i] > -1 ) {
        fprintf(fp_c_pts, "%d %d\n", i+1, closest_pt[i]+1); /* for MATLAB index */
      }
    }
    fclose(fp_c_pts);
  }
}


/*
 * compute rotation matrix and translation vector
 */
void ComputeRotationAndTranslation(void)
{
  register int  i;
  double        Sxx, Sxy, Sxz,
                Syx, Syy, Syz,
                Szx, Szy, Szz; 
  double        max_eval;
  double        wp_sum, wq_sum;
  double        wpq, wpq_sum;
  point_xyz     mean_corres_p,
                mean_corres_q;        


  // initialize values
  Sxx = Sxy = Sxz = Syx = Syy = Syz = Szx = Szy = Szz = 0.0;
  mean_corres_p.x = mean_corres_p.y = mean_corres_p.z = 0.0;
  mean_corres_q.x = mean_corres_q.y = mean_corres_q.z = 0.0;
  wp_sum = wq_sum = wpq_sum = 0.0;


  /*
   * With weight
   */
  if (weight_flag) {
    for (i=0; i<num_p; i++) {
      if (closest_pt[i] > -1) {

        wpq = wp[i] * wq[closest_pt[i]];
      
        mean_corres_p.x += wpq * p[i].x;
        mean_corres_p.y += wpq * p[i].y;
        mean_corres_p.z += wpq * p[i].z;
      
        mean_corres_q.x += wpq * q[closest_pt[i]].x;
        mean_corres_q.y += wpq * q[closest_pt[i]].y;
        mean_corres_q.z += wpq * q[closest_pt[i]].z;

        wpq_sum += wpq;
      }
    }
    mean_corres_p.x /= wpq_sum;
    mean_corres_p.y /= wpq_sum;
    mean_corres_p.z /= wpq_sum;
    mean_corres_q.x /= wpq_sum;
    mean_corres_q.y /= wpq_sum;
    mean_corres_q.z /= wpq_sum;

    for (i=0; i<num_p; i++) {
      if (closest_pt[i] > -1) {
    
        wpq = wp[i] * wq[closest_pt[i]];
        
        Sxx += wpq * (p[i].x             - mean_corres_p.x) * 
                     (q[closest_pt[i]].x - mean_corres_q.x);
        Sxy += wpq * (p[i].x             - mean_corres_p.x) * 
                     (q[closest_pt[i]].y - mean_corres_q.y);
        Sxz += wpq * (p[i].x             - mean_corres_p.x) * 
                     (q[closest_pt[i]].z - mean_corres_q.z);
        Syx += wpq * (p[i].y             - mean_corres_p.y) * 
                     (q[closest_pt[i]].x - mean_corres_q.x);
        Syy += wpq * (p[i].y             - mean_corres_p.y) * 
                     (q[closest_pt[i]].y - mean_corres_q.y);
        Syz += wpq * (p[i].y             - mean_corres_p.y) * 
                     (q[closest_pt[i]].z - mean_corres_q.z);
        Szx += wpq * (p[i].z             - mean_corres_p.z) * 
                     (q[closest_pt[i]].x - mean_corres_q.x);
        Szy += wpq * (p[i].z             - mean_corres_p.z) * 
                     (q[closest_pt[i]].y - mean_corres_q.y);
        Szz += wpq * (p[i].z             - mean_corres_p.z) * 
                     (q[closest_pt[i]].z - mean_corres_q.z);
      }
    }
  }


  /*
   * Without weight
   */  
  else {
    for (i=0; i<num_p; i++) {
      if (closest_pt[i] > -1) {
        mean_corres_p.x += p[i].x;
        mean_corres_p.y += p[i].y;
        mean_corres_p.z += p[i].z;
      
        mean_corres_q.x += q[closest_pt[i]].x;
        mean_corres_q.y += q[closest_pt[i]].y;
        mean_corres_q.z += q[closest_pt[i]].z;

        Sxx += (p[i].x) * (q[closest_pt[i]].x);
        Sxy += (p[i].x) * (q[closest_pt[i]].y);
        Sxz += (p[i].x) * (q[closest_pt[i]].z);
        Syx += (p[i].y) * (q[closest_pt[i]].x);
        Syy += (p[i].y) * (q[closest_pt[i]].y);
        Syz += (p[i].y) * (q[closest_pt[i]].z);
        Szx += (p[i].z) * (q[closest_pt[i]].x);
        Szy += (p[i].z) * (q[closest_pt[i]].y);
        Szz += (p[i].z) * (q[closest_pt[i]].z);
      }
    }
    mean_corres_p.x /= (double)num_closest_pt;
    mean_corres_p.y /= (double)num_closest_pt;
    mean_corres_p.z /= (double)num_closest_pt;
    mean_corres_q.x /= (double)num_closest_pt;
    mean_corres_q.y /= (double)num_closest_pt;
    mean_corres_q.z /= (double)num_closest_pt;

    Sxx = Sxx / (double)num_closest_pt - (mean_corres_p.x * mean_corres_q.x);
    Sxy = Sxy / (double)num_closest_pt - (mean_corres_p.x * mean_corres_q.y);
    Sxz = Sxz / (double)num_closest_pt - (mean_corres_p.x * mean_corres_q.z);
    Syx = Syx / (double)num_closest_pt - (mean_corres_p.y * mean_corres_q.x);
    Syy = Syy / (double)num_closest_pt - (mean_corres_p.y * mean_corres_q.y);
    Syz = Syz / (double)num_closest_pt - (mean_corres_p.y * mean_corres_q.z);
    Szx = Szx / (double)num_closest_pt - (mean_corres_p.z * mean_corres_q.x);
    Szy = Szy / (double)num_closest_pt - (mean_corres_p.z * mean_corres_q.y);
    Szz = Szz / (double)num_closest_pt - (mean_corres_p.z * mean_corres_q.z);
  }
/*
  printf("%lf %lf %lf\n", mean_corres_p.x, mean_corres_p.y, mean_corres_p.z);
  printf("%lf %lf %lf\n", mean_corres_q.x, mean_corres_q.y, mean_corres_q.z);
  printf("%lf %lf %lf\n", Sxx, Sxy, Sxz);
  printf("%lf %lf %lf\n", Syx, Syy, Syz);
  printf("%lf %lf %lf\n", Szx, Szy, Szz);
  exit(0);
*/

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
}
  


/*
 * update p_new by applying rotation and translation
 */
void ApplyRotationAndTranslation()
{
  int i;
 
  for (i=0; i<num_p; i++) {
    p_new[i].x = (R->entry[0][0] * p[i].x) +
                 (R->entry[0][1] * p[i].y) +
                 (R->entry[0][2] * p[i].z) + t->entry[0];

    p_new[i].y = (R->entry[1][0] * p[i].x) +
                 (R->entry[1][1] * p[i].y) +
                 (R->entry[1][2] * p[i].z) + t->entry[1];

    p_new[i].z = (R->entry[2][0] * p[i].x) +
                 (R->entry[2][1] * p[i].y) +
                 (R->entry[2][2] * p[i].z) + t->entry[2];
  }
}

