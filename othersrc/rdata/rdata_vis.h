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

#ifndef _RDATA_VIS_H_
#define _RDATA_VIS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* 
 * Stuructures
 */
typedef struct rdata_vis
{
  char    filename[256];  // full path filename
  char    name[64];      // filename without full path
  int     num_pt,
          num_tri;
  float   *xyz;
  float   *nor;
  float   *rgb;
  int     *tri;
  double  M[16];
  int     flag_rgb;
  int     flag_display;
} rdata_vis;

int ReadRdataVis(char *name, rdata_vis *rd);


#endif
