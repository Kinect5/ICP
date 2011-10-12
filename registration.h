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

#ifndef _REGISTRATION_H_
#define _REGISTRATION_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "glut.h"
#include "matrix.h"
#include "point3D.h"
#include "rdata_vis.h"
#include "ICP.h"

void CourseRegistration(void);
void FineRegistration(void);
void SaveRegistration(void);

#endif

