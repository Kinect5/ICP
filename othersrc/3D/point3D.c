//
//   Copyright (C) 2009 ViGIR Lab (http://vigir.missouri.edu)
//   Written by Guilherme N. DeSouza <desouzag@missouri.edu>
//
//
/*
 * Some useful functions on 3D points
 * 
 * June 2001, Johnny Park
 */

#include "point3D.h"

/*
 * Eucleadian distance between two point_xyz's
 */
double E_distance(point_xyz *a, point_xyz *b)
{ 
  return(sqrt(SQ(a->x - b->x) + SQ(a->y - b->y) + SQ(a->z - b->z)));
}

/*
 * Angle between two unit normals
 */
double Angle(point_normal *a, point_normal *b)
{
  return( (180/PI) * acos((a->nx)*(b->nx)+(a->ny)*(b->ny)+(a->nz)*(b->nz)) );
}

/*
 * Finds center of mass of point_xyz's
 */
point_xyz CenterOfMass(point_xyz *a, int num_point)
{
  register int i;
  point_xyz center;

  center.x = 0; center.y = 0; center.z = 0;

  for (i=0; i<num_point; i++) {
    center.x += a[i].x;
    center.y += a[i].y;
    center.z += a[i].z;
  }
  center.x /= (double)num_point;
  center.y /= (double)num_point;
  center.z /= (double)num_point;

  return (center);
}

