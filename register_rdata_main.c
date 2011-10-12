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
 * Interactive registration 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <glut.h>
#include "trackball.h"
#include "rdata_vis.h"
#include "pick_corres_pt.h"
#include "registration.h"

#define NUM_COLOR 18
#define MAX_NUM_CORRES  30

/*
 * Function Declaration
 */
void DefineInitialView(void);
void Initialize(void);
void InitializeViewAnc(void);
void InitializeViewMov(void);
void InitializeViewRes(void);
void RecalculateModelViewAnc(void);
void RecalculateModelViewMov(void);
void RecalculateModelViewRes(void);
void ReadAllData(char **name);
void DisplayMainAnc(void);
void DisplayMainMov(void);
void DisplayMainRes(void);
void DisplayMenuAnc(void);
void DisplayMenuMov(void);
void DrawMenuAnc(GLenum mode);
void DrawMenuMov(GLenum mode);
void ReshapeMainAnc(int width, int height);
void ReshapeMainMov(int width, int height);
void ReshapeMainRes(int width, int height);
void ReshapeMenuAnc(int width, int height);
void ReshapeMenuMov(int width, int height);
void MouseMainAnc(int button, int state, int x, int y);
void MouseMainMov(int button, int state, int x, int y);
void MouseMainRes(int button, int state, int x, int y);
void MouseMenuAnc(int button, int state, int x, int y);
void MouseMenuMov(int button, int state, int x, int y);
int  ProcessHits(int hits, GLuint buffer[]);
void MouseMotionMainAnc(int x, int y);
void MouseMotionMainMov(int x, int y);
void MouseMotionMainRes(int x, int y);
void VisMainAnc(int visible);
void VisMainMov(int visible);
void VisMainRes(int visible);
void SpinAnc(void);
void SpinMov(void);
void MakeMenu(void);
void RegistrationMenu(int value);
void ViewMenu(int value);
void ShadeModelMenu(int value);
void ShininessMenu(int value);
void MainMenu(int value);

/*
 * Global Variables
 */
rdata_vis *rd_anc;            // array of "anchor" range data 
rdata_vis *rd_mov;            // array of "moving" range data

int     num_rdata_anc,
        num_rdata_mov,          
        display_option_anc = 1,
        display_option_mov = 1,
        draw_axes_anc = 1,
        draw_axes_mov = 1,
        draw_bbox_anc = 1,
        draw_bbox_mov = 1,
        lighting_anc = 1,
        lighting_mov = 1,
        common_color_anc = 0,
        common_color_mov = 0,
        spinning_anc = 0,
        spinning_mov = 0,
        spinning_res = 0,
        rotating_anc = 0,
        rotating_mov = 0,
        rotating_res = 0,
        scaling_anc = 0,
        scaling_mov = 0,
        scaling_res = 0,
        panning_anc = 0,
        panning_mov = 0,
        panning_res = 0;
        
float   curquat_anc[4],
        curquat_mov[4],
        curquat_res[4],
        lastquat_anc[4],
        lastquat_mov[4],
        lastquat_res[4],
        xpan_anc, ypan_anc, dist_anc,
        xpan_mov, ypan_mov, dist_mov,
        xpan_res, ypan_res, dist_res,
        xc_anc, yc_anc, zc_anc, r_anc,
        xc_mov, yc_mov, zc_mov, r_mov,
        xmin_anc, ymin_anc, zmin_anc, 
        xmax_anc, ymax_anc, zmax_anc,
        xmin_mov, ymin_mov, zmin_mov, 
        xmax_mov, ymax_mov, zmax_mov;
int     W_main_anc = 500, H_main_anc = 500,
        W_main_mov = 500, H_main_mov = 500,
        W_main_res = 500, H_main_res = 500,
        Wmenu_anc  , Hmenu_anc,
        Wmenu_mov  , Hmenu_mov,
        beginx, beginy;
int     new_model_anc, new_model_mov, new_model_res;
int     window_main_anc, window_menu_anc,
        window_main_mov, window_menu_mov,
        window_main_res;

float   ambient[4]={0.0, 0.0, 0.0, 1.0};         // ambient 
float   diffuse[4]={0.0, 0.0, 0.0, 1.0};         // diffuse
float   specular[4]={0.0, 0.0, 0.0, 1.0};        // specular
float   shininess_anc,
        shininess_mov;      // material shininess (0.0 - 128.0)

int     *corres_rd_anc,
        *corres_rd_mov,
        *corres_pt_anc,
        *corres_pt_mov;
int     num_corres_anc = 0,
        num_corres_mov = 0;
double  new_M[] = {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};



int main(int argc, char **argv)
{
  int     i;  

  // check argument
  if (argc != 3) {
    fprintf(stderr, "Argument Error!\n");
    exit(1);
  }

  ReadAllData(argv);

  // initialize flag_display
  for (i=0; i<num_rdata_anc; i++) rd_anc[i].flag_display = 1;
  for (i=0; i<num_rdata_mov; i++) rd_mov[i].flag_display = 1;

  // set menu size
  Wmenu_anc = Wmenu_mov = 400;
  Hmenu_anc = (num_rdata_anc+1)*25;
  Hmenu_mov = (num_rdata_mov+1)*25;

  // allocate memory
  corres_rd_anc = (int *) malloc (MAX_NUM_CORRES * sizeof(int));
  corres_rd_mov = (int *) malloc (MAX_NUM_CORRES * sizeof(int));
  corres_pt_anc = (int *) malloc (MAX_NUM_CORRES * sizeof(int));
  corres_pt_mov = (int *) malloc (MAX_NUM_CORRES * sizeof(int));

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
  
  // main window for anchor range data
  glutInitWindowPosition(50, 50);
  glutInitWindowSize(W_main_anc, H_main_anc);
  window_main_anc = glutCreateWindow("anchor range data");
  glutDisplayFunc(DisplayMainAnc);
  glutReshapeFunc(ReshapeMainAnc);
  glutVisibilityFunc(VisMainAnc);
  glutMouseFunc(MouseMainAnc);
  glutMotionFunc(MouseMotionMainAnc);
  MakeMenu();

  // menu window for anchor range data
  glutInitWindowPosition(50, 80+H_main_anc);
  glutInitWindowSize(Wmenu_anc, Hmenu_anc);
  window_menu_anc = glutCreateWindow("anchor range data menu");
  glutDisplayFunc(DisplayMenuAnc);
  glutReshapeFunc(ReshapeMenuAnc);
	glutMouseFunc(MouseMenuAnc);
  
  // main window for moving range data
  glutInitWindowPosition(60+W_main_anc, 50);
  glutInitWindowSize(W_main_mov, H_main_mov);
  window_main_mov = glutCreateWindow("moving range data");
  glutDisplayFunc(DisplayMainMov);
  glutReshapeFunc(ReshapeMainMov);
  glutVisibilityFunc(VisMainMov);
  glutMouseFunc(MouseMainMov);
  glutMotionFunc(MouseMotionMainMov);
  MakeMenu();
 
  // menu window for moving range data
  glutInitWindowPosition(60+W_main_anc, 80+H_main_anc);
  glutInitWindowSize(Wmenu_mov, Hmenu_mov);
  window_menu_anc = glutCreateWindow("moving range data menu");
  glutDisplayFunc(DisplayMenuMov);
  glutReshapeFunc(ReshapeMenuMov);
	glutMouseFunc(MouseMenuMov);

  // main window for registration result
  glutInitWindowPosition(70+W_main_anc+W_main_mov, 50);
  glutInitWindowSize(W_main_res, H_main_res);
  window_main_res = glutCreateWindow("registration result");
  glutDisplayFunc(DisplayMainRes);
  glutReshapeFunc(ReshapeMainRes);
  glutVisibilityFunc(VisMainRes);
  glutMouseFunc(MouseMainRes);
  glutMotionFunc(MouseMotionMainRes);
  
  DefineInitialView();
  Initialize();

  glutMainLoop();
  return 0;

}

void DefineInitialView(void)
{
  register int i, j;
  float x, y, z;

  // find bounding box
  xmin_anc = ymin_anc = zmin_anc = 999999;
  xmax_anc = ymax_anc = zmax_anc = -999999;
  xmin_mov = ymin_mov = zmin_mov = 999999;
  xmax_mov = ymax_mov = zmax_mov = -999999;
  for (j=0; j<num_rdata_anc; j++) {
    for (i=0; i<rd_anc[j].num_pt; i++) {
      // transform
      x = rd_anc[j].M[0] * rd_anc[j].xyz[3*i] +
          rd_anc[j].M[4] * rd_anc[j].xyz[3*i+1] +
          rd_anc[j].M[8] * rd_anc[j].xyz[3*i+2] + rd_anc[j].M[12];
      y = rd_anc[j].M[1] * rd_anc[j].xyz[3*i] +
          rd_anc[j].M[5] * rd_anc[j].xyz[3*i+1] +
          rd_anc[j].M[9] * rd_anc[j].xyz[3*i+2] + rd_anc[j].M[13];
      z = rd_anc[j].M[2] * rd_anc[j].xyz[3*i] +
          rd_anc[j].M[6] * rd_anc[j].xyz[3*i+1] +
          rd_anc[j].M[10]* rd_anc[j].xyz[3*i+2] + rd_anc[j].M[14];
      if (x < xmin_anc) xmin_anc = x;
      if (y < ymin_anc) ymin_anc = y;
      if (z < zmin_anc) zmin_anc = z;
    }
    for (i=0; i<rd_anc[j].num_pt; i++) {
      x = rd_anc[j].M[0] * rd_anc[j].xyz[3*i] +
          rd_anc[j].M[4] * rd_anc[j].xyz[3*i+1] +
          rd_anc[j].M[8] * rd_anc[j].xyz[3*i+2] + rd_anc[j].M[12];
      y = rd_anc[j].M[1] * rd_anc[j].xyz[3*i] +
          rd_anc[j].M[5] * rd_anc[j].xyz[3*i+1] +
          rd_anc[j].M[9] * rd_anc[j].xyz[3*i+2] + rd_anc[j].M[13];
      z = rd_anc[j].M[2] * rd_anc[j].xyz[3*i] +
          rd_anc[j].M[6] * rd_anc[j].xyz[3*i+1] +
          rd_anc[j].M[10]* rd_anc[j].xyz[3*i+2] + rd_anc[j].M[14];
      if (x > xmax_anc) xmax_anc = x;
      if (y > ymax_anc) ymax_anc = y;
      if (z > zmax_anc) zmax_anc = z;
    }
  }
  for (j=0; j<num_rdata_mov; j++) {
    for (i=0; i<rd_mov[j].num_pt; i++) {
      // transform
      x = rd_mov[j].M[0] * rd_mov[j].xyz[3*i] +
          rd_mov[j].M[4] * rd_mov[j].xyz[3*i+1] +
          rd_mov[j].M[8] * rd_mov[j].xyz[3*i+2] + rd_mov[j].M[12];
      y = rd_mov[j].M[1] * rd_mov[j].xyz[3*i] +
          rd_mov[j].M[5] * rd_mov[j].xyz[3*i+1] +
          rd_mov[j].M[9] * rd_mov[j].xyz[3*i+2] + rd_mov[j].M[13];
      z = rd_mov[j].M[2] * rd_mov[j].xyz[3*i] +
          rd_mov[j].M[6] * rd_mov[j].xyz[3*i+1] +
          rd_mov[j].M[10]* rd_mov[j].xyz[3*i+2] + rd_mov[j].M[14];
      if (x < xmin_mov) xmin_mov = x;
      if (y < ymin_mov) ymin_mov = y;
      if (z < zmin_mov) zmin_mov = z;
    }
    for (i=0; i<rd_mov[j].num_pt; i++) {
      x = rd_mov[j].M[0] * rd_mov[j].xyz[3*i] +
          rd_mov[j].M[4] * rd_mov[j].xyz[3*i+1] +
          rd_mov[j].M[8] * rd_mov[j].xyz[3*i+2] + rd_mov[j].M[12];
      y = rd_mov[j].M[1] * rd_mov[j].xyz[3*i] +
          rd_mov[j].M[5] * rd_mov[j].xyz[3*i+1] +
          rd_mov[j].M[9] * rd_mov[j].xyz[3*i+2] + rd_mov[j].M[13];
      z = rd_mov[j].M[2] * rd_mov[j].xyz[3*i] +
          rd_mov[j].M[6] * rd_mov[j].xyz[3*i+1] +
          rd_mov[j].M[10]* rd_mov[j].xyz[3*i+2] + rd_mov[j].M[14];
      if (x > xmax_mov) xmax_mov = x;
      if (y > ymax_mov) ymax_mov = y;
      if (z > zmax_mov) zmax_mov = z;
    }
  }

  // compute center of bounding box
  xc_anc = 0.5*(xmin_anc + xmax_anc);
  yc_anc = 0.5*(ymin_anc + ymax_anc);
  zc_anc = 0.5*(zmin_anc + zmax_anc);
  xc_mov = 0.5*(xmin_mov + xmax_mov);
  yc_mov = 0.5*(ymin_mov + ymax_mov);
  zc_mov = 0.5*(zmin_mov + zmax_mov);

  // radius of a sphere that contains the bounding box
  r_anc = sqrt( (xmax_anc - xc_anc) * (xmax_anc - xc_anc) + 
                (ymax_anc - yc_anc) * (ymax_anc - yc_anc) + 
                (zmax_anc - zc_anc) * (zmax_anc - zc_anc)  );
  r_mov = sqrt( (xmax_mov - xc_mov) * (xmax_mov - xc_mov) + 
                (ymax_mov - yc_mov) * (ymax_mov - yc_mov) + 
                (zmax_mov - zc_mov) * (zmax_mov - zc_mov)  );
}


void Initialize(void) 
{				
  float light[4];
  float lmodel[4];

  // set default reflectance properties
  shininess_anc = shininess_mov = 20.0;
  
  glutSetWindow(window_main_anc);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  glShadeModel(GL_SMOOTH);
  glutSetWindow(window_main_mov);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  glShadeModel(GL_SMOOTH);
  glutSetWindow(window_main_res);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  glShadeModel(GL_SMOOTH);
  // light properties:
  light[0] = 1.0; light[1] = 1.0; light[2] = 1.0; light[3] = 1.0;
  glLightfv(GL_LIGHT0, GL_AMBIENT, light);
  light[0] = 1.0; light[1] = 1.0; light[2] = 1.0; light[3] = 1.0;
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light);
  light[0] = 1.0; light[1] = 1.0; light[2] = 1.0; light[3] = 1.0;
  glLightfv(GL_LIGHT0, GL_SPECULAR, light);
  // light is always at the origin of the eye coordinate
  light[0] = 0.0; light[1] = 0.0; light[2] = 0.0; light[3] = 1.0;
  glLightfv(GL_LIGHT0, GL_POSITION, light);
  // ambient RGBA intensity of the entire scene
  lmodel[0] = 0.1; lmodel[1] = 0.1; lmodel[2] = 0.1;
  glutSetWindow(window_main_anc);
  glEnable(GL_LIGHT0);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel);
  glutSetWindow(window_main_mov);
  glEnable(GL_LIGHT0);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel);
  glutSetWindow(window_main_res);
  glEnable(GL_LIGHT0);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel);

  // set projection matrix
  glutSetWindow(window_main_anc);
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPerspective( 45.0, 1.0, r_anc/10.0, r_anc*5);
  glutSetWindow(window_main_mov);
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPerspective( 45.0, 1.0, r_mov/10.0, r_mov*5);
  glutSetWindow(window_main_res);
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPerspective( 45.0, 1.0, r_anc/10.0, r_anc*5);

  // initialize view
  InitializeViewAnc();
  InitializeViewMov();
  InitializeViewRes();
}


void InitializeViewAnc(void)
{
  trackball(curquat_anc, 0.0, 0.0, 0.0, 0.0);
  xpan_anc = ypan_anc = dist_anc = 0;
  RecalculateModelViewAnc();
}

void InitializeViewMov(void)
{
  trackball(curquat_mov, 0.0, 0.0, 0.0, 0.0);
  xpan_mov = ypan_mov = dist_mov = 0;
  RecalculateModelViewMov();
}

void InitializeViewRes(void)
{
  trackball(curquat_res, 0.0, 0.0, 0.0, 0.0);
  xpan_res = ypan_res = dist_res = 0;
  RecalculateModelViewRes();
}

void RecalculateModelViewAnc(void)
{
  GLfloat m[4][4];

  glutSetWindow(window_main_anc);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(xpan_anc, ypan_anc, -dist_anc);
  glTranslatef(xc_anc, -zc_anc, -3.0*r_anc);
  build_rotmatrix(m, curquat_anc);
  glMultMatrixf(&m[0][0]);
  glRotatef(180, 0.0, 1.0, 0.0);
  glRotatef(-90, 1.0, 0.0, 0.0);

  new_model_anc = 0;
}

void RecalculateModelViewMov(void)
{
  GLfloat m[4][4];

  glutSetWindow(window_main_mov);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(xpan_mov, ypan_mov, -dist_mov);
  glTranslatef(xc_mov, -zc_mov, -3.0*r_mov);
  build_rotmatrix(m, curquat_mov);
  glMultMatrixf(&m[0][0]);
  glRotatef(180, 0.0, 1.0, 0.0);
  glRotatef(-90, 1.0, 0.0, 0.0);

  new_model_mov = 0;
}

void RecalculateModelViewRes(void)
{
  GLfloat m[4][4];

  glutSetWindow(window_main_res);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(xpan_res, ypan_res, -dist_res);
  glTranslatef(xc_anc, -zc_anc, -3.0*r_anc);
  build_rotmatrix(m, curquat_res);
  glMultMatrixf(&m[0][0]);
  glRotatef(180, 0.0, 1.0, 0.0);
  glRotatef(-90, 1.0, 0.0, 0.0);

  new_model_res = 0;
}

void ReadAllData(char **name)
{
  int     i;
  FILE    *fp;
  char    filename[256];

 // anchor data
  if ( (fp = fopen(name[1], "r")) != NULL) { // multiple range data
    fscanf(fp, "%d\n", &num_rdata_anc);
    rd_anc = (rdata_vis *) malloc (num_rdata_anc * sizeof(rdata_vis));
    for (i=0; i<num_rdata_anc; i++) {
      fscanf(fp, "%s\n", filename);
      if (!ReadRdataVis(filename, &rd_anc[i])) {
        fprintf(stderr, "Exiting...\n");
        exit(1);
      }
    }
  }
  else { // single range data
    num_rdata_anc = 1;
    rd_anc = (rdata_vis *) malloc (sizeof(rdata_vis));  // I have no idea why, but there seems to be
    rd_anc = (rdata_vis *) malloc (sizeof(rdata_vis));  // some bug with malloc.  So, call twice here!!
    if (!ReadRdataVis(name[1], &rd_anc[0])) {
      fprintf(stderr, "Exiting...\n");
      exit(1);
    }
  }  

 // moving data
  if ( (fp = fopen(name[2], "r")) != NULL) { // multiple range data
    fscanf(fp, "%d\n", &num_rdata_mov);
    rd_mov = (rdata_vis *) malloc (num_rdata_mov * sizeof(rdata_vis));
    for (i=0; i<num_rdata_mov; i++) {
      fscanf(fp, "%s\n", filename);
      if (!ReadRdataVis(filename, &rd_mov[i])) {
        fprintf(stderr, "Exiting...\n");
        exit(1);
      }
    }
  }
  else { // single range data
    num_rdata_mov = 1;
    rd_mov = (rdata_vis *) malloc (sizeof(rdata_vis));
    if (!ReadRdataVis(name[2], &rd_mov[0])) {
      fprintf(stderr, "Exiting...\n");
      exit(1);
    }
  }  
}


void DisplayMainAnc(void)
{
  int     i, j;  
  float   v[3];
  char    num[2];
  double  MV[16];
  double  MVt[16];

  if (new_model_anc) {
    RecalculateModelViewAnc();
  }

  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);

  if (lighting_anc) glEnable(GL_LIGHTING);
  else  glDisable(GL_LIGHTING);

  // assign reflectance
  diffuse[0] = 0.2; diffuse[1] = 0.2; diffuse[2] = 0.8;
  for (i=0; i<3; i++) ambient[i] = 0.3 * diffuse[i];
  glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
  glMaterialf(GL_FRONT, GL_SHININESS, shininess_anc);	
  glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
  glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
  glColor4fv(diffuse);

  for (j=0; j<num_rdata_anc; j++) {
    if (!rd_anc[j].flag_display) continue;

    // apply M matrix
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrixd(rd_anc[j].M);

    // draw data with triangles mesh
    if (display_option_anc == 1) {
      // reference array elements
      glEnableClientState(GL_VERTEX_ARRAY);
      glEnableClientState(GL_NORMAL_ARRAY);
      glDisableClientState(GL_COLOR_ARRAY);	
	
      glVertexPointer(3, GL_FLOAT, 0, rd_anc[j].xyz);
      glNormalPointer(GL_FLOAT, 0, rd_anc[j].nor);

      glPolygonMode(GL_FRONT, GL_FILL);
      glDrawElements(GL_TRIANGLES, rd_anc[j].num_tri*3, 
                     GL_UNSIGNED_INT, rd_anc[j].tri);
    }

    // draw data with triangle lines
    else if (display_option_anc == 2) { 
      // reference array elements
      glEnableClientState(GL_VERTEX_ARRAY);
      glEnableClientState(GL_NORMAL_ARRAY);
      glDisableClientState(GL_COLOR_ARRAY);
	
      glVertexPointer(3, GL_FLOAT, 0, rd_anc[j].xyz);
      glNormalPointer(GL_FLOAT, 0, rd_anc[j].nor);

      glPolygonMode(GL_FRONT, GL_LINE);
      glDrawElements(GL_TRIANGLES, rd_anc[j].num_tri*3, 
                     GL_UNSIGNED_INT, rd_anc[j].tri);
    }

    // draw data with points
    else if (display_option_anc == 3) { 	
      // reference array elements
      glEnableClientState(GL_VERTEX_ARRAY);
      glEnableClientState(GL_NORMAL_ARRAY);
      glDisableClientState(GL_COLOR_ARRAY);
	
      glVertexPointer(3, GL_FLOAT, 0, rd_anc[j].xyz);
      glNormalPointer(GL_FLOAT, 0, rd_anc[j].nor);

      glDrawArrays(GL_POINTS, 0, rd_anc[j].num_pt*3);
    }
  
    // draw data with its normal as color
    else if (display_option_anc == 4) { 
      glDisable(GL_LIGHTING);
	
      // reference array elements
      glEnableClientState(GL_VERTEX_ARRAY);
      glDisableClientState(GL_NORMAL_ARRAY);
      glEnableClientState(GL_COLOR_ARRAY);
	
      glVertexPointer(3, GL_FLOAT, 0, rd_anc[j].xyz);
      glColorPointer(3, GL_FLOAT, 0, rd_anc[j].nor);

      glPolygonMode(GL_FRONT, GL_FILL);
      glDrawElements(GL_TRIANGLES, rd_anc[j].num_tri*3, 
                     GL_UNSIGNED_INT, rd_anc[j].tri);
    }

    else if (display_option_anc == 5) {
      if (rd_anc[j].flag_rgb) {
        glDisable(GL_LIGHTING);

        // reference array elements
        glEnableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
        glEnableClientState(GL_COLOR_ARRAY);
	
        glVertexPointer(3, GL_FLOAT, 0, rd_anc[j].xyz);
        glColorPointer(3, GL_FLOAT, 0, rd_anc[j].rgb);

        glPolygonMode(GL_FRONT, GL_FILL);
        glDrawElements(GL_TRIANGLES, rd_anc[j].num_tri*3, 
                       GL_UNSIGNED_INT, rd_anc[j].tri);
      }
      // if rgb is not available, draw with triable mesh
      else {
        // reference array elements
        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_COLOR_ARRAY);	
	
        glVertexPointer(3, GL_FLOAT, 0, rd_anc[j].xyz);
        glNormalPointer(GL_FLOAT, 0, rd_anc[j].nor);

        glPolygonMode(GL_FRONT, GL_FILL);
        glDrawElements(GL_TRIANGLES, rd_anc[j].num_tri*3, 
                       GL_UNSIGNED_INT, rd_anc[j].tri);
      }
    }
    
    glPopMatrix();
  }

  // draw axes
  if (draw_axes_anc) {
    glDisable(GL_LIGHTING);
    glBegin( GL_LINES );
    // x axis
    glColor3d( 1, 0, 0 );
    glVertex3f(  0.0, 0.0, 0.0 );
    glVertex4f(  1.0, 0.0, 0.0, 0.0 );

    // y axis
    glColor3d( 0, 1, 0 );
    glVertex3f( 0.0, 0.0, 0.0 );
    glVertex4f( 0.0, 1.0, 0.0, 0.0 );

    // z axis
    glColor3d( 0, 0, 1);
    glVertex3f( 0.0, 0.0, 0.0 );
    glVertex4f( 0.0, 0.0, 1.0, 0.0 );
    glEnd();
  }

  // draw bounding box
  if (draw_bbox_anc) {
    glDisable(GL_LIGHTING);
    glColor3d(1, 1, 0);

    glLineStipple(4, 0xAAAA);
    glEnable(GL_LINE_STIPPLE);

    glBegin(GL_LINE_LOOP);
    glVertex3f(xmin_anc, ymax_anc, zmax_anc);
    glVertex3f(xmax_anc, ymax_anc, zmax_anc);
    glVertex3f(xmax_anc, ymax_anc, zmin_anc);
    glVertex3f(xmin_anc, ymax_anc, zmin_anc);
    glEnd();
  
    glBegin(GL_LINE_LOOP);
    glVertex3f(xmin_anc, ymin_anc, zmax_anc);
    glVertex3f(xmax_anc, ymin_anc, zmax_anc);
    glVertex3f(xmax_anc, ymin_anc, zmin_anc);
    glVertex3f(xmin_anc, ymin_anc, zmin_anc);
    glEnd();

    glBegin(GL_LINES);
    glVertex3f(xmin_anc, ymin_anc, zmax_anc);
    glVertex3f(xmin_anc, ymax_anc, zmax_anc);
    glVertex3f(xmax_anc, ymin_anc, zmax_anc);
    glVertex3f(xmax_anc, ymax_anc, zmax_anc);
    glVertex3f(xmax_anc, ymin_anc, zmin_anc);
    glVertex3f(xmax_anc, ymax_anc, zmin_anc);
    glVertex3f(xmin_anc, ymin_anc, zmin_anc);
    glVertex3f(xmin_anc, ymax_anc, zmin_anc);
    glEnd();
    
    glDisable(GL_LINE_STIPPLE);
  }

  // draw corresponding points
  glDisable(GL_DEPTH_TEST); 
  glDisable(GL_LIGHTING);
  glGetDoublev(GL_MODELVIEW_MATRIX, MV);
  for (i=0; i<16; i++) MVt[i] = MV[i];
  MVt[1] = MV[4]; MVt[4] = MV[1];
  MVt[2] = MV[8]; MVt[8] = MV[2];
  MVt[6] = MV[9]; MVt[9] = MV[6];
  MVt[12] = MVt[13] = MVt[14] = 0;
  glColor3f(1.0, 1.0, 1.0);
  for (i=0; i<num_corres_anc; i++) {
    v[0] = (rd_anc[corres_rd_anc[i]].xyz[3*corres_pt_anc[i]] *
            rd_anc[corres_rd_anc[i]].M[0]) +
           (rd_anc[corres_rd_anc[i]].xyz[3*corres_pt_anc[i]+1] *
            rd_anc[corres_rd_anc[i]].M[4]) +
           (rd_anc[corres_rd_anc[i]].xyz[3*corres_pt_anc[i]+2] *
            rd_anc[corres_rd_anc[i]].M[8]) +
            rd_anc[corres_rd_anc[i]].M[12];
    v[1] = (rd_anc[corres_rd_anc[i]].xyz[3*corres_pt_anc[i]] *
            rd_anc[corres_rd_anc[i]].M[1]) +
           (rd_anc[corres_rd_anc[i]].xyz[3*corres_pt_anc[i]+1] *
            rd_anc[corres_rd_anc[i]].M[5]) +
           (rd_anc[corres_rd_anc[i]].xyz[3*corres_pt_anc[i]+2] *
            rd_anc[corres_rd_anc[i]].M[9]) +
            rd_anc[corres_rd_anc[i]].M[13];
    v[2] = (rd_anc[corres_rd_anc[i]].xyz[3*corres_pt_anc[i]] *
            rd_anc[corres_rd_anc[i]].M[2]) +
           (rd_anc[corres_rd_anc[i]].xyz[3*corres_pt_anc[i]+1] *
            rd_anc[corres_rd_anc[i]].M[6]) +
           (rd_anc[corres_rd_anc[i]].xyz[3*corres_pt_anc[i]+2] *
            rd_anc[corres_rd_anc[i]].M[10]) +
            rd_anc[corres_rd_anc[i]].M[14];
    glPointSize(6.0);
    glBegin(GL_POINTS);
    glVertex4f(v[0], v[1], v[2], 1.0);
    glEnd();

    glPushMatrix();
    glLoadIdentity();
    glMultMatrixd(MV);
    glTranslatef(v[0], v[1], v[2]);
    glMultMatrixd(MVt);
    glScalef(0.05, 0.05, 0.05);
    
    sprintf(num, "%d", i);
    glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, num[0]);
    if (i >= 10)
      glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, num[1]);
    
    glPopMatrix();
  }

  glutSwapBuffers();
}


void DisplayMainMov(void)
{
  int     i, j;  
  float   v[3];
  char    num[2];
  double  MV[16], MVt[16];

  if (new_model_mov) {
    RecalculateModelViewMov();
  }

  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);

  if (lighting_mov) glEnable(GL_LIGHTING);
  else  glDisable(GL_LIGHTING);

  // assign reflectance
  diffuse[0] = 0.8; diffuse[1] = 0.2; diffuse[2] = 0.2;
  for (i=0; i<3; i++) ambient[i] = 0.3 * diffuse[i];
  glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
  glMaterialf(GL_FRONT, GL_SHININESS, shininess_mov );	
  glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
  glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
  glColor4fv(diffuse);
 
  for (j=0; j<num_rdata_mov; j++) {
    if (!rd_mov[j].flag_display) continue;

    // apply M matrix
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrixd(rd_mov[j].M);

    // draw data with triangles mesh
    if (display_option_mov == 1) {
      // reference array elements
      glEnableClientState(GL_VERTEX_ARRAY);
      glEnableClientState(GL_NORMAL_ARRAY);
      glDisableClientState(GL_COLOR_ARRAY);	
	
      glVertexPointer(3, GL_FLOAT, 0, rd_mov[j].xyz);
      glNormalPointer(GL_FLOAT, 0, rd_mov[j].nor);

      glPolygonMode(GL_FRONT, GL_FILL);
      glDrawElements(GL_TRIANGLES, rd_mov[j].num_tri*3, 
                     GL_UNSIGNED_INT, rd_mov[j].tri);
    }

    // draw data with triangle lines
    else if (display_option_mov == 2) { 
      // reference array elements
      glEnableClientState(GL_VERTEX_ARRAY);
      glEnableClientState(GL_NORMAL_ARRAY);
      glDisableClientState(GL_COLOR_ARRAY);
	
      glVertexPointer(3, GL_FLOAT, 0, rd_mov[j].xyz);
      glNormalPointer(GL_FLOAT, 0, rd_mov[j].nor);

      glPolygonMode(GL_FRONT, GL_LINE);
      glDrawElements(GL_TRIANGLES, rd_mov[j].num_tri*3, 
                     GL_UNSIGNED_INT, rd_mov[j].tri);
    }

    // draw data with points
    else if (display_option_mov == 3) { 	
      // reference array elements
      glEnableClientState(GL_VERTEX_ARRAY);
      glEnableClientState(GL_NORMAL_ARRAY);
      glDisableClientState(GL_COLOR_ARRAY);
	
      glVertexPointer(3, GL_FLOAT, 0, rd_mov[j].xyz);
      glNormalPointer(GL_FLOAT, 0, rd_mov[j].nor);

      glDrawArrays(GL_POINTS, 0, rd_mov[j].num_pt*3);
    }
  
    // draw data with its normal as color
    else if (display_option_mov == 4) { 
      glDisable(GL_LIGHTING);
	
      // reference array elements
      glEnableClientState(GL_VERTEX_ARRAY);
      glDisableClientState(GL_NORMAL_ARRAY);
      glEnableClientState(GL_COLOR_ARRAY);
	
      glVertexPointer(3, GL_FLOAT, 0, rd_mov[j].xyz);
      glColorPointer(3, GL_FLOAT, 0, rd_mov[j].nor);

      glPolygonMode(GL_FRONT, GL_FILL);
      glDrawElements(GL_TRIANGLES, rd_mov[j].num_tri*3, 
                     GL_UNSIGNED_INT, rd_mov[j].tri);
    }

    else if (display_option_mov == 5) {
      if (rd_mov[j].flag_rgb) {
        glDisable(GL_LIGHTING);

        // reference array elements
        glEnableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
        glEnableClientState(GL_COLOR_ARRAY);
	
        glVertexPointer(3, GL_FLOAT, 0, rd_mov[j].xyz);
        glColorPointer(3, GL_FLOAT, 0, rd_mov[j].rgb);

        glPolygonMode(GL_FRONT, GL_FILL);
        glDrawElements(GL_TRIANGLES, rd_mov[j].num_tri*3, 
                       GL_UNSIGNED_INT, rd_mov[j].tri);
      }
      // if rgb is not available, draw with triable mesh
      else {
        // reference array elements
        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_COLOR_ARRAY);	
	
        glVertexPointer(3, GL_FLOAT, 0, rd_mov[j].xyz);
        glNormalPointer(GL_FLOAT, 0, rd_mov[j].nor);

        glPolygonMode(GL_FRONT, GL_FILL);
        glDrawElements(GL_TRIANGLES, rd_mov[j].num_tri*3, 
                       GL_UNSIGNED_INT, rd_mov[j].tri);
      }
    }
    glPopMatrix();
  }

  // draw axes
  if (draw_axes_mov) {
    glDisable(GL_LIGHTING);
    glBegin( GL_LINES );
    // x axis
    glColor3d( 1, 0, 0 );
    glVertex3f(  0.0, 0.0, 0.0 );
    glVertex4f(  1.0, 0.0, 0.0, 0.0 );

    // y axis
    glColor3d( 0, 1, 0 );
    glVertex3f( 0.0, 0.0, 0.0 );
    glVertex4f( 0.0, 1.0, 0.0, 0.0 );

    // z axis
    glColor3d( 0, 0, 1);
    glVertex3f( 0.0, 0.0, 0.0 );
    glVertex4f( 0.0, 0.0, 1.0, 0.0 );
    glEnd();
  }

  // draw bounding box
  if (draw_bbox_mov) {
    glDisable(GL_LIGHTING);
    glColor3d(1, 1, 0);

    glLineStipple(4, 0xAAAA);
    glEnable(GL_LINE_STIPPLE);

    glBegin(GL_LINE_LOOP);
    glVertex3f(xmin_mov, ymax_mov, zmax_mov);
    glVertex3f(xmax_mov, ymax_mov, zmax_mov);
    glVertex3f(xmax_mov, ymax_mov, zmin_mov);
    glVertex3f(xmin_mov, ymax_mov, zmin_mov);
    glEnd();
  
    glBegin(GL_LINE_LOOP);
    glVertex3f(xmin_mov, ymin_mov, zmax_mov);
    glVertex3f(xmax_mov, ymin_mov, zmax_mov);
    glVertex3f(xmax_mov, ymin_mov, zmin_mov);
    glVertex3f(xmin_mov, ymin_mov, zmin_mov);
    glEnd();

    glBegin(GL_LINES);
    glVertex3f(xmin_mov, ymin_mov, zmax_mov);
    glVertex3f(xmin_mov, ymax_mov, zmax_mov);
    glVertex3f(xmax_mov, ymin_mov, zmax_mov);
    glVertex3f(xmax_mov, ymax_mov, zmax_mov);
    glVertex3f(xmax_mov, ymin_mov, zmin_mov);
    glVertex3f(xmax_mov, ymax_mov, zmin_mov);
    glVertex3f(xmin_mov, ymin_mov, zmin_mov);
    glVertex3f(xmin_mov, ymax_mov, zmin_mov);
    glEnd();
    
    glDisable(GL_LINE_STIPPLE);
  }

  // draw corresponding points
  glDisable(GL_DEPTH_TEST); 
  glDisable(GL_LIGHTING);
  glGetDoublev(GL_MODELVIEW_MATRIX, MV);
  for (i=0; i<16; i++) MVt[i] = MV[i];
  MVt[1] = MV[4]; MVt[4] = MV[1];
  MVt[2] = MV[8]; MVt[8] = MV[2];
  MVt[6] = MV[9]; MVt[9] = MV[6];
  MVt[12] = MVt[13] = MVt[14] = 0;
  glColor3f(1.0, 1.0, 1.0);
  for (i=0; i<num_corres_mov; i++) {
    v[0] = (rd_mov[corres_rd_mov[i]].xyz[3*corres_pt_mov[i]] *
            rd_mov[corres_rd_mov[i]].M[0]) +
           (rd_mov[corres_rd_mov[i]].xyz[3*corres_pt_mov[i]+1] *
            rd_mov[corres_rd_mov[i]].M[4]) +
           (rd_mov[corres_rd_mov[i]].xyz[3*corres_pt_mov[i]+2] *
            rd_mov[corres_rd_mov[i]].M[8]) +
            rd_mov[corres_rd_mov[i]].M[12];
    v[1] = (rd_mov[corres_rd_mov[i]].xyz[3*corres_pt_mov[i]] *
            rd_mov[corres_rd_mov[i]].M[1]) +
           (rd_mov[corres_rd_mov[i]].xyz[3*corres_pt_mov[i]+1] *
            rd_mov[corres_rd_mov[i]].M[5]) +
           (rd_mov[corres_rd_mov[i]].xyz[3*corres_pt_mov[i]+2] *
            rd_mov[corres_rd_mov[i]].M[9]) +
            rd_mov[corres_rd_mov[i]].M[13];
    v[2] = (rd_mov[corres_rd_mov[i]].xyz[3*corres_pt_mov[i]] *
            rd_mov[corres_rd_mov[i]].M[2]) +
           (rd_mov[corres_rd_mov[i]].xyz[3*corres_pt_mov[i]+1] *
            rd_mov[corres_rd_mov[i]].M[6]) +
           (rd_mov[corres_rd_mov[i]].xyz[3*corres_pt_mov[i]+2] *
            rd_mov[corres_rd_mov[i]].M[10]) +
            rd_mov[corres_rd_mov[i]].M[14];
    glPointSize(6.0);
    glBegin(GL_POINTS);
    glVertex4f(v[0], v[1], v[2], 1.0);
    glEnd();

    glPushMatrix();
    glLoadIdentity();
    glMultMatrixd(MV);
    glTranslatef(v[0], v[1], v[2]);
    glMultMatrixd(MVt);
    glScalef(0.05, 0.05, 0.05);
    sprintf(num, "%d", i);
    glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, num[0]);
    if (i >= 10)
      glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, num[1]);
    
    glPopMatrix();
  }

  glutSwapBuffers();
}

void DisplayMainRes(void)
{
  int     i, j;

  if (new_model_res) {
    RecalculateModelViewRes();
  }

  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);

  glEnable(GL_LIGHTING);

  // assign reflectance
  diffuse[0] = 0.2; diffuse[1] = 0.2; diffuse[2] = 0.8; 
  for (i=0; i<3; i++) ambient[i] = 0.3 * diffuse[i];
  glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
  glMaterialf(GL_FRONT, GL_SHININESS, shininess_anc );	
  glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
  glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
  glColor4fv(diffuse);

  // draw anchor range data
  for (j=0; j<num_rdata_anc; j++) {

    // apply M matrix
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrixd(rd_anc[j].M);

    // reference array elements
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);	
	
    glVertexPointer(3, GL_FLOAT, 0, rd_anc[j].xyz);
    glNormalPointer(GL_FLOAT, 0, rd_anc[j].nor);

    glPolygonMode(GL_FRONT, GL_FILL);
    glDrawElements(GL_TRIANGLES, rd_anc[j].num_tri*3, 
                   GL_UNSIGNED_INT, rd_anc[j].tri);

    glPopMatrix();
  }

  // assign reflectance
  diffuse[0] = 0.8; diffuse[1] = 0.2; diffuse[2] = 0.2; 
  for (i=0; i<3; i++) ambient[i] = 0.3 * diffuse[i];
  glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
  glMaterialf(GL_FRONT, GL_SHININESS, shininess_mov );	
  glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
  glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
  glColor4fv(diffuse);

  // draw moving range data
  for (j=0; j<num_rdata_mov; j++) {
    // apply new_M and M matrix
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrixd(new_M);
    glMultMatrixd(rd_mov[j].M);

    // reference array elements
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);	
	
    glVertexPointer(3, GL_FLOAT, 0, rd_mov[j].xyz);
    glNormalPointer(GL_FLOAT, 0, rd_mov[j].nor);

    glPolygonMode(GL_FRONT, GL_FILL);
    glDrawElements(GL_TRIANGLES, rd_mov[j].num_tri*3, 
                   GL_UNSIGNED_INT, rd_mov[j].tri);

    glPopMatrix();
  }

  // draw axes
  glDisable(GL_LIGHTING);
  glBegin( GL_LINES );
  // x axis
  glColor3d( 1, 0, 0 );
  glVertex3f(  0.0, 0.0, 0.0 );
  glVertex4f(  1.0, 0.0, 0.0, 0.0 );

  // y axis
  glColor3d( 0, 1, 0 );
  glVertex3f( 0.0, 0.0, 0.0 );
  glVertex4f( 0.0, 1.0, 0.0, 0.0 );

  // z axis
  glColor3d( 0, 0, 1);
  glVertex3f( 0.0, 0.0, 0.0 );
  glVertex4f( 0.0, 0.0, 1.0, 0.0 );
  glEnd();
 
  glutSwapBuffers();
}


void DisplayMenuAnc(void)
{
  glClearColor(0.7, 0.7, 0.7, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);

  DrawMenuAnc(GL_RENDER);

  glutSwapBuffers();
}


void DisplayMenuMov(void)
{
  glClearColor(0.7, 0.7, 0.7, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);

  DrawMenuMov(GL_RENDER);

  glutSwapBuffers();
}


void DrawMenuAnc(GLenum mode)
{
  int   i, j;
  char  string[256];
  int   len;
  float v_offset, dv,
        h_offset;

  v_offset = 10.0 / (float)Hmenu_anc;
  h_offset = 10.0 / (float)Wmenu_anc;
  dv = 1.0 / (float)(num_rdata_anc+1);
 
  glColor3f(0.0, 0.0, 0.5);
  glRasterPos2f(0.4, 1.0-dv+v_offset);
  strcpy(string, "num_pt");
  for (i = 0; i < 6; i++) {
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, string[i]);
  }
  glRasterPos2f(0.6, 1.0-dv+v_offset);
  strcpy(string, "num_tri");
  for (i = 0; i < 7; i++) {
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, string[i]);
  }
  glRasterPos2f(0.8, 1.0-dv+v_offset);
  strcpy(string, "display");
  for (i = 0; i < 7; i++) {
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, string[i]);
  }

  for (i=0; i<num_rdata_anc; i++) {
    // print name
    glColor3f(0.5, 0.0, 0.0);
    glRasterPos2f(h_offset, 1.0-dv*(i+2)+v_offset);
    len = strlen(rd_anc[i].name);
    for (j=0; j<len; j++) {
      glutBitmapCharacter(GLUT_BITMAP_9_BY_15, rd_anc[i].name[j]);
    }
    // print number of points
    glColor3f(0.0, 0.0, 0.0);
    glRasterPos2f(0.4, 1.0-dv*(i+2)+v_offset);
    sprintf(string, "%d", rd_anc[i].num_pt);
    len = strlen(string);
    for (j=0; j<len; j++) {
      glutBitmapCharacter(GLUT_BITMAP_9_BY_15, string[j]);
    }
    // print number of triangles
    glRasterPos2f(0.6, 1.0-dv*(i+2)+v_offset);
    sprintf(string, "%d", rd_anc[i].num_tri);
    len = strlen(string);
    for (j=0; j<len; j++) {
      glutBitmapCharacter(GLUT_BITMAP_9_BY_15, string[j]);
    }
    // draw displaying status
    if (mode == GL_SELECT) glLoadName(i);
    glBegin(GL_QUADS);
    if (rd_anc[i].flag_display) {
      glColor3f(0.2, 0.2, 0.8);
    }
    else glColor3f(0, 0, 0);
    glVertex3f(0.86, 1.0-dv*(i+2)+v_offset, 0.0);
    glVertex3f(0.89, 1.0-dv*(i+2)+v_offset, 0.0);
    glVertex3f(0.89, 1.0-dv*(i+1.6)+v_offset, 0.0);
    glVertex3f(0.86, 1.0-dv*(i+1.6)+v_offset, 0.0);
    glEnd();
  }
}


void DrawMenuMov(GLenum mode)
{
  int   i, j;
  char  string[256];
  int   len;
  float v_offset, dv,
        h_offset;

  v_offset = 10.0 / (float)Hmenu_mov;
  h_offset = 10.0 / (float)Wmenu_mov;
  dv = 1.0 / (float)(num_rdata_mov+1);
 
  glColor3f(0.0, 0.0, 0.5);
  glRasterPos2f(0.4, 1.0-dv+v_offset);
  strcpy(string, "num_pt");
  for (i = 0; i < 6; i++) {
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, string[i]);
  }
  glRasterPos2f(0.6, 1.0-dv+v_offset);
  strcpy(string, "num_tri");
  for (i = 0; i < 7; i++) {
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, string[i]);
  }
  glRasterPos2f(0.8, 1.0-dv+v_offset);
  strcpy(string, "display");
  for (i = 0; i < 7; i++) {
    glutBitmapCharacter(GLUT_BITMAP_9_BY_15, string[i]);
  }

  for (i=0; i<num_rdata_mov; i++) {
    // print name
    glColor3f(0.5, 0.0, 0.0);
    glRasterPos2f(h_offset, 1.0-dv*(i+2)+v_offset);
    len = strlen(rd_mov[i].name);
    for (j=0; j<len; j++) {
      glutBitmapCharacter(GLUT_BITMAP_9_BY_15, rd_mov[i].name[j]);
    }
    // print number of points
    glColor3f(0.0, 0.0, 0.0);
    glRasterPos2f(0.4, 1.0-dv*(i+2)+v_offset);
    sprintf(string, "%d", rd_mov[i].num_pt);
    len = strlen(string);
    for (j=0; j<len; j++) {
      glutBitmapCharacter(GLUT_BITMAP_9_BY_15, string[j]);
    }
    // print number of triangles
    glRasterPos2f(0.6, 1.0-dv*(i+2)+v_offset);
    sprintf(string, "%d", rd_mov[i].num_tri);
    len = strlen(string);
    for (j=0; j<len; j++) {
      glutBitmapCharacter(GLUT_BITMAP_9_BY_15, string[j]);
    }
    // draw displaying status
    if (mode == GL_SELECT) glLoadName(i);
    glBegin(GL_QUADS);
    if (rd_mov[i].flag_display) {
      glColor3f(0.8, 0.2, 0.2);
    }
    else glColor3f(0, 0, 0);
    glVertex3f(0.86, 1.0-dv*(i+2)+v_offset, 0.0);
    glVertex3f(0.89, 1.0-dv*(i+2)+v_offset, 0.0);
    glVertex3f(0.89, 1.0-dv*(i+1.6)+v_offset, 0.0);
    glVertex3f(0.86, 1.0-dv*(i+1.6)+v_offset, 0.0);
    glEnd();
  }
}


void ReshapeMainAnc(int width, int height)
{
  double aspect;

  // reset gluPerspective to avoid distortion in W->V mapping
  aspect = (double)width / (double)height;
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPerspective(45.0, aspect, r_anc/10.0, r_anc*5);

  // reset viewport transformation
  glViewport(0, 0, width, height);

  W_main_anc = width; H_main_anc = height;
}


void ReshapeMainMov(int width, int height)
{
  double aspect;

  // reset gluPerspective to avoid distortion in W->V mapping
  aspect = (double)width / (double)height;
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPerspective(45.0, aspect, r_mov/10.0, r_mov*5);

  // reset viewport transformation
  glViewport(0, 0, width, height);

  W_main_mov = width; H_main_mov = height;
}


void ReshapeMainRes(int width, int height)
{
  double aspect;

  // reset gluPerspective to avoid distortion in W->V mapping
  aspect = (double)width / (double)height;
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPerspective(45.0, aspect, r_anc/10.0, r_anc*5);

  // reset viewport transformation
  glViewport(0, 0, width, height);

  W_main_res = width; H_main_res = height;
}

void ReshapeMenuAnc(int width, int height)
{
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluOrtho2D(0, 1, 0, 1);

  // reset viewport transformation
  glViewport(0, 0, Wmenu_anc, Hmenu_anc);
}

void ReshapeMenuMov(int width, int height)
{
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluOrtho2D(0, 1, 0, 1);

  // reset viewport transformation
  glViewport(0, 0, Wmenu_mov, Hmenu_mov);
}


void MouseMainAnc(int button, int state, int x, int y)
{
  double  MV[16], P[16];
  int     viewport[4];
  double  x1, y1, z1, x2, y2, z2;

  spinning_anc = panning_anc = scaling_anc = rotating_anc = 0;
  
  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
    beginx = x;
    beginy = y;
    glutIdleFunc(NULL);
    // shift button is down
    if (glutGetModifiers() & GLUT_ACTIVE_SHIFT) 
      scaling_anc = 1;
    else
      rotating_anc = 1;
  }
  else if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN) {
    // shift button is down --> picking mode
    if (glutGetModifiers() & GLUT_ACTIVE_SHIFT) {
      if (num_corres_anc >= MAX_NUM_CORRES) 
        return;
      glGetIntegerv(GL_VIEWPORT, viewport);
      glGetDoublev(GL_MODELVIEW_MATRIX, MV);
      glGetDoublev(GL_PROJECTION_MATRIX, P);
      gluUnProject((double)x, (double)(viewport[3] - y - 1), 0.0, 
                   MV, P, viewport, &x1, &y1, &z1);
      gluUnProject((double)x, (double)(viewport[3] - y - 1), 1.0, 
                   MV, P, viewport, &x2, &y2, &z2);
  
      PickCorresAnc(x1, y1, z1, x2, y2, z2);  
      glutPostRedisplay();
    }
    else {
      panning_anc = 1;
      beginx = x;
      beginy = y;
      glutIdleFunc(NULL);
    }
  }
  //else if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
  //  spinning = 1;
  //  glutIdleFunc(Spin);
  //}
}

void MouseMainMov(int button, int state, int x, int y)
{
  double  MV[16], P[16];
  int     viewport[4];
  double  x1, y1, z1, x2, y2, z2;

  spinning_mov = panning_mov = scaling_mov = rotating_mov = 0;
  
  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
    beginx = x;
    beginy = y;
    glutIdleFunc(NULL);
    // shift button is down
    if (glutGetModifiers() & GLUT_ACTIVE_SHIFT) 
      scaling_mov = 1;
    else
      rotating_mov = 1;
  }
  else if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN) {
    // shift button is down --> picking mode
    if (glutGetModifiers() & GLUT_ACTIVE_SHIFT) {
      if (num_corres_mov >= MAX_NUM_CORRES) 
        return;
      glGetIntegerv(GL_VIEWPORT, viewport);
      glGetDoublev(GL_MODELVIEW_MATRIX, MV);
      glGetDoublev(GL_PROJECTION_MATRIX, P);
      gluUnProject((double)x, (double)(viewport[3] - y - 1), 0.0, 
                   MV, P, viewport, &x1, &y1, &z1);
      gluUnProject((double)x, (double)(viewport[3] - y - 1), 1.0, 
                   MV, P, viewport, &x2, &y2, &z2);
  
      PickCorresMov(x1, y1, z1, x2, y2, z2);  
      glutPostRedisplay();
    }
    else {
      panning_mov = 1;
      beginx = x;
      beginy = y;
      glutIdleFunc(NULL);
    }
  }
  //else if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
  //  spinning = 1;
  //  glutIdleFunc(Spin);
  //}
}


void MouseMainRes(int button, int state, int x, int y)
{
  spinning_res = panning_res = scaling_res = rotating_res = 0;
  
  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
    beginx = x;
    beginy = y;
    glutIdleFunc(NULL);
    // shift button is down
    if (glutGetModifiers() & GLUT_ACTIVE_SHIFT) 
      scaling_res = 1;
    else
      rotating_res = 1;
  }
  else if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN) {
    panning_res = 1;
    beginx = x;
    beginy = y;
    glutIdleFunc(NULL);
  }
  //else if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
  //  spinning = 1;
  //  glutIdleFunc(Spin);
  //}
}


void MouseMenuAnc(int button, int state, int x, int y)
{
  GLint   viewport[4];
  GLuint  select_buf[512];
  GLint   hits;
  int     picked_rdata;
  
  if (button != GLUT_LEFT_BUTTON || state != GLUT_DOWN) return;

  glGetIntegerv(GL_VIEWPORT, viewport);

  glSelectBuffer(512, select_buf);
  glRenderMode(GL_SELECT);

  glInitNames();
  glPushName(0);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  // create 5x5 pixel picking region near cursor locatoin
  gluPickMatrix((double)x, (double)(viewport[3]-y), 5.0, 5.0, viewport);
  gluOrtho2D(0, 1, 0, 1);
  DrawMenuAnc(GL_SELECT);

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glFlush();

  hits = glRenderMode(GL_RENDER);
  picked_rdata = ProcessHits(hits, select_buf);
  if (picked_rdata > -1)
    rd_anc[picked_rdata].flag_display = !rd_anc[picked_rdata].flag_display;
  glutPostRedisplay();
  glutSetWindow(window_main_anc);
  glutPostRedisplay();
}


void MouseMenuMov(int button, int state, int x, int y)
{
  GLint   viewport[4];
  GLuint  select_buf[512];
  GLint   hits;
  int     picked_rdata;
  
  if (button != GLUT_LEFT_BUTTON || state != GLUT_DOWN) return;

  glGetIntegerv(GL_VIEWPORT, viewport);

  glSelectBuffer(512, select_buf);
  glRenderMode(GL_SELECT);

  glInitNames();
  glPushName(0);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  // create 5x5 pixel picking region near cursor locatoin
  gluPickMatrix((double)x, (double)(viewport[3]-y), 5.0, 5.0, viewport);
  gluOrtho2D(0, 1, 0, 1);
  DrawMenuMov(GL_SELECT);

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glFlush();

  hits = glRenderMode(GL_RENDER);
  picked_rdata = ProcessHits(hits, select_buf);
  if (picked_rdata > -1)
    rd_mov[picked_rdata].flag_display = !rd_mov[picked_rdata].flag_display;
  glutPostRedisplay();
  glutSetWindow(window_main_mov);
  glutPostRedisplay();
}


int ProcessHits(int hits, GLuint buffer[])
{
  unsigned int i;
  GLuint  names, *ptr;

  ptr = (GLuint *)buffer;
  for (i=0; i<hits; i++) {
    names = *ptr;
    //printf("number of names for this hit = %d\n", names); 
    ptr++;
    //printf("z1 is %g; ", (float)*ptr/0x7fffffff); 
    ptr++;
    //printf("z2 is %g\n", (float)*ptr/0x7fffffff); 
    ptr++;
    fflush(stdout);
    return ((int)(*ptr));
  }
  return -1;
} 

void MouseMotionMainAnc(int x, int y)
{
  if (scaling_anc) {
    dist_anc += (beginy - y);
    beginx = x;
    beginy = y;
    new_model_anc = 1;
    glutPostRedisplay();
  }
  else if (panning_anc) {
    xpan_anc -= (beginx - x);
    ypan_anc += (beginy - y);
    beginx = x;
    beginy = y;
    new_model_anc = 1;
    glutPostRedisplay();
  }
  else if (rotating_anc) {
    trackball(lastquat_anc, (2.0*beginx - W_main_anc) / W_main_anc, 
                            (H_main_anc - 2.0*beginy) / H_main_anc,
                            (2.0*x - W_main_anc) / W_main_anc,      
                            (H_main_anc - 2.0*y) / H_main_anc);
    add_quats(lastquat_anc, curquat_anc, curquat_anc);
    beginx = x;
    beginy = y;
    new_model_anc = 1;
    glutPostRedisplay();
  }
}


void MouseMotionMainMov(int x, int y)
{
  if (scaling_mov) {
    dist_mov += (beginy - y);
    beginx = x;
    beginy = y;
    new_model_mov = 1;
    glutPostRedisplay();
  }
  else if (panning_mov) {
    xpan_mov -= (beginx - x);
    ypan_mov += (beginy - y);
    beginx = x;
    beginy = y;
    new_model_mov = 1;
    glutPostRedisplay();
  }
  else if (rotating_mov) {
    trackball(lastquat_mov, (2.0*beginx - W_main_mov) / W_main_mov, 
                            (H_main_mov - 2.0*beginy) / H_main_mov,
                            (2.0*x - W_main_mov) / W_main_mov,      
                            (H_main_mov - 2.0*y) / H_main_mov);
    add_quats(lastquat_mov, curquat_mov, curquat_mov);
    beginx = x;
    beginy = y;
    new_model_mov = 1;
    glutPostRedisplay();
  }
}

void MouseMotionMainRes(int x, int y)
{
  if (scaling_res) {
    dist_res += (beginy - y);
    beginx = x;
    beginy = y;
    new_model_res = 1;
    glutPostRedisplay();
  }
  else if (panning_res) {
    xpan_res -= (beginx - x);
    ypan_res += (beginy - y);
    beginx = x;
    beginy = y;
    new_model_res = 1;
    glutPostRedisplay();
  }
  else if (rotating_res) {
    trackball(lastquat_res, (2.0*beginx - W_main_res) / W_main_res, 
                            (H_main_res - 2.0*beginy) / H_main_res,
                            (2.0*x - W_main_res) / W_main_res,      
                            (H_main_res - 2.0*y) / H_main_res);
    add_quats(lastquat_res, curquat_res, curquat_res);
    beginx = x;
    beginy = y;
    new_model_res = 1;
    glutPostRedisplay();
  }
}


void VisMainAnc(int visible)
{
  if (visible == GLUT_VISIBLE) {
    if (spinning_anc) glutIdleFunc(SpinAnc);
  }
  else {
    if (spinning_anc) glutIdleFunc(NULL);
  }
}


void VisMainMov(int visible)
{
  if (visible == GLUT_VISIBLE) {
    if (spinning_mov) glutIdleFunc(SpinMov);
  }
  else {
    if (spinning_mov) glutIdleFunc(NULL);
  }
}

void VisMainRes(int visible)
{
  if (visible == GLUT_VISIBLE) {
    if (spinning_res) glutIdleFunc(SpinMov);
  }
  else {
    if (spinning_res) glutIdleFunc(NULL);
  }
}

void SpinAnc(void)
{
  // function from trackball.c
  add_quats(lastquat_anc, curquat_anc, curquat_anc);
  new_model_anc = 1;
  glutPostRedisplay(); 
}

void SpinMov(void)
{
  // function from trackball.c
  add_quats(lastquat_mov, curquat_mov, curquat_mov);
  new_model_mov = 1;
  glutPostRedisplay(); 
}


void MakeMenu(void)
{
  int regist_menu_id;
	int view_menu_id;
  int shademodel_menu_id;
  int shininess_menu_id;

  regist_menu_id = glutCreateMenu(RegistrationMenu);
  glutAddMenuEntry("Clear last point", 1);
	glutAddMenuEntry("Clear all points", 2);
	glutAddMenuEntry("Course registration", 3);
	glutAddMenuEntry("Fine registration", 4);
  glutAddMenuEntry("Save result", 5);

	view_menu_id = glutCreateMenu(ViewMenu);
	glutAddMenuEntry("Mesh", 1);
	glutAddMenuEntry("Wireframe", 2);
	glutAddMenuEntry("Point", 3);
	glutAddMenuEntry("Normal", 4);
  glutAddMenuEntry("Texture", 5);

  shademodel_menu_id = glutCreateMenu(ShadeModelMenu);
  glutAddMenuEntry("Flat", 1);
  glutAddMenuEntry("Smooth", 2);

  shininess_menu_id = glutCreateMenu(ShininessMenu);
  glutAddMenuEntry("Dull", 1);
  glutAddMenuEntry("Medium", 2);
  glutAddMenuEntry("Shiny", 3);

  glutCreateMenu(MainMenu);
  glutAddSubMenu("Registration", regist_menu_id);
  glutAddSubMenu("Viewing Option", view_menu_id);
  glutAddSubMenu("Shade Model", shademodel_menu_id);
  glutAddSubMenu("Shininess", shininess_menu_id);
  glutAddMenuEntry("Toggle light", 1);
  glutAddMenuEntry("Toggle axes", 2);
  glutAddMenuEntry("Toggle bounding box", 3);
  glutAddMenuEntry("Reset view", 4);
  glutAddMenuEntry("QUIT", 0);

	glutAttachMenu(GLUT_RIGHT_BUTTON);
}

void RegistrationMenu(int value)
{
  int win;

  win = glutGetWindow();

  switch (value) {
  case 1: // clear last point
    if (win == window_main_anc && num_corres_anc > 0) 
      num_corres_anc--;
    else if (win == window_main_mov && num_corres_mov > 0)
      num_corres_mov--;
    break;
  case 2: // clear all points
    if (win == window_main_anc)
      num_corres_anc = 0;
    else
      num_corres_mov = 0;
    break;
  case 3: // course registration
    if (num_corres_anc < 4 || num_corres_mov < 4) {
      printf("\nNumber of corresponding points must be at least 3\n");
      fflush(stdout);
    }
    else {
      CourseRegistration();
      glutSetWindow(window_main_res);
    }
    break;
  case 4: // fine registration
    FineRegistration();
    glutSetWindow(window_main_res);
    break;
  case 5: // save registration result
    SaveRegistration();
    DefineInitialView();
    break;
  }  
  glutPostRedisplay();
}

void ViewMenu(int value)
{
  int win;

  win = glutGetWindow();

  if (win == window_main_anc) {
    switch (value) {
    case 1: // triangle mesh
      display_option_anc = 1;
		  break;
    case 2:	// triangle edge
      display_option_anc = 2;
      break;
    case 3: // point
      display_option_anc = 3; 
      break;
    case 4: // normal
      display_option_anc = 4;
      break;
    case 5: // texture
      display_option_anc = 5;
      break;
    }
    DisplayMainAnc();
  }
  else {
     switch (value) {
    case 1: // triangle mesh
      display_option_mov = 1;
		  break;
    case 2:	// triangle edge
      display_option_mov = 2;
      break;
    case 3: // point
      display_option_mov = 3; 
      break;
    case 4: // normal
      display_option_mov = 4;
      break;
    case 5: // texture
      display_option_mov = 5;
      break;
    }
    DisplayMainMov();
  } 
}

void ShadeModelMenu(int value)
{
  switch (value) {
  case 1: // flat shading
    glShadeModel(GL_FLAT);
    break;
  case 2: // smooth shading
    glShadeModel(GL_SMOOTH);
    break;
  }
  DisplayMainAnc();
  DisplayMainMov();
}

void ShininessMenu(int value)
{
  int win;

  win = glutGetWindow();
  if (win == window_main_anc) {
    switch (value) {
    case 1: // dull
      shininess_anc = 2.0;
      break;
    case 2: // medium
      shininess_anc = 20.0;
      break;
    case 3: // shiny
      shininess_anc = 100.0;
      break;
    }
    DisplayMainAnc();
  }
  else {
    switch (value) {
    case 1: // dull
      shininess_mov = 2.0;
      break;
    case 2: // medium
      shininess_mov = 20.0;
      break;
    case 3: // shiny
      shininess_mov = 100.0;
      break;
    }
    DisplayMainMov();
  }
}

void MainMenu(int value)
{
  int win;

  win = glutGetWindow();

  if (win == window_main_anc) {
    switch (value) {
    case 0: 
      exit(0);
      break;
    case 1: // toggle lighting
      lighting_anc = !lighting_anc;
      break;
    case 2: // toggle axes
      draw_axes_anc = !draw_axes_anc;
      break;
    case 3: // toggle bonding box
      draw_bbox_anc = !draw_bbox_anc;
      break;
    case 4: // reset view
      InitializeViewAnc();
      break;
    }
  }
  else {
    switch (value) {
    case 0: 
      exit(0);
      break;
    case 1: // toggle lighting
      lighting_mov = !lighting_mov;
      break;
    case 2: // toggle axes
      draw_axes_mov = !draw_axes_mov;
      break;
    case 3: // toggle bonding box
      draw_bbox_mov = !draw_bbox_mov;
      break;
    case 4: // reset view
      InitializeViewMov();
      break;
    }
  }
}


