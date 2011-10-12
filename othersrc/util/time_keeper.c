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
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "time_keeper.h"

typedef struct time_keeper {
  clock_t begin_clock;
  clock_t save_clock;
  time_t begin_time;
  time_t save_time;
} time_keeper;


static time_keeper  tk;

void StartTime(void)
{
  tk.begin_clock = tk.save_clock = clock();
  tk.begin_time = tk.save_time = time(NULL);
}

void PrintTime(void)
{
  //char    s1[MAXSTRING], s2[MAXSTRING];
  //int     field_width, n1, n2;
  double  clocks_per_second = (double) CLOCKS_PER_SEC,
          user_time, real_time;

  user_time = (clock() - tk.save_clock) / clocks_per_second;
  real_time = difftime(time(NULL), tk.save_time);
  tk.save_clock = clock();
  tk.save_time = time(NULL);

/*
  n1 = sprintf(s1, "%.2f", user_time);
  n2 = sprintf(s2, "%.2f", real_time);
  field_width = (n1 > n2) ? n1 : n2;
  printf("%s%*.2f%s\n%s%*.2f%s\n\n",
         "User time: ", field_width, user_time, " seconds",
         "Real time: ", field_width, real_time, " seconds" );
*/         
  printf("*** User time: %.2f seconds ***\n\n", user_time);
  fflush(stdout);
}
