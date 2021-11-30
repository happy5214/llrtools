/********************************************
 *  av_time.c                               *
 *  estimate average time per LLR test      *
 *  for fixed k and given (nmin,nmax) range *
 ********************************************/
 
#include <stdio.h>
#include <string.h>

#include "llrtools.h"

int main(int argc, char *argv[])
{
  long k, nmin, nmax;
  double average_time;
  int n_times;      /* number of (fftlen, msecs) pairs */
  int i;
  llrtools_data_t data;

  char maxlen_filename[30];
  char times_filename[30];
  if (argc >= 2) {
    snprintf(maxlen_filename, 30, "maxlen/%s", argv[1]);
    snprintf(times_filename, 30, "times/%s", argv[1]);
  } else {
    strcpy(maxlen_filename, "maxlen/default");
    strcpy(times_filename, "times/default");
  }

  if (!read_maxlen_file(&data, maxlen_filename) ||    /* init fftlen and nmers arrays */
      !read_msecs_file(&data, times_filename)) {      /* read timings from file */
    return 1;
  }

  printf(IT_ON "k" IT_OFF " = "); scanf("%ld", &k);

  printf(IT_ON "n" IT_OFF " (min) = "); scanf("%ld", &nmin);
  printf(IT_ON "n" IT_OFF " (max) = "); scanf("%ld", &nmax);

  average_time = compute_average_time(&data, k, nmin, nmax);

  printf("average time per LLR test: %.3f sec\n", average_time);

  if (k > 1048576)          /* in case of zero-padding */
    printf("Warning: estimated time may be wrong by about 10 percent!\n");

  return 0;
}
