/********************************************
 *  av_time.c                               *
 *  estimate average time per LLR test      *
 *  for fixed k and given (nmin,nmax) range *
 ********************************************/
 
#include <stdio.h>
#include "llrtools.h"

int main(void)
{
  long k, nmin, nmax;
  double average_time;
  int n_times;      /* number of (fftlen, msecs) pairs */
  int i;
  llrtools_data_t data;

  if (!read_maxlen_file(&data, "maxlen.txt") ||    /* init fftlen and nmers arrays */
      !read_msecs_file(&data, "times.txt")) {      /* read timings from file */
    return 1;
  }

  printf("k = "); scanf("%ld", &k);

  printf("n(min) = "); scanf("%ld", &nmin);
  printf("n(max) = "); scanf("%ld", &nmax);

  average_time = compute_average_time(&data, k, nmin, nmax);

  printf("average time per LLR test: %.3f sec\n", average_time);

  if (k > 1048576)          /* in case of zero-padding */
    printf("Warning: estimated time may be wrong by about 10 percent!\n");

  return (0);
}
