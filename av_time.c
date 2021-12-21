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
  long k, n_min, n_max;
  double average_time;
  llrtools_data_t data;

  data.errors[0] = '\0';

  char maxlen_filename[30];
  char times_filename[30];
  generate_data_file_name(maxlen_filename, 30, "maxlen", argc - 1, argv + 1);
  generate_data_file_name(times_filename, 30, "times", argc - 1, argv + 1);

  if (!read_maxlen_file(&data, maxlen_filename) ||    /* init fftlen and nmers arrays */
      !read_msecs_file(&data, times_filename)) {      /* read timings from file */
    fputs(data.errors, stderr);
    return 1;
  }

  printf(IT_ON "k" IT_OFF " = "); scanf("%ld", &k);

  printf(IT_ON "n" IT_OFF " (min) = "); scanf("%ld", &n_min);
  printf(IT_ON "n" IT_OFF " (max) = "); scanf("%ld", &n_max);

  average_time = compute_average_time(&data, k, n_min, n_max);
  if (average_time < 0.0) {
    fputs(data.errors, stderr);
    return 2;
  }

  printf("average time per LLR test: %.3f sec\n", average_time);

  if (k > 1048576)          /* in case of zero-padding */
    printf("Warning: estimated time may be wrong by about 10 percent!\n");

  return 0;
}
