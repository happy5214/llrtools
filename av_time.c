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

  nfft = read_maxlen_file();               /* init fftlen and nmers arrays */
  n_times = read_msecs_file();             /* read timings from file */
  
  printf("k = "); scanf("%ld", &k);

  printf("n(min) = "); scanf("%ld", &nmin);
  printf("n(max) = "); scanf("%ld", &nmax);

  average_time = compute_average_time(k, nmin, nmax);

  printf("average time per LLR test: %.3f sec\n", average_time);
  
  if (k > 1048576)          /* in case of zero-padding */
    printf("Warning: estimated time may be wrong by about 10 percent!\n");
    
  return (0);
}
