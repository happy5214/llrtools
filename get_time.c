/**********************************************
 *  get_time.c                                *
 *  estimate total and average time           *
 *  for LLR input file                        *
 *                                            *
 *  input file header information is ignored  *
 *  we simply assume k*2^n-1                  *
 **********************************************/

#include <stdio.h>
#include "llrtools.h"

int main(int argc, char *argv[])
{
  FILE *infile;
  double average_time, total_time;
  long k, n;
  long count;        /* counter for (k,n) pairs */
  int n_times;       /* number of (fftlen, msecs) pairs */
  int i;
  char string[100];

  nfft = read_maxlen_file();
  n_times = read_msecs_file();

  if (argc < 2)
  {
     printf("usage: %s <LLR-inputfile>\n", argv[0]);
     exit(1);
  }
  infile = fopen(argv[1], "rt");
  if (infile == NULL)
  {
    printf("file %s not in current directory !\n", argv[1]);
    exit(2);
  }

  fgets(string, 80, infile);     /* dummy read header information */

  average_time = 0.0;
  total_time = 0.0;
  count = 0L;
  
  while (!feof(infile))
  {
    fgets(string, 80, infile);
    if ((strlen(string) > 1) && (!feof(infile)))
    {
      sscanf(string,"%ld %ld", &k, &n);
      i = fftlen_from_k_and_n(k, n);
      if (k < 1048576)                                      /* no zero-padding */
        total_time += (double)n * msecs[i]/1000.0;
      else                                                  /* with zero-padding */
        total_time += (double)n * 1.085 * msecs[i]/1000.0;  /* approx. scaling factor */
      count++;
    }
  }
  
  average_time = total_time / count;

  printf("number of (k,n) pairs in file: %ld\n", count);
  printf("estimated total time for LLR testing the whole file: %.3f sec\n", total_time);
  printf("average time per LLR test: %.3f sec\n", average_time);

  if (k > 1048576)          /* in case of zero-padding */
    printf("Warning: estimated times may be wrong by about 10 percent!\n");

  return (0);
}
