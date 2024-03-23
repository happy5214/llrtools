/**********************************************
 *  get_percentage.c                          *
 *  report percentage of FFT length usage     *
 *  for LLR input file                        *
 *                                            *
 *  input file header information is ignored  *
 *  we simply assume k*2^n-1                  *
 **********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "llrtools.h"

int main(int argc, char *argv[])
{
  FILE *infile;
  long k, n;
  long count;        /* counter for (k,n) pairs */
  long counts[100];
  int fcount;        /* counter for FFT lengths */
  int i;
  char string[100];

  nfft = read_maxlen_file();

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

  count = 0L;
  for (i=0; i<nfft; i++)
    counts[i] = 0L;
  
  while (!feof(infile))
  {
    fgets(string, 80, infile);
    if ((strlen(string) > 1) && (!feof(infile)))
    {
      sscanf(string,"%ld %ld", &k, &n);
      count++;
      i = fftlen_from_k_and_n(k, n);
      counts[i]++;
    }
  }
  
  fcount = 0;

  for (i=0; i<nfft; i++)
  {
    if (counts[i] > 0)
    {
      if (fcount > 0)
        printf(", ");
      printf("%ik: %4.1f%%", fftlen[i]/1024, 100*(double)(counts[i])/(double)(count));
      fcount++;
    }
  }
  printf("\n");

  return (0);
}
