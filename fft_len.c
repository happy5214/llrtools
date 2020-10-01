/*********************************************
 *  fft_len.c                                *
 *  generate a list of the used FFT lengths  *
 *  for k fixed and (nmin,nmax) range        * 
 *********************************************/

#include <stdio.h>
#include "llrtools.h"

int main(void)
{
  long k, nmin, nmax;
  int i;
  llrtools_data_t data;

  if (!read_maxlen_file(&data, "maxlen.txt")) {    /* init fftlen and nmers arrays */
    return 1;
  }

  printf("k = "); scanf("%ld", &k);

  printf("n(min) = "); scanf("%ld", &nmin);
  printf("n(max) = "); scanf("%ld", &nmax);

  generate_list(&data, k, nmin, nmax);

  return (0);
}
