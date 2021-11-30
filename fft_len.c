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

  printf(IT_ON "k" IT_OFF " = "); scanf("%ld", &k);

  printf(IT_ON "n" IT_OFF " (min) = "); scanf("%ld", &nmin);
  printf(IT_ON "n" IT_OFF " (max) = "); scanf("%ld", &nmax);

  generate_list(&data, k, nmin, nmax);

  return (0);
}
