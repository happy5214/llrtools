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

  nfft = read_maxlen_file();          /* initialize fftlen and nmers arrays */

  printf("k = "); scanf("%ld", &k);

  printf("n(min) = "); scanf("%ld", &nmin);
  printf("n(max) = "); scanf("%ld", &nmax);

  generate_list(k, nmin, nmax);

  return (0);
}
