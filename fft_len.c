/*********************************************
 *  fft_len.c                                *
 *  generate a list of the used FFT lengths  *
 *  for k fixed and (nmin,nmax) range        * 
 *********************************************/

#include <stdio.h>
#include "llrtools.h"

int main(int argc, char *argv[])
{
  long k, nmin, nmax;
  int i;
  llrtools_data_t data;

  data.errors[0] = '\0';

  char maxlen_filename[30];
  generate_data_file_name(maxlen_filename, 30, "maxlen", argc - 1, argv + 1);

  if (!read_maxlen_file(&data, maxlen_filename)) {    /* init fftlen and nmers arrays */
    return 1;
  }

  printf(IT_ON "k" IT_OFF " = "); scanf("%ld", &k);

  printf(IT_ON "n" IT_OFF " (min) = "); scanf("%ld", &nmin);
  printf(IT_ON "n" IT_OFF " (max) = "); scanf("%ld", &nmax);

  generate_list(&data, k, nmin, nmax);

  return (0);
}
