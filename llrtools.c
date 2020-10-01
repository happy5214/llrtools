/*********************************
 *  llrtools.c                   *
 *  v0.02  2005-10-02            *
 *********************************/

#include <stdio.h>
#include <math.h>
#include "llrtools.h"

int read_msecs_file()      /* read timings from file */
{
  FILE *ms_file;
  char string[100];
  long fft;
  float ms;
  int i, j;

  ms_file = fopen("times.txt","rt");
  fgets(string, 80, ms_file);

  printf("%s", string);     /* print header from file times.txt */

  for (i=0; i<nfft; i++)    /* initialize the msecs array */
    msecs[i] = -1.0;        /* set all to -1.0 == no timings available */

  i = 0;
  while (!feof(ms_file))
  {
    fgets(string, 80, ms_file);
    if ((strlen(string) > 1) && (!feof(ms_file)))
    {
      sscanf(string,"%ld %f", &fft, &ms);
      j = 0;
      while ((fft != fftlen[j]) && (j < nfft))  /* look for matching fftlen */
        j++;
      if (j == nfft)
        printf("*** WARNING: timings on fftlen = %ld will be ignored!\n", fft);
      else
        msecs[j] = (double) ms;
      i++;
    }
  }
  fclose(ms_file);
  return (i);
}

int read_maxlen_file()    /* read a list of FFT lengths and n_max (Mers.) */
{
  FILE *mx_file;
  char string[100];
  long fft, max_n;
  int i;

  i = 0;
  mx_file = fopen("maxlen.txt","rt");
  while (!feof(mx_file))
  {
    fgets(string, 80, mx_file);
    if ((strlen(string) > 1) && (!feof(mx_file)))
    {
      sscanf(string,"%ld %ld", &fft, &max_n);
      fftlen[i] = fft;
      n_mers[i] = max_n;
      i++;
    }
  }
  fclose(mx_file);
  return (i);
}

int fftlen_from_k_and_n(long k, long n)
/* find the appropriate fftlen for given k and n */
/* returns the index in array fftlen */
{
  long nmax;
  int i;
  i = 0;
  nmax = nmax_from_fftlen(k, fftlen[i], n_mers[i]);
  while (nmax < n)
  {
    i++;
    nmax = nmax_from_fftlen(k, fftlen[i], n_mers[i]);
  }
  return (i);
}

long nmax_from_fftlen(long k, long fftlen, long n_mersenne)
/* find the max. allowed n for given k and fftlen */
/* This is adapted from George Woltmans gwnum v24.14 */
{
  long nmax;
  double log2k;
  log2k = log((double) k) / log(2.0);

/* We need to decide whether zero-padded FFT is used or not.
   The algorithm used by gwnum is a bit tricky and needs
   more information than we already have here. 
   For simplicity we assume zero-padding for k > 2^20. 
   This might be wrong for a few cases in the range k = 1-1.3M */

  if (k < 1048576)        /* is k < 2^20 ??? */
    nmax = n_mersenne -= (long)(log2k + log2k*(double)(fftlen/2));
  else                    /* zero-padded FFT is used, if k > 2^20 */
    nmax = (long)(((double)n_mersenne + (double)fftlen*0.3)/2.0);

  nmax--;                 /* decrement nmax by one */
  return (nmax);
}

/* This is no longer used */
long nmax_from_fftlen_zeropad(long k, long fftlen, long n_mersenne)
{
  long nmax;
  nmax = (long)(((double)n_mersenne + (double)fftlen*0.3)/2.0);
  nmax--;
  return (nmax);
}

double compute_average_time(long k, long nmin, long nmax)
/* Get the average time for the range (nmin,nmax) with k fixed. */
{
  double t_average, T_tot;
  int i;
  int imin, imax;
  long n;
/* First, find the min. and max. FFT lengths */
  i = 0;
  n = nmax_from_fftlen(k, fftlen[i], n_mers[i]);
  while (n<nmin)
  {
    i++;
    n = nmax_from_fftlen(k, fftlen[i], n_mers[i]);
  }
  imin = i;
  while (n<nmax)
  {
    i++;
    n = nmax_from_fftlen(k, fftlen[i], n_mers[i]);
  }
  imax = i;
/*  printf ("%ld %ld\n", fftlen[imin], fftlen[imax]); */

/* Check, if we have timings for each of the necessary FFT lengths */
  for (i=imin; i<=imax; i++)
  {
    if (msecs[i] < 0.0)
    {
      printf("Sorry, there is no timing information for fftlen = %ld !\n", fftlen[i]);
      exit(2);
    }
  }

/* The total time is computed by integrating the function t(n).
   
       t(n) = n * iteration_time_at_fftlen(n)

   This function increases linearly with n for constant fftlen,
   but has steps at the FFT switching points.

   T_tot = (msecs(n_max)*n_max*n_max - msecs(n_min*n_min*n_min) / 2.0
          - Sum(n_i*n_i*(msecs(n_(i+1)) - msecs(n_i)) / 2.0) ,

   where the Sum runs over the FFT intervals and 
   msecs(n) is the time per iteration for given n.

   T_tot would be the total time, if we test 
   every n out of the interval (nmin,nmax).
   Therefore:     t_average = T_tot / (nmax - nmin)   */
  
  T_tot  = msecs[imax]*(double)nmax*(double)nmax;
  T_tot -= msecs[imin]*(double)nmin*(double)nmin;
  T_tot /= 2.0;

  for (i=imin; i<imax; i++)
  {
    n = nmax_from_fftlen(k, fftlen[i], n_mers[i]);
    T_tot -= (double)n*(double)n*(msecs[i+1]-msecs[i])/2.0;
  }
  T_tot /= 1000.0;      /*  convert msec -> sec */

  t_average = T_tot / (double)(nmax-nmin);

  if (k > 1048576)     /* different times in case of zero-padded FFT's */
    t_average *= 1.085;        /* use empirical scaling factor */
  return(t_average);
}

void generate_list(long k, long nmin, long nmax)
/* list the FFT lengths and switching points 
   for the interval (nmin,nmax) with k fixed */
{
  int i;
  int imin, imax;
  long n, nz;
  i = 0;
  n = nmax_from_fftlen(k, fftlen[i], n_mers[i]);
  while (n<nmin)
  {
    i++;
    n = nmax_from_fftlen(k, fftlen[i], n_mers[i]);
  }
  imin = i;
  while (n<nmax)
  {
    i++;
    n = nmax_from_fftlen(k, fftlen[i], n_mers[i]);
  }
  imax = i;
/*  printf ("%ld %ld\n", fftlen[imin], fftlen[imax]); */

  printf("\n");
  printf("The following FFT lengths would be used:\n\n");
  printf("    fftlen       nmax\n");
  printf("-----------------------\n");
  for (i=imin; i<=imax; i++)
  {
    n = nmax_from_fftlen(k, fftlen[i], n_mers[i]);
    printf("%10ld %10ld\n", fftlen[i], n);
  }
}
