/*********************************
 *  llrtools.c                   *
 *  v0.03  2008-12-03            *
 *********************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "llrtools.h"

int read_msecs_file(llrtools_data_t *data, char *file_name)      /* read timings from file */
{
  FILE *ms_file;
  char string[100];
  long fft;
  float ms;
  int i, j;

  ms_file = fopen(file_name, "rt");
  if (ms_file == NULL)
  {
    printf("msecs file %s not readable!\n", file_name);
    return 0;
  }

  fgets(string, 80, ms_file);

  //printf("%s", string);     /* print header from file times.txt */

  for (i=0; i<data->n_fft; i++)    /* initialize the msecs array */
    data->msecs[i] = -1.0;        /* set all to -1.0 == no timings available */

  i = 0;
  while (!feof(ms_file))
  {
    fgets(string, 80, ms_file);
    if ((strlen(string) > 1) && (!feof(ms_file)))
    {
      sscanf(string,"%ld %f", &fft, &ms);
      j = 0;
      while ((fft != data->fftlen[j]) && (j < data->n_fft))  /* look for matching fftlen */
        j++;
      if (j == data->n_fft)
        printf("*** WARNING: timings on fftlen = %ld will be ignored!\n", fft);
      else
        data->msecs[j] = (double) ms;
      i++;
    }
  }
  fclose(ms_file);
  data->n_times = i;
  return 1;
}

int read_maxlen_file(llrtools_data_t *data, char *file_name)    /* read a list of FFT lengths and n_max (Mers.) */
{
  FILE *mx_file;
  char string[100];
  long fft, max_n;
  int i;

  i = 0;
  mx_file = fopen(file_name, "rt");
  if (mx_file == NULL)
  {
    printf("maxlen file %s not readable!\n", file_name);
    return 0;
  }

  while (!feof(mx_file))
  {
    fgets(string, 80, mx_file);
    if ((strlen(string) > 1) && (!feof(mx_file)))
    {
      sscanf(string,"%ld %ld", &fft, &max_n);
      data->fftlen[i] = fft;
      data->n_mers[i] = max_n;
      i++;
    }
  }
  fclose(mx_file);
  data->n_fft = i;
  return 1;
}

static long nmax_from_fftlen(long k, long fftlen, long n_mersenne)
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
    nmax = n_mersenne - (long)(log2k + log2k*(double)(fftlen/2));
  else                    /* zero-padded FFT is used, if k > 2^20 */
    nmax = (long)(((double)n_mersenne + (double)fftlen*0.3)/2.0);

  nmax--;                 /* decrement nmax by one */
  return (nmax);
}

int fftlen_from_k_and_n(llrtools_data_t *data, long k, long n)
/* find the appropriate fftlen for given k and n */
/* returns the index in array fftlen */
{
  long nmax;
  int i;
  i = 0;
  nmax = nmax_from_fftlen(k, data->fftlen[i], data->n_mers[i]);
  while (nmax < n)
  {
    i++;
    nmax = nmax_from_fftlen(k, data->fftlen[i], data->n_mers[i]);
  }
  return (i);
}

double compute_average_time(llrtools_data_t *data, long k, long nmin, long nmax)
/* Get the average time for the range (nmin,nmax) with k fixed. */
{
  double t_average, T_tot;
  int i;
  int imin, imax;
  long n;
/* First, find the min. and max. FFT lengths */
  i = 0;
  n = nmax_from_fftlen(k, data->fftlen[i], data->n_mers[i]);
  while (n<nmin)
  {
    i++;
    n = nmax_from_fftlen(k, data->fftlen[i], data->n_mers[i]);
  }
  imin = i;
  while (n<nmax)
  {
    i++;
    n = nmax_from_fftlen(k, data->fftlen[i], data->n_mers[i]);
  }
  imax = i;
/*  printf ("%ld %ld\n", fftlen[imin], fftlen[imax]); */

/* Check, if we have timings for each of the necessary FFT lengths */
  for (i=imin; i<=imax; i++)
  {
    if (data->msecs[i] < 0.0)
    {
      printf("Sorry, there is no timing information for fftlen = %ld !\n", data->fftlen[i]);
      exit(2);
    }
  }

/* The total time is computed by integrating the function t(n).
   
       t(n) = n * iteration_time_at_fftlen(n)

   This function increases linearly with n for constant fftlen,
   but has steps at the FFT switching points.

   T_tot = (msecs(n_max)*n_max*n_max - msecs(n_min)*n_min*n_min) / 2.0
          - Sum(n_i*n_i*(msecs(n_(i+1)) - msecs(n_i)) / 2.0) ,

   where the Sum runs over the FFT intervals and 
   msecs(n) is the time per iteration for given n.

   T_tot would be the total time, if we test 
   every n out of the interval (nmin,nmax).
   Therefore:     t_average = T_tot / (nmax - nmin)   */
  
  T_tot  = data->msecs[imax]*(double)nmax*(double)nmax;
  T_tot -= data->msecs[imin]*(double)nmin*(double)nmin;
  T_tot /= 2.0;

  for (i=imin; i<imax; i++)
  {
    n = nmax_from_fftlen(k, data->fftlen[i], data->n_mers[i]);
    T_tot -= (double)n*(double)n*(data->msecs[i+1]-data->msecs[i])/2.0;
  }
  T_tot /= 1000.0;      /*  convert msec -> sec */

  t_average = T_tot / (double)(nmax-nmin);

  if (k > 1048576)     /* different times in case of zero-padded FFT's */
    t_average *= 1.085;        /* use empirical scaling factor */
  return(t_average);
}

void generate_list(llrtools_data_t *data, long k, long nmin, long nmax)
/* list the FFT lengths and switching points 
   for the interval (nmin,nmax) with k fixed */
{
  int i;
  int imin, imax;
  long n, nz;
  i = 0;
  n = nmax_from_fftlen(k, data->fftlen[i], data->n_mers[i]);
  while (n<nmin)
  {
    i++;
    n = nmax_from_fftlen(k, data->fftlen[i], data->n_mers[i]);
  }
  imin = i;
  while (n<nmax)
  {
    i++;
    n = nmax_from_fftlen(k, data->fftlen[i], data->n_mers[i]);
  }
  imax = i;
/*  printf ("%ld %ld\n", fftlen[imin], fftlen[imax]); */

  printf("\n");
  printf("The following FFT lengths would be used:\n\n");
  printf("    fftlen       " IT_ON "n_max" IT_OFF "\n");
  printf("-----------------------\n");
  for (i=imin; i<=imax; i++)
  {
    n = nmax_from_fftlen(k, data->fftlen[i], data->n_mers[i]);
    printf("%10ld %10ld\n", data->fftlen[i], n);
  }
}

llrtools_times_t get_times(char *llr_file_name, char *maxlen_file_name, char *msecs_file_name)
{
  FILE *infile;
  long k, n;
  int i;
  char string[100];
  llrtools_data_t data;
  llrtools_times_t time_data;

  time_data.term_count = -1L;

  if (!read_maxlen_file(&data, maxlen_file_name) || !read_msecs_file(&data, msecs_file_name))
    return time_data;

  time_data.term_count = 0L;
  infile = fopen(llr_file_name, "rt");
  if (infile == NULL)
    return time_data;

  fgets(string, 80, infile);     /* dummy read header information */

  time_data.average_time = 0.0;
  time_data.total_time = 0.0;

  while (!feof(infile))
  {
    fgets(string, 80, infile);
    if ((strlen(string) > 1) && (!feof(infile)))
    {
      sscanf(string,"%ld %ld", &k, &n);
      i = fftlen_from_k_and_n(&data, k, n);
      if (k < 1048576)                                      /* no zero-padding */
        time_data.total_time += (double)n * data.msecs[i]/1000.0;
      else                                                  /* with zero-padding */
        time_data.total_time += (double)n * 1.085 * data.msecs[i]/1000.0;  /* approx. scaling factor */
      time_data.term_count++;
    }
  }

  time_data.average_time = time_data.total_time / time_data.term_count;

  return time_data;
}
