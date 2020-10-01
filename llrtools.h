/* llrtools.h */

long fftlen[100];   /* list of FFT lengths */
long n_mers[100];   /* list of n_max (Mersenne) for FFT lengths */
double msecs[100];  /* list of times per iteration (msecs) for FFT lengths */

int nfft;           /* number of FFT lengths in file maxlen.txt */

int read_msecs_file(void);
int read_maxlen_file(void);
long nmax_from_fftlen(long k, long fftlen, long n_mersenne);
long nmax_from_fftlen_zeropad(long k, long fftlen, long n_mersenne);
double compute_average_time(long k, long nmin, long nmax);
void generate_list(long k, long nmin, long nmax);
