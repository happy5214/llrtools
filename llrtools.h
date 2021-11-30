/* llrtools.h */

#ifndef LLRTOOLS_H
#define LLRTOOLS_H 1

#define MAX(a,b) ((a) > (b)? (a) : (b))

#define FFTLEN_SIZE 100
#define ERROR_SIZE 500

#define IT_ON "\033[3m"
#define IT_OFF "\033[m"

typedef struct {
	long fftlen[FFTLEN_SIZE];   /* list of FFT lengths */
	long n_mers[FFTLEN_SIZE];   /* list of n_max (Mersenne) for FFT lengths */
	double msecs[FFTLEN_SIZE];  /* list of times per iteration (msecs) for FFT lengths */
	int n_fft;                  /* number of FFT lengths in file maxlen.txt */
	int n_times;                /* number of (fftlen, msecs) pairs */
	char errors[ERROR_SIZE];    /* error messages */
} llrtools_data_t;

typedef struct {
	double average_time;
	double total_time;
	long term_count;
} llrtools_times_t;

int read_msecs_file(llrtools_data_t *data, char *file_name);
int read_maxlen_file(llrtools_data_t *data, char *file_name);
int fftlen_from_k_and_n(llrtools_data_t *data, long k, long n);
double compute_average_time(llrtools_data_t *data, long k, long nmin, long nmax);
void generate_list(llrtools_data_t *data, long k, long nmin, long nmax);
llrtools_times_t get_times(char *llr_file_name, char *maxlen_file_name, char *msecs_file_name);

#endif // LLRTOOLS_H
