README file for LLRtools by Thomas Ritschel
-------------------------------------------

This is just a collection of routines to get information about FFT lengths used by LLR
and to estimate running times (total and averaged).  Originally made for the 321search
project, it has been modified to work also for general k.  Please feel free to further
development.  Perhaps Jean Penne and/or George Woltman could add some of the features
to future versions of LLR/PRP.


General information:
--------------------

There are the following files:

   llrtools.c  - the tools collection
   llrtools.h  - header file for llrtools.c

   fft_len.c   - lists the used FFT lengths for given (nmin,nmax) range (fixed k)
   av_time.c   - estimates the average time per LLR test (fixed k)
   get_time.c  - estimates total and average times for given LLR input file

The following two files are also required:

   maxlen.txt  - table of FFT lengths and associated max. Mersenne exponents (nmers)
   times.txt   - list of timings (msecs) for different FFT lengths

Depending on your CPU (SSE2 or non-SSE2), you will need to copy either "maxlen_P4.txt"
or "maxlen_Athlon.txt" to "maxlen.txt".  (The other two files "maxlen_P4_proth.txt" and
"maxlen_Athlon_proth.txt" contain information about the FFT lengths and assouiated
max. Mersenne exponents for Proth tests, e.g. k*2^n+1, but this isn't tested so far...)

To get correct timings you will need to parametrize the file times.txt for your
specific hardware.  There are examples given for two Athlons (1 and 2 GHz) and Opteron.
There is no need to get timings for any kind of FFT length.  It's sufficient to have
data on those FFT lengths covered by your range of n.

For small and medium sized k (typically k < 2^20), you should get very accurate results
on the FFT lengths and the timings.  In cases where LLR uses zero-padding (which is the
case for very large k), the results may be wrong by about 10 percent, since the algorithm
used by George Woltmans gwnum is much more sophisticated and needs additional information,
we do not have here.


On the computation of total and average times:
----------------------------------------------

If we assume uniform sampling of our candidates in a given range of n, the total time
can be approximated by integrating the function t(n).

                    1           n_max
  t_average = ------------- * Integral (t(n) dn)
              n_max - n_min    n=n_min

Generally t(n) is neither a parabola nor a linear function.  Rather than it is a 
non-steady function, partially linear but having steps at the FFT turnover points.
We can get the area under this function by integrating it partially.

For fixed FFT length we'll have:

  t(n) = n * msecs_per_iteration(n) ,

where msecs_per_iteration(n) is a constant value, depending on the FFT length.
We will substitute this by m(i), where i is an index for the FFT length.

By splitting the integral into parts at the FFT turnover points we get the following
(after integration):

                    1
  t_average = ----------------- * [ m(i_max)*n_max^2 - m(i_min)*n_min^2 
              2 (n_max - n_min)

                                       i_max
                                    -   Sum    [ (m(i+1) - m(i)) * n(i)^2 ] ]
                                      i=i_min

Here, n(i) is the max. n for FFT length i, i_min and i_max are indexes of the min. and
max. FFT lengths used.

Therefore we need information at which n the FFT length will change as well as the 
times per iteration for each of the FFT lengths.


Typical usage:
--------------

- These are command line tools, so run them manually from the command line rather than
  double clicking on them.  You would'nt see any screen output in the later case.

- First run fft_len to get some insight on the FFT lengths for the (nmin,nmax) range
  of your interest.

     fft_len

  Sample screen output:

     k = 3
     n(min) = 2000000
     n(max) = 5000000

     The following FFT lengths would be used:

         fftlen       nmax
     -----------------------
         114688    2233110
         131072    2560126
         163840    3180158
         196608    3777190
         229376    4411222
         262144    5056254


- Then check, if all the necessary FFT lengths are covered by the file "times.txt".
  If not, get the msecs/iteration from runing LLR for a few seconds or minutes on some
  appropriate (k,n) pairs.

- Run av_time, to get the average time for a given range of n.  If there is insufficient
  information on the timings you will get an error.

     av_time

  Sample screen output:

     --- Athlon 1GHz ---
     k = 3
     n(min) = 2000000
     n(max) = 5000000

     average time per LLR test: 90019.175 sec


- If you already have an input file (perhaps partially sieved) for LLR, you may want to 
  see, how long it will take to run LLR on the whole file.
  This information can be obtained by running get_time.  This program needs the file name
  of your sieve file as a command line parameter.  So simply run it like the following:

     get_time input.txt          (replace "input.txt" by your file name)

  Sample screen output:

     --- Athlon 1GHz ---
     number of (k,n) pairs in file: 9357
     estimated total time for LLR testing the whole file: 100376901.970 sec
     average time per LLR test: 10727.466 sec
