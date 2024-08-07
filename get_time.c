/**********************************************
 *  get_time.c                                *
 *  estimate total and average time           *
 *  for LLR input file                        *
 *                                            *
 *  input file header information is ignored  *
 *  we simply assume k*2^n-1                  *
 **********************************************/

#include <stdio.h>
#include <stdlib.h>
#include "llrtools.h"

int main(int argc, char *argv[])
{
  llrtools_times_t time_data;

  if (argc < 2)
  {
     printf("usage: %s <LLR-inputfile> [device]\n", argv[0]);
     exit(1);
  }

  char maxlen_filename[30];
  char times_filename[30];
  generate_data_file_name(maxlen_filename, 30, "maxlen", argc - 2, argv + 2);
  generate_data_file_name(times_filename, 30, "times", argc - 2, argv + 2);

  time_data = get_times(argv[1], maxlen_filename, times_filename);
  if (time_data.term_count < -1) {
    return 1;
  } else if (time_data.term_count == 0) {
    fprintf(stderr, "empty input file %s\n", argv[1]);
    return 2;
  }

  printf("number of (" IT_ON "k,n" IT_OFF ") pairs in file: %ld\n", time_data.term_count);
  printf("estimated total time for LLR testing the whole file: %.3f sec\n", time_data.total_time);
  printf("average time per LLR test: %.3f sec\n", time_data.average_time);

  return 0;
}
