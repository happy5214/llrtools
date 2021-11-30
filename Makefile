# Makefile for llrtools

CC = gcc
LIBS = -lm
FLAGS = -O2 -fPIC

objs = llrtools.o

libllrtools_so = libllrtools.so
libllrtools_objs = $(objs)

av_time_objs = av_time.o $(objs)
fft_len_objs = fft_len.o $(objs)
get_time_objs = get_time.o $(objs)

.PHONY: all libllrtools

all: av_time fft_len get_time libllrtools

libllrtools: $(libllrtools_so)

libllrtools.so: $(libllrtools_objs)
	$(CC) -shared -o $@ $(libllrtools_objs)

av_time: $(av_time_objs)
	$(CC) -o $@ $(av_time_objs) $(LIBS)

fft_len: $(fft_len_objs)
	$(CC) -o $@ $(fft_len_objs) $(LIBS)

get_time: $(get_time_objs)
	$(CC) -o $@ $(get_time_objs) $(LIBS)

%.o: %.c llrtools.h
	$(CC) -c -o $@ $< $(FLAGS)
