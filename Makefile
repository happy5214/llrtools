# Makefile for llrtools

CC = gcc
LIBS = -lm

all: av_time fft_len get_time

av_time: llrtools.o av_time.o
	$(CC) -o av_time llrtools.o av_time.o $(LIBS)

fft_len: llrtools.o fft_len.o
	$(CC) -o fft_len llrtools.o fft_len.o $(LIBS)

get_time: llrtools.o get_time.o
	$(CC) -o get_time llrtools.o get_time.o $(LIBS)

.c.o:
	$(CC) -c $<
