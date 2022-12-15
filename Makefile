CC = gcc

CFLAGS = -Wall -Ofast -funroll-loops

LIBS = -lm

us: main.c
	$(CC) $(CFLAGS) main.c -o output $(LIBS)
