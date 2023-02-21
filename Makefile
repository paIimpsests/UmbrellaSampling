CC = gcc

CFLAGS = -Wall -Ofast

LIBS = -lm

all: umbrella energy

umbrella: umbrella.c Makefile
	$(CC) $(CFLAGS) umbrella.c -o umbrella $(LIBS)

energy: computeDGn.c Makefile
	$(CC) $(CFLAGS) computeDGn.c -o energy $(LIBS)

