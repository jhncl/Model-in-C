# Makefile for building C stuff with GSL
CAPL=5
BURN=1
ITERS=1000
THIN=1

CFLAGS=-g -ansi -Wall -pedantic -I /usr/include/gsl
LDFLAGS=-lgsl -lgslcblas
CC=gcc 


EXEC=main

test.tab:main
	./main 1 1 1 CAPL> test.tab
	echo "Done.G^"

main:main.o headers.h datain.o datain.h functions.o functions.h print.o print.h 
	$(CC) $(CFLAGS) main.o datain.o functions.o print.o -o main $(LDFLAGS)

clean:
	rm -f core *.o *~ *.tab *.out *.ind *.dat $(EXEC)

edit:
	gnuclient README.txt Makefile *.h *.c &

print:
	a2ps README.txt Makefile *.h *.c
