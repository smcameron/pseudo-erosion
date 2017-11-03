
CC=gcc

all:	pseudo-erosion

CFLAGS=-g

open-simplex-noise.o:	open-simplex-noise.c open-simplex-noise.h
	${CC} ${CFLAGS} -c open-simplex-noise.c

png_utils.o:	png_utils.c png_utils.h
	${CC} ${CFLAGS} -c png_utils.c

pseudo-erosion:	pseudo-erosion.c png_utils.o open-simplex-noise.o
	${CC} ${CFLAGS} -o pseudo-erosion png_utils.o open-simplex-noise.o pseudo-erosion.c -lm -lpng

clean:
	rm -f *.o pseudo-erosion

