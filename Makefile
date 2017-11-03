
CC=gcc

all:	pseudo-erosion

png_utils.o:	png_utils.c png_utils.h
	${CC} ${CFLAGS} -c png_utils.c

pseudo-erosion:	pseudo-erosion.c png_utils.o
	${CC} ${CFLAGS} -o pseudo-erosion png_utils.o pseudo-erosion.c -lpng

