# This is the makefile for kpar which is ksl's io routine library
#
# Notes - This make package assumes ksl's standard directory environment.
#   for code development. 
#
# make installs the subroutine library and the header files in LIB and
# INCLUDE respectively.  A separate install step is not used.
# 
# History
# 	07jul	ksl 	Separate making prototypes from the standardd
# 			compilations so that it could be used on
# 			systems without the utility cproto.  

LIB = ../../lib
INCLUDE = ../../include


CC = gcc
CFLAGS = -g -Wall


FILES = \
rdpar.o log.o lineio.o

SOURCE = \
rdpar.c log.c lineio.c


libkpar.a: $(FILES)
	ar ru libkpar.a $(FILES)
	ranlib libkpar.a
	mv libkpar.a $(LIB)
	cp log.h $(INCLUDE)

prototypes: $(SOURCE)
	cproto $(SOURCE) >log.h

clean:
	rm *.o *~
# Test routimes

t_lineio: $(FILES) t_lineio.o
	gcc t_lineio.o $(FILES) -o Test/t_lineio

t_elog: t_elog.o rdpar.o log.o lineio.o
	gcc t_elog.o rdpar.o log.o lineio.o -o Test/t_elog
