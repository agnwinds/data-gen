INCLUDE = ../../include
LIB = ../../lib
BIN = ../../bin

CC = gcc
#CFLAGS = -xpg -xCC
# CFLAGS = -g -xCC -I$(INCLUDE)    
CFLAGS = -g -pg  -I$(INCLUDE)    
LDFLAGS= -L$(LIB) -lm -lkpar 

prototypes:
	cproto -I$(INCLUDE) py_read_kurucz.c indexx.c nrutil.c  > foo.h
	cproto -I$(INCLUDE) py_top_phot.c py_read_verner.c >>foo.h
	cp foo.h templates.h
	rm foo.h

#Read Kurucz line list and produce python line list and intermediate level file
py_read_kurucz: py_read_kurucz.o indexx.o nrutil.o
	${CC} ${CFLAGS} py_read_kurucz.o indexx.o nrutil.o \
	${LDFLAGS} -o py_read_kurucz
	mv $@ $(BIN)/py_read_kurucz

py_link: py_link.o
	${CC} ${CFLAGS} py_link.o \
	${LDFLAGS} -o py_link
	mv $@ $(BIN)/py_link

#Read Verner line list and produce python line list and intermediate level file
py_read_verner: py_read_verner.o indexx.o nrutil.o 
	${CC} ${CFLAGS} py_read_verner.o indexx.o nrutil.o $(LDFLAGS) -o py_read_verner
	mv $@ $(BIN)/py_read_verner


#Read topbase photoionization x sections, smooth and write in python format
py_top_phot: py_top_phot.o
	${CC} ${CFLAGS} py_top_phot.o ${LDFLAGS} -o py_top_phot
	mv $@ $(BIN)/py_top_phot

#Read topbase levels files and write in python format
py_top_levels: py_top_levels.o indexx.o nrutil.o
	${CC} ${CFLAGS} py_top_levels.o indexx.o nrutil.o ${LDFLAGS} -o py_top_levels
	mv $@ $(BIN)/py_top_levels


clean :
	rm -f *.o  *~

all: clean py_read_kurucz py_link py_read_verner py_top_phot py_top_levels
