CC = gcc
CFLAGS = -g  -I$$HOME/include
FFLAGS = -g
LDFLAGS= -L$$HOME/lib -lm -lkpar 
BIN = $$HOME/bin

py_link: py_link.o
	${CC} ${CFLAGS} py_link.o \
	${LDFLAGS} -o py_link

py_read_verner: py_read_verner.o indexx.o nrutil.o 
	${CC} ${CFLAGS} py_read_verner.o indexx.o nrutil.o $(LDFLAGS) -o py_read_verner

permitted: permitted.o
	${CC} ${CFLAGS} permitted.o $(LDFLAGS) -o permitted

