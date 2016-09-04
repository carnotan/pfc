TARGET = paralelo
TARGET_DBG = paralelo_dbg

CC = gcc
NVCC = nvcc

CFLAGS = -Wall -O3 -std=c++11 -D_MWAITXINTRIN_H_INCLUDED
SRC = $(wildcard *.cu)
OBJ = ${SRC:.cu=.o}
OBJ_DBG = ${SRC:.cu=_dbg.o}

${TARGET}: $(OBJ)
	${NVCC} -o $@ $(OBJ)

${TARGET_DBG}: $(OBJ_DBG)
	${NVCC} -o $@ $(OBJ_DBG)

.SUFFIXES: .o .cu _dbg.o

.cu.o:
	$(NVCC) -ccbin=${CC} --compiler-options ${CFLAGS} -c $<

.cu_dbg.o:
	$(NVCC) -ccbin=${CC} --compiler-options ${CFLAGS} ${DBGFLAGS} -c -o $*_dbg.o $<

tidy:
	rm -f *.o *~

clean: tidy
	rm -f ${TARGET} ${TARGET_DBG}
