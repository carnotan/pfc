TARGET = serie
TARGET_DBG = serie_dbg

CC = g++
NVCC = nvcc

CFLAGS = -Wall -O3 -std=c++11
SRC = $(wildcard *.cpp)
OBJ = ${SRC:.cpp=.o}
OBJ_DBG = ${SRC:.cpp=_dbg.o}

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
