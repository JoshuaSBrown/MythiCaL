CC=g++
VAR=
CFLAGS = -Wall -Wextra -std=c++11 -pedantic -g  $(VAR)

PROG = testCluster
SRC = $(PROG).cpp libCluster.cpp
HDR = cluster.h
OBJ = $(SRC:.cpp=.o)

$(PROG): $(OBJ)
$(OBJ): $(SRC)

.PHONY: clean
clean:
	 $(RM) $(OBJ)

all: testCluster

testCluster: testCluster.cpp libCluster.o
	$(CC) $(CFLAGS) -o testCluster libCluster.o testCluster.cpp

libCluster.o: libCluster.cpp
	$(CC) $(CFLAGS) -c libCluster.cpp
