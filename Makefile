CC=g++
CXXFLAGS = -Wall -Wextra -std=c++11 -pedantic -g

PROG = testCluster
SRC = $(PROG).cpp libCluster.cpp
HDR = cluster.h
OBJ = $(SRC:.c=.o)

$(PROG): $(OBJ)
$(OBJ): $(HDR)
TAGS: $(SRC) $(HDR)
	etags &^


.PHONY: clean
clean:
	$(RM) $(OBJ)


