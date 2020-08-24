CC=g++
DEBUG=-g
#FLAGS=-Wall -O2
FLAGS=-Wall
LIBS=-lm
OBJS=utils.o string.o


all: run_hypergraphkernel

utils.o: utils.h utils.cc config.h
	$(CC) $(DEBUG) $(FLAGS) -c utils.cc

string.o: string.h string.cc
	$(CC) $(DEBUG) $(FLAGS) -c string.cc


simhypergraph.o: simhypergraph.h simhypergraph.cc string.h
	$(CC) $(DEBUG) $(FLAGS) -c simhypergraph.cc



run_hypergraphkernel.o: hypergraphkernel.h hypergraphkernel.cc config.h
	$(CC) $(DEBUG) $(FLAGS) -c hypergraphkernel.cc

run_hypergraphkernel: run_hypergraphkernel.cc hypergraphkernel.o simhypergraph.o $(OBJS) 
	$(CC) $(DEBUG) $(FLAGS) $(LIBS) run_hypergraphkernel.cc hypergraphkernel.o simhypergraph.o $(OBJS) -o $@



clean:
	rm -rf run_hypergraphkernel *.o *.*~ *~ core *.dSYM

