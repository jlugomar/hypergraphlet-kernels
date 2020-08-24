CC=g++
DEBUG=-g
#FLAGS=-Wall -O2
FLAGS=-Wall
LIBS=-lm
OBJS=utils.o string.o


all: run_hyperkernel

utils.o: utils.h utils.cc config.h
	$(CC) $(DEBUG) $(FLAGS) -c utils.cc

string.o: string.h string.cc
	$(CC) $(DEBUG) $(FLAGS) -c string.cc


simhypergraph.o: simhypergraph.h simhypergraph.cc string.h
	$(CC) $(DEBUG) $(FLAGS) -c simhypergraph.cc



run_hyperkernel.o: hyperkernel.h hyperkernel.cc config.h
	$(CC) $(DEBUG) $(FLAGS) -c hyperkernel.cc

run_hyperkernel: run_hyperkernel.cc hyperkernel.o simhypergraph.o $(OBJS) 
	$(CC) $(DEBUG) $(FLAGS) $(LIBS) run_hyperkernel.cc hyperkernel.o simhypergraph.o $(OBJS) -o $@



clean:
	rm -rf run_hyperkernel *.o *.*~ *~ core *.dSYM

