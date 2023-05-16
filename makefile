SHELL := /bin/bash
#PATH  := node_modules/.bin:$(PATH)

compile:
	gcc -c src/drop/drop.c
	gcc -c src/drop/commands.c
	gcc -c src/include/readProtein/readProtein.c
	gcc -c src/dropanalysis/trial.c
	gcc -o drop drop.o commands.o readProtein.o trial.o
	rm *.o

all: compile
