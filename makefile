SHELL := /bin/bash
#PATH  := node_modules/.bin:$(PATH)
buildDir = $(pwd)

compile:
	gcc -c src/drop/drop.c
	gcc -c src/drop/commands.c
	gcc -c src/include/readProtein/readProtein.c
	gcc -c src/include/vectorCalculus/vectorCalculus.c
	gcc -c src/include/dihedralRotation/dihedralRotation.c
	gcc -c src/include/stericClash/stericClash.c
	gcc -c src/include/fileHandling/fileHandling.c
	gcc -c src/dropanalysis/setDihedral.c
	gcc -c src/dropanalysis/measureDihedral.c
	gcc -c src/dropanalysis/stericClashes.c
	gcc -o drop drop.o commands.o readProtein.o setDihedral.o measureDihedral.o dihedralRotation.o vectorCalculus.o fileHandling.o stericClashes.o stericClash.o -lm
	rm *.o

all: compile
