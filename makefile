SHELL := /bin/bash
#PATH  := node_modules/.bin:$(PATH)
buildDir = $(pwd)

compile:
	gcc -c src/drop/drop.c
	gcc -c src/drop/commands.c
	gcc -c src/include/exceptions/exceptions.c
	gcc -c src/include/readProtein/readProtein.c
	gcc -c src/include/vectorCalculus/vectorCalculus.c
	gcc -c src/include/dihedralRotation/dihedralRotation.c
	gcc -c src/include/vdwEnergy/vdwEnergy.c
	gcc -c src/include/stericClash/stericClash.c
	gcc -c src/include/fileHandling/fileHandling.c
	gcc -c src/dropanalysis/setDihedral.c
	gcc -c src/dropanalysis/setDihedralList.c
	gcc -c src/dropanalysis/measureDihedrals.c
	gcc -c src/dropanalysis/stericClashes.c
	gcc -c src/dropanalysis/stericScan.c
	gcc -c src/dropanalysis/vdwScan.c
	gcc -o drop drop.o commands.o exceptions.o readProtein.o setDihedral.o setDihedralList.o measureDihedrals.o dihedralRotation.o vectorCalculus.o fileHandling.o vdwEnergy.o stericClashes.o stericClash.o stericScan.o vdwScan.o -lm
	rm *.o

all: compile
