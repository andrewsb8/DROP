SHELL := /bin/bash
#PATH  := node_modules/.bin:$(PATH)
buildDir = $(pwd)

all: compile

compile:
	gcc -c src/dropinfo/dropinfo.c
	gcc -c src/drop/drop.c
	gcc -c src/drop/commands.c
	gcc -c src/utils/logging/logging.c
	gcc -c src/utils/readProtein/readProtein.c
	gcc -c src/utils/vectorCalculus/vectorCalculus.c
	gcc -c src/utils/dihedralRotation/dihedralRotation.c
	gcc -c src/utils/vdwEnergy/vdwEnergy.c
	gcc -c src/utils/stericClash/stericClash.c
	gcc -c src/utils/fileHandling/fileHandling.c
	gcc -c src/dropanalysis/setDihedral.c
	gcc -c src/dropanalysis/setDihedralList.c
	gcc -c src/dropanalysis/measureDihedrals.c
	gcc -c src/dropanalysis/stericClashes.c
	gcc -c src/dropanalysis/stericScan.c
	gcc -c src/dropanalysis/vdwScan.c
	gcc -c src/dropanalysis/vdwScanSC.c
	gcc -o drop dropinfo.o drop.o commands.o logging.o readProtein.o setDihedral.o setDihedralList.o measureDihedrals.o dihedralRotation.o vectorCalculus.o fileHandling.o vdwEnergy.o stericClashes.o stericClash.o stericScan.o vdwScan.o vdwScanSC.o -lm
	rm *.o

test:
	gcc -c src/dropinfo/dropinfo.c
	gcc -c src/utils/vectorCalculus/vectorCalculus.c
	gcc -c src/utils/readProtein/readProtein.c
	gcc -c src/utils/logging/logging.c
	gcc -c src/utils/dihedralRotation/dihedralRotation.c
	gcc -o test_drop tests/test_drop.c readProtein.o logging.o dihedralRotation.o vectorCalculus.o dropinfo.o -lcheck -lm -lrt -lsubunit
	rm *.o
	-./test_drop
	rm test_drop
