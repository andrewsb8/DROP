SHELL := /bin/bash
#PATH  := node_modules/.bin:$(PATH)
buildDir = $(pwd)

all: compile

compile:
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
	gcc -o drop drop.o commands.o logging.o readProtein.o setDihedral.o setDihedralList.o measureDihedrals.o dihedralRotation.o vectorCalculus.o fileHandling.o vdwEnergy.o stericClashes.o stericClash.o stericScan.o vdwScan.o -lm
	rm *.o


TESTFLAGS=-Wall -Wextra -Wconversion -Wredundant-decls -Wshadow -Wno-unused-parameter -Wno-unused-variable -Wno-unused-result -O3

test:
	gcc $(TESTFLAGS) -c tests/test_main.c
	gcc $(TESTFLAGS) -c tests/test_cases/mytests.c
	gcc -c src/utils/readProtein/readProtein.c
	gcc test_main.o mytests.o readProtein.o -o test_binary -lm
	-./test_binary
	rm *.o test_binary
