CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

makeall: mainparallel.exe main.exe

mainparallel.exe : mainparallel.o random.o classes.o
	$(CC) random.o classes.o mainparallel.o -o mainparallel.exe
main.exe : main.o random.o classes.o
	g++ random.o classes.o main.o -o main.exe
classes.o: classes.cpp classes.h random.h vectoralgebra.h
	g++ -c classes.cpp -o classes.o $(CFLAGS)
mainparallel.o : mainparallel.cpp vectoralgebra.h classes.h
	$(CC) -c mainparallel.cpp -o mainparallel.o $(CFLAGS)
main.o : main.cpp vectoralgebra.h classes.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean_dat:
	find . -maxdepth 1 -name '*.dat' ! -name 'American_capitals.dat' -type f -delete
clean_o: 
	rm *.o 