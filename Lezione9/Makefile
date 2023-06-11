CC = g++
CFLAGS = -Wall -O3 --std=c++11

makeall: main.exe

main.exe : main.o random.o classes.o
	$(CC) random.o classes.o main.o -o main.exe
classes.o: classes.cpp classes.h random.h vectoralgebra.h
	g++ -c classes.cpp -o classes.o $(CFLAGS)
main.o : main.cpp vectoralgebra.h classes.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean_dat: 
	rm *.dat
clean_o: 
	rm *.o 