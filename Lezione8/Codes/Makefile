CC = g++
CFLAGS = -Wall -O3 --std=c++11

makeall: main.exe main2.exe

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
main.exe : main.o random.o WaveFunctionBase.h
	$(CC) random.o main.o -o main.exe
main.o : main.cpp main.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
main2.exe : main2.o random.o WaveFunctionBase.h vectoralgebra.h
	$(CC) random.o main2.o -o main2.exe
main2.o: main2.cpp main2.h
	$(CC) -c main2.cpp -o main2.o $(CFLAGS)
clean_dat: 
	rm *.dat
clean_o: 
	rm *.o 