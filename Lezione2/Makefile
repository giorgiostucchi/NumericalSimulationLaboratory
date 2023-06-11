CC = g++
CFLAGS = -Wall -O3 --std=c++11

makeall: Esercizio1.exe Esercizio2.exe Esercizio3.exe

Esercizio1.exe : Esercizio1.o random.o
	$(CC) random.o Esercizio1.o -o Esercizio1.exe
Esercizio1.o : Esercizio1.cpp
	$(CC) -c Esercizio1.cpp -o Esercizio1.o $(CFLAGS)
Esercizio2.exe : Esercizio2.o random.o
	$(CC) random.o Esercizio2.o -o Esercizio2.exe
Esercizio2.o : Esercizio2.cpp
	$(CC) -c Esercizio2.cpp -o Esercizio2.o $(CFLAGS)
Esercizio3.exe : Esercizio3.o random.o
	$(CC) random.o Esercizio3.o -o Esercizio3.exe
Esercizio3.o : Esercizio3.cpp
	$(CC) -c Esercizio3.cpp -o Esercizio3.o $(CFLAGS)
main.exe : main.o random.o
	$(CC) random.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe *.txt
