lab2: run

run: comp
	./lab2

comp: lab2.cpp Matrix.cpp Matrix.h
	g++ lab2.cpp Matrix.cpp -o lab2

