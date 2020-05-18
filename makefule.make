all: strassen

randmst: strassen.cpp
	g++ -O2 -Wall -g -o  strassen strassen.cpp

clean:
	rm -f strassen