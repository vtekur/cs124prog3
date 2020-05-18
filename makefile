all: partition

randmst: partition.cpp
	g++ -O2 -Wall -g -o  partition partition.cpp

clean:
	rm -f partition