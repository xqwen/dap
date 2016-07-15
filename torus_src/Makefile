main: main.o controller.o logistic.o 
	g++  -O3 main.o controller.o logistic.o  -lm -lgsl -lgslcblas -lboost_iostreams -lz  -o torus 
static: main.o controller.o logistic.o
	g++  -O3 main.o controller.o logistic.o  -lm -lgsl -lgslcblas -lboost_iostreams -lz -static -o torus.static
main.o: main.cc 
	g++ -c  main.cc	
controller.o: controller.cc classdef.h
	g++    -c  controller.cc 
logistic.o: logistic.cc logistic.h
	g++ -c logistic.cc -fpermissive
clean:
	rm *.o torus	
