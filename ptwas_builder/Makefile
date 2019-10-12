main: main.o parser.o 
	g++   -O3 main.o parser.o -lm -L /usr/local/lib -lgsl -lgslcblas -o ptwas_builder
static:	main.o parser.o
	g++   -O3 main.o parser.o -lm -L /usr/local/lib -lgsl -lgslcblas -static -o ptwas_builder.static
main.o: main.cc
	g++ -c -O3  main.cc
parser.o: parser.cc parser.h
	g++ -c -O3 parser.cc
clean:
	rm *.o ptwas_builder
