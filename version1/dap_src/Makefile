main: main.o controller.o parser.o SSLR.o
	g++  -fopenmp  -O3 main.o parser.o controller.o SSLR.o -lm -L /usr/local/lib -lgsl -lgslcblas -o dap
static:	main.o controller.o parser.o SSLR.o 
	g++  -fopenmp  -O3 main.o parser.o controller.o SSLR.o -lm -L /usr/local/lib -lgsl -lgslcblas -static -o dap.static
main.o: main.cc
	g++ -c -O3  main.cc
parser.o: parser.cc parser.h
	g++ -c -O3 parser.cc
controller.o: controller.cc controller.h
	g++ -c -O3 -fopenmp  controller.cc
SSLR.o: SSLR.h SSLR.cc
	g++ -c -O3 SSLR.cc
clean:
	rm *.o dap
