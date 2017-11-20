main: main.o controller.o parser.o MLR.o
	g++  -fopenmp  -O3 main.o parser.o controller.o MLR.o -lm -L /usr/local/lib -lgsl -lgslcblas -o dap-g
static:	main.o controller.o parser.o MLR.o 
	g++  -fopenmp  -O3 main.o parser.o controller.o M.o -lm -L /usr/local/lib -lgsl -lgslcblas -static -o  dap-g.static
main.o: main.cc
	g++ -c -O3  main.cc
parser.o: parser.cc parser.h
	g++ -c -O3 parser.cc
controller.o: controller.cc controller.h
	g++ -c -O3 -fopenmp  controller.cc
MLR.o: MLR.cc MLR.h
	g++ -c -O3 MLR.cc
clean:
	rm *.o dap-g
