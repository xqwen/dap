main: main.o controller.o  SSLR_SS.o
	g++   -fopenmp  -O3  main.o  controller.o SSLR_SS.o -lm -L /usr/local/lib -lgsl -lgslcblas -o dap_ss
main.o: main.cc
	g++ -c -O3  main.cc
controller.o: controller.cc controller.h
	g++ -c  -fopenmp  -O3    controller.cc
SSLR.o: SSLR_SS.h SSLR_SS.cc
	g++ -c -O3  SSLR_SS.cc
clean:
	rm *.o dap_ss		
