main: main.o controllerZ.o  
	g++   -fopenmp  -O3  main.o  controllerZ.o -lm  -o dap_ss
main.o: main.cc
	g++ -c -O3  main.cc
controllerZ.o: controllerZ.cc controllerZ.h
	g++ -c   -O3    controllerZ.cc
clean:
	rm *.o dap_ss		
