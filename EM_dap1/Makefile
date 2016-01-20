CPPFLAGS= -O3 -I/wsu/home/groups/piquelab/include -fopenmp
LDFLAGS= -L/wsu/apps/gnu-4.7.2/python/Anaconda/Anaconda-2.1.0/lib -L/wsu/home/groups/piquelab/lib -L/wsu/home/groups/piquelab/lib64 -lm -lgsl -lgslcblas -lboost_iostreams -lz

main: main.o controller.o logistic.o
	g++ $(CPPFLAGS) $^ $(LDFLAGS) -o em_dap1

%.o: %.cc
	g++ -c $(CPPFLAGS) $^

clean:
	rm *.o em_dap1	



# main: main.o controller.o logistic.o
# 	g++  -O3 main.o controller.o logistic.o -lm -lgsl -lgslcblas -lboost_iostreams -lz -o em_dap1 
# main.o: main.cc 
# 	g++ -c  main.cc	
# controller.o: controller.cc classdef.h
# 	g++    -c  controller.cc
# logistic.o: logistic.cc logistic.h
# 	g++ -c logistic.cc 
# clean:
# 	rm *.o em_dap1	
