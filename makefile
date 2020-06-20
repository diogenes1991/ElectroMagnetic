
CXX = g++-8
CUU = nvcc

OBJECTS=$(patsubst %.cu,%.o,$(wildcard *.cu))


all: compile run

compile: $(OBJECTS)
	$(CUU) $(OBJECTS) -o Test

run:
	@./Test

%.o:%.cu
	$(CUU) -c $< -o $@
%.o:%.cc
	$(CXX) -c $< -c $@

