CXX      = /usr/lib64/openmpi/bin/mpic++
CXXFLAGS = -Wall -g

dft2d:	dft2d.o Complex.o InputImage.o
	$(CXX) -g -o dft2d dft2d.o Complex.o InputImage.o

clean:
	@rm *.o dft2d
