galsim: galsim.c
	gcc -w -O3 -o galsim galsim.c -lm -ffast-math -fopenmp

clean:
	rm -f galsim