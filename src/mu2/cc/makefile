gslinc = $(shell gsl-config --cflags)
gsllib = $(shell gsl-config --libs)

libs: libkcd.so

libkcd.so : scatter.c
	gcc ./scatter.c -Ofast -fPIC -shared $(gslinc) $(gsllib) -o $(HOME)/libkcd.so
