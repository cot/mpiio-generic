CC=mpicc

all: iread_double.c
	$(CC) iread_double.c -o iread_double

clean:
	rm -f iread_double coord*
