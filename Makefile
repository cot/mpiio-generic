CC=mpicc

check: CC+= -DCHECK
check: debug

debug: CC+= -g -pg
debug: all

all: iread_double.c
	$(CC) iread_double.c -o iread_double

clean:
	rm -f iread_double coord*
