#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define bigsize (65536*16) // 16*1024*65536
/* Uses asynchronous I/O. Each process writes to separate files and
   reads them back. The file name is taken as a command-line argument,
   and the process rank is appended to it.*/

int main(int argc, char *argv[]) {
	int size;
	int lastsize;
	int i, rank, computeprocs, totalprocs, ndble, len, ndbleperproc, ndblelastproc;
	double *buf;
	char *_coordX, *_coordY, *_coordZ;
	char *tmp;
	int errs=0;
	MPI_File fh;
	MPI_Status status;
	MPI_Request request;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &totalprocs);
/* Ajustement des tailles du problème */
	/* totalprocs = computeprocs + master */
	computeprocs = totalprocs - 1;

/* Allocation des structures de données  */
	/* _coordX _coordY _coordZ */
	if(rank == 0) printf("bigsize is = %i bytes => %i doubles\n",bigsize,bigsize/sizeof(double));
	MPI_Barrier(MPI_COMM_WORLD);
	size = bigsize / computeprocs;

	len = 8;
	_coordX = (char *)malloc(len + 10);
	strcpy( _coordX, "coordX" );
	_coordY = (char *)malloc(len + 10);
	strcpy( _coordY, "coordY" );
	_coordZ = (char *)malloc(len + 10);
	strcpy( _coordZ, "coordZ" );

	if((rank!=computeprocs-1) || (rank!=0)) {
		buf = (double *)malloc(size);
	}
	else if (rank==computeprocs-1) {
		buf = (double *)malloc(bigsize - (computeprocs-1)*size);
	}
	ndble = size/sizeof(double);
	ndbleperproc = ndble / computeprocs;
	ndblelastproc= ndble - ndbleperproc * (computeprocs - 1);
	printf("size = %i and ndble = %i ndbleperproc = %i ndblelastproc = %i \n",size,ndble,ndbleperproc,ndblelastproc);
	if(rank!=0) {
		for (i=0; i<ndble; i++) {
			buf[i] = rank*100000 + sqrt(i);
		}
	}
/* Remplissage des structures de donnees */
	/* Un fichier _coordX.'myrank' par processus (hors master) */
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank!=0) {
		tmp = (char *) malloc(len+10);
		strcpy(tmp, _coordX);
		/* ------- */
		sprintf(_coordX, "%s.%d", tmp, rank);
		MPI_File_open(MPI_COMM_SELF, _coordX, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
		MPI_File_set_view(fh, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_iwrite(fh, buf, ndble, MPI_DOUBLE, &request);
//		MPI_Wait( &request, &status );
		MPI_File_close(&fh);
		/* ------- */
		strcpy(tmp, _coordY);
		sprintf(_coordY, "%s.%d", tmp, rank);
		MPI_File_open(MPI_COMM_SELF, _coordY, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
		MPI_File_set_view(fh, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_iwrite(fh, buf, ndble, MPI_DOUBLE, &request);
//		MPI_Wait( &request, &status );
		MPI_File_close(&fh);
		/* ------- */
		strcpy(tmp, _coordZ);
		sprintf(_coordZ, "%s.%d", tmp, rank);
		MPI_File_open(MPI_COMM_SELF, _coordZ, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
		MPI_File_set_view(fh, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_iwrite(fh, buf, ndble, MPI_DOUBLE, &request);
//		MPI_Wait( &request, &status );
		MPI_File_close(&fh);

	}
/* Boucle iterative */
	/* reopen the file and read the data back */
	for (i=0; i<ndble; i++) {
		buf[i] = 0;
	}
	if(rank!=0) {
		MPI_File_open(MPI_COMM_SELF, _coordX, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
		MPI_File_set_view(fh, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_iread(fh, buf, ndble, MPI_DOUBLE, &request);
//		MPI_Wait( &request, &status );
		MPI_File_close(&fh);
	}
	/* check if the data read is correct */
	if(rank != 0) {
		for (i=0; i<ndble; i++) {
			if ( buf[i] != (rank*100000 + sqrt(i)) ) {
				errs++;
				fprintf(stderr, "Process %d: error, read %d, should be %d\n", rank, buf[i], rank*100000+i);fflush(stderr);
			}
		}
		free(buf);
		free(_coordX);
		free(tmp);
	}
	MPI_Finalize();
	return errs;
}
