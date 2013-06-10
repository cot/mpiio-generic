#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define bigsize (65536*16) // 16*1024*65536
#define iteration 1

int main(int argc, char *argv[]) {
	int size;
	int lastsize;
	int i,j,k,l,iter;
	int rank;
	int computeprocs, totalprocs, ndbleperproc, ndblelastproc;
	int ndble, len;

	double *bufX; /* buffer des positions */
	double *bufY; /* buffer des positions */
	double *bufZ; /* buffer des positions */
	double _partX, _partY, _partZ; /* positions de la particule etudiee  */
	double _fX, _fY, _fZ; /* champ de force */

	char *_coordX, *_coordY, *_coordZ;
	char *tmp;
	int errs=0;
	MPI_File fhX;
	MPI_File fhY;
	MPI_File fhZ;
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

	if(rank!=computeprocs-1) {
		bufX = (double *)malloc(size);
		bufY = (double *)malloc(size);
		bufZ = (double *)malloc(size);
	}
	else if (rank==computeprocs-1) {
		bufX = (double *)malloc(bigsize - (computeprocs-1)*size);
		bufY = (double *)malloc(bigsize - (computeprocs-1)*size);
		bufZ = (double *)malloc(bigsize - (computeprocs-1)*size);
	}
	ndble = size/sizeof(double);
	ndbleperproc = ndble / computeprocs;
	ndblelastproc= ndble - ndbleperproc * (computeprocs - 1);
	printf("size = %i and ndble = %i ndbleperproc = %i ndblelastproc = %i \n",size,ndble,ndbleperproc,ndblelastproc);
	if(rank!=0) {
		for (i=0; i<ndble; i++) {
			bufX[i] = rank*100000 + sqrt(i);
			bufY[i] = rank*100000 + sqrt(i);
			bufZ[i] = rank*100000 + sqrt(i);
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
		MPI_File_open(MPI_COMM_SELF, _coordX, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fhX);
		MPI_File_set_view(fhX, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_iwrite(fhX, bufX, ndble, MPI_DOUBLE, &request);
		MPI_Wait( &request, &status );
		MPI_File_close(&fhX);
		/* ------- */
		strcpy(tmp, _coordY);
		sprintf(_coordY, "%s.%d", tmp, rank);
		MPI_File_open(MPI_COMM_SELF, _coordY, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fhY);
		MPI_File_set_view(fhY, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_iwrite(fhY, bufY, ndble, MPI_DOUBLE, &request);
		MPI_Wait( &request, &status );
		MPI_File_close(&fhY);
		/* ------- */
		strcpy(tmp, _coordZ);
		sprintf(_coordZ, "%s.%d", tmp, rank);
		MPI_File_open(MPI_COMM_SELF, _coordZ, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fhZ);
		MPI_File_set_view(fhZ, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
		MPI_File_iwrite(fhZ, bufZ, ndble, MPI_DOUBLE, &request);
		MPI_Wait( &request, &status );
		MPI_File_close(&fhZ);

	}
/* Boucle iterative */
	for(iter=0;iter<iteration;iter++) {
		/* ouverture des fichiers de donnees pour acces aux positions */
		for (i=0; i<ndble; i++) {
			bufX[i] = 0;
			bufY[i] = 0;
			bufZ[i] = 0;
		}
		/* MAJ des structures de donnees sur les esclaves */
		if(rank!=0) {

			MPI_File_open(MPI_COMM_SELF, _coordX, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fhX);
			MPI_File_set_view(fhX, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
			MPI_File_iread(fhX, bufX, ndble, MPI_DOUBLE, &request);
			MPI_Wait( &request, &status );
			MPI_File_close(&fhX);

			MPI_File_open(MPI_COMM_SELF, _coordY, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fhY);
			MPI_File_set_view(fhY, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
			MPI_File_iread(fhY, bufY, ndble, MPI_DOUBLE, &request);
			MPI_Wait( &request, &status );
			MPI_File_close(&fhY);

			MPI_File_open(MPI_COMM_SELF, _coordZ, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fhZ);
			MPI_File_set_view(fhZ, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
			MPI_File_iread(fhZ, bufZ, ndble, MPI_DOUBLE, &request);
			MPI_Wait( &request, &status );
			MPI_File_close(&fhZ);

			_partX = 0.0;
			_partY = 0.0;
			_partZ = 0.0;
			_fX    = 0.0;
			_fY    = 0.0;
			_fZ    = 0.0;

		}

		MPI_Barrier(MPI_COMM_WORLD);
		for(i=1; i<computeprocs ; i++) {
			/* MAJ des structures de donnees sur le master, bloc par bloc */
			if(rank==0) {

				sprintf(_coordX, "coordX.%d",  i);
				sprintf(_coordY, "coordY.%d",  i);
				sprintf(_coordZ, "coordZ.%d",  i);

		                tmp = (char *) malloc(len+10);
				strcpy(tmp, _coordX);
//				printf("tmp = %s \n",tmp);

				MPI_File_open(MPI_COMM_SELF, _coordX, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fhX);
				MPI_File_set_view(fhX, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
				MPI_File_iread(fhX, bufX, ndble, MPI_DOUBLE, &request);
				MPI_Wait( &request, &status );
				MPI_File_close(&fhX);

				MPI_File_open(MPI_COMM_SELF, _coordY, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fhY);
				MPI_File_set_view(fhY, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
				MPI_File_iread(fhY, bufY, ndble, MPI_DOUBLE, &request);
				MPI_Wait( &request, &status );
				MPI_File_close(&fhY);

				MPI_File_open(MPI_COMM_SELF, _coordZ, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fhZ);
				MPI_File_set_view(fhZ, 0, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
				MPI_File_iread(fhZ, bufZ, ndble, MPI_DOUBLE, &request);
				MPI_Wait( &request, &status );
				MPI_File_close(&fhZ);
				

			}
		}

	


	}

/* Verification des donnees */
	if(rank != 0) {
		for (i=0; i<ndble; i++) {
			if ( bufX[i] != (rank*100000 + sqrt(i)) ) {
				errs++;
				fprintf(stderr, "Process %d: error, read %d, should be %d\n", rank, bufX[i], rank*100000+i);fflush(stderr);
			}
		}
		free(bufX);
		free(bufY);
		free(bufZ);
		free(_coordX);
		free(_coordY);
		free(_coordZ);
		free(tmp);
	}
	MPI_Finalize();
	return errs;
}
