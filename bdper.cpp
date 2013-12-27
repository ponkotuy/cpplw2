
#ifndef PERIODICAL_BOUNDARY_CONDITION
#define PERIODICAL_BOUNDARY_CONDITION

#include <omp.h>

#define MAT(i,j,ix) ((j)*(ix) + (i))

void bdperx(float* u, int ix, int jx){
#pragma omp parallel for schedule(static)
	for(int j=0; j<jx; j++){
        u[MAT(0,j,ix)] = u[MAT(ix-2,j,ix)];
        u[MAT(ix-1,j,ix)] = u[MAT(1,j,ix)];
	}
	return;
}

void bdpery(float* u, int ix, int jx){
#pragma omp parallel for schedule(static)
	for(int i=0; i<ix; i++){
        u[MAT(i,0,ix)] = u[MAT(i,jx-2,ix)];
        u[MAT(i,jx-1,ix)] = u[MAT(i,1,ix)];
	}
	return;
}

#endif	/* #ifndef PERIODICAL_BOUNDARY_CONDITION */
