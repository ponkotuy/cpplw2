#include <omp.h>

#define MAT(i,j,ix) ((j)*(ix) + (i))

void bdfrex(float* u, int ix, int jx){
#pragma omp parallel for schedule(static)
    for(int j=0; j<jx; j++){
        u[MAT(0,j,ix)] = u[MAT(1,j,ix)];
        u[MAT(ix-1,j,ix)] = u[MAT(ix-2,j,ix)];
    }
    return;
}

void bdfrey(float* u, int ix, int jx){
#pragma omp parallel for schedule(static)
    for(int i=0; i<ix; i++){
        u[MAT(i,0,ix)] = u[MAT(i,1,ix)];
        u[MAT(i,jx-1,ix)] = u[MAT(i,jx-2,ix)];
    }
    return;
}
