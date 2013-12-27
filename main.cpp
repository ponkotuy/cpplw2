// -*- coding:utf-8 -*-

// Calc : 278*(ix-1)*(jx-1)*step

#include <iostream>
#include <cmath>
#include <cfloat>
#include <omp.h>
#define MAT(i,j,ix) ((j)*(ix) + (i))

void DefaultValue(float* ro, float* pr, float* vx, float* vy,
                  int ix, int jx);
float cfl(const float* ro, const float* pr,
          const float* vx, const float* vy,
          float dx, float gamma, float safety, int ix, int jx);
void lw2(float* ro, float* pr, float* vx, float* vy,
         float dt, float dx, float dy, float qav, float gamma,
         int ix, int jx);
void output(const float* ro, const float* pr,
            const float* vx, const float* vy,
            float dx, float dy, int OutNum, int ix, int jx);
void bnd(float* ro, float* pr, float* vx, float* vy, int ix, int jx);
void chkdav(float* u, float* vx, float* vy, float floor,
            int ix, int jx);

int main(){
    const int ix = 1024;        // Grid Size X
    const int jx = 1024;         // Grid Size Y
    const float tMax = 0.2;	// Time Step
    const float OutStep = 0.05;	// Output per step
    const int OutNum = (int)(tMax/OutStep) + 1;
    const float dx = 0.004;
    const float dy = 0.004;
    const float gamma = 5.0/3.0;
    const float floor = 1.0e-9;	// Threshold Density&Pressure
    const float safety = 0.4;	// Courant Number Safety
    const float qv = 3.0;	// Artificial Viscosity
    
    float* ro = new float[ix*jx]; // Density
    float* pr = new float[ix*jx]; // Pressure
    float* vx = new float[ix*jx]; // Velocity x
    float* vy = new float[ix*jx]; // Velocity y

    DefaultValue(ro, pr, vx, vy, ix, jx);

    // Main
    float t = 0.0;
    int step = 0;
    while(t<tMax){
	
        // Output
        static float OutTime = 0.0;
        if(OutTime<=t){
            std::cout << "t = " << t << std::endl;
            output(ro, pr, vx, vy, dx, dy, OutNum, ix, jx);
            OutTime += OutStep;
        }
	
        // determine dt 9
        float dt = cfl(ro, pr, vx, vy,
                       fminf(dx,dy),
                       gamma, safety, ix, jx);
        if(dt<1.0e-9) return 1;

        // Lax-Wendroff 269
        lw2(ro, pr, vx, vy, dt, dx, dy, qv, gamma, ix, jx);

        // Boundary Condition
        bnd(ro, pr, vx, vy, ix, jx);

        // Check ro&pr
        chkdav(ro, vx, vy, floor, ix, jx);
        chkdav(pr, vx, vy, floor, ix, jx);

        // Next Step
        t += dt;
        ++step;
    }
    std::cout << "t = " << t << std::endl;
    output(ro, pr, vx, vy, dx, dy, OutNum, ix, jx);
    std::cout << "step = " << step << std::endl;
    
    delete[] ro;
    delete[] pr;
    delete[] vx;
    delete[] vy;
    return 0;
}

float cfl(const float* ro, const float* pr,
          const float* vx, const float* vy,
          float dx, float gamma, float safety, int ix, int jx)
{
    float* dtq = new float[ix*jx];
    float dt = FLT_MAX;

#pragma omp parallel for schedule(guided)
    for(int j=1; j<jx-1; j++){
        for(int i=1; i<ix-1; i++){
            const int id = MAT(i,j,ix);
            const float v2 = vx[id]*vx[id] + vy[id]*vy[id];
            dtq[id] = safety*dx / sqrtf(v2 + gamma*pr[id]/ro[id]);
        }
    }

    for(int j=1; j<jx-1; j++){
        for(int i=1; i<ix-1; i++){
            dt = fminf( dt, dtq[MAT(i,j,ix)] );
        }
    }
    return dt;
}

void chkdav(float* u, float* vx, float* vy, float floor,
            int ix, int jx)
{
    const int ixjx = ix*jx;
#pragma omp parallel for schedule(static)
    for(int id=0; id<ixjx; id++){
        if(u[id]<floor){
            u[id] = floor;
            vx[id] = 0.0;
            vy[id] = 0.0;
        }
    }
    return;
}
