// -*- coding:utf-8 -*-

#include <iostream>
#include <cmath>
#include <omp.h>

#define MAT(i,j,ix) ((j)*(ix) + (i))

inline void cal_u(float * __restrict__ u, float ro, float pr, float vx, float vy, float gamma);
inline void cal_fg(float * __restrict__ f, float * __restrict__ g,
            float ro, float pr, float vx, float vy, float gamma);
inline void u_to_fg(float * __restrict__ f, float * __restrict__ g, const float * __restrict__ u, float gamma);

void lw2(float * __restrict__ ro, float * __restrict__ pr, float * __restrict__ vx, float * __restrict__ vy,
         float dt, float dx, float dy, float qv, float gamma,
         int ix, int jx)
{
    const float dxi = 1.0/dx;
    const float dyi = 1.0/dy;
    const int ixjx = ix*jx;
    float* f = new float[ixjx*4];
    float* g = new float[ixjx*4];
    float* u = new float[ixjx*4];
    float* u_ = new float[ixjx*4];

    float* qx = new float[ixjx];
    float* qy = new float[ixjx];

#pragma omp parallel
    {

        // Initialize 11+23
#pragma omp for schedule(guided)
        for(int id=0; id<ixjx; id++){
            cal_u(&u[id*4], ro[id], pr[id], vx[id], vy[id], gamma);
            cal_fg(&f[id*4], &g[id*4],
                   ro[id], pr[id], vx[id], vy[id], gamma);
        }

        // First Step 18*4
#pragma omp for schedule(guided)
        for(int j=0; j<jx-1; ++j){
            for(int i=0; i<ix-1; ++i){
                const int id = MAT(i,j,ix);
                for(int n=0; n<4; n++){
                    u_[id*4+n] = 0.25*( u[id*4+n] + u[(id+1)*4+n]
                                        + u[(id+ix)*4+n] + u[(id+ix+1)*4+n])
                        - 0.25*dt*dxi*( f[(id+1)*4+n] - f[id*4+n]
                                        + f[(id+ix+1)*4+n] - f[(id+ix)*4+n])
                        - 0.25*dt*dyi*( g[(id+ix)*4+n] - g[id*4+n]
                                        + g[(id+ix+1)*4+n] - g[(id+1)*4+n]);
                }
            }
        }

        // u -> F,G 25
#pragma omp for schedule(guided)
        for(int j=0; j<jx-1; ++j){
            for(int i=0; i<ix-1; ++i){
                const int id = MAT(i,j,ix);
                u_to_fg(&f[id*4], &g[id*4], &u_[id*4], gamma);
            }
        }

        // Second Step 14*4
#pragma omp for schedule(guided)
        for(int j=1; j<jx-1; ++j){
            for(int i=1; i<ix-1; ++i){
                const int id = MAT(i,j,ix);
                for(int n=0; n<4; n++){
                    u_[id*4+n] = u[id*4+n]
                        - 0.5*dt*dxi*( f[id*4+n] - f[(id-1)*4+n]
                                       + f[(id-ix)*4+n] - f[(id-ix-1)*4+n] )
                        - 0.5*dt*dyi*( g[id*4+n] - g[(id-ix)*4+n]
                                       + g[(id-1)*4+n] - g[(id-ix-1)*4+n]);
                }
            }
        }

        // Termination 4 + 4 + 16*4+10
#pragma omp for schedule(static)
        for(int j=0; j<jx; ++j){
            for(int i=0; i<ix-1; ++i){
                const int id = MAT(i,j,ix);
                qx[id] = qv * fmaxf(0.0, fabs(vx[id+1] - vx[id]) - 1.0e-4);
            }
        }
#pragma omp for schedule(static)
        for(int j=0; j<jx-1; ++j){
            for(int i=0; i<ix; ++i){
                const int id = MAT(i,j,ix);
                qy[id] = qv * fmaxf(0.0, fabs(vy[id+ix] - vy[id]) - 1.0e-4);
            }
        }
#pragma omp for schedule(guided)
        for(int j=1; j<jx-1; ++j){
            for(int i=1; i<ix-1; ++i){
                const int id = MAT(i,j,ix);
                for(int n=0; n<4; n++){
                    u_[id*4+n] += dt * dxi
                        * ( qx[id-1] * ( u[(id-1)*4+n] - u[id*4+n] )
                            + qx[id] * ( u[(id+1)*4+n] - u[id*4+n] ) );
                    u_[id*4+n] += dt * dyi
                        * ( qy[id-ix] * ( u[(id-ix)*4+n] - u[id*4+n] )
                            + qy[id] * ( u[(id+ix)*4+n] - u[id*4+n] ) );
                }
                const float ro_ = u_[id*4];
                const float vx_ = u_[id*4+1]/ro_;
                const float vy_ = u_[id*4+2]/ro_;
                const float v2 = vx_*vx_ + vy_*vy_;
                const float e_ = u_[id*4+3]; // e_= Energy*Density
                ro[id] = ro_;
                pr[id] = (gamma - 1.0) * (e_ - 0.5*v2*ro_);
                vx[id] = vx_;
                vy[id] = vy_;
            }
        }
    }
    
    delete[] f; delete[] g; delete[] u; delete[] u_;
    delete[] qx; delete[] qy;
    return;
}

// 物理量->U 11
inline void cal_u(float * __restrict__ u, float ro, float pr, float vx, float vy, float gamma){
    const float v2 = vx*vx + vy*vy;
    const float e = 0.5*v2 + pr / ((gamma-1.0)*ro);
    u[0] = ro;
    u[1] = ro*vx;
    u[2] = ro*vy;
    u[3] = ro*e;
    return;
}

// 物理量->F,G 23
inline void cal_fg(float* __restrict__ f, float* __restrict__ g,
            float ro, float pr, float vx, float vy,
            float gamma)
{
    const float v2 = vx*vx + vy*vy;
    const float h = 0.5*v2 + gamma*pr / ((gamma-1.0)*ro);
    f[0] = ro*vx;
    f[1] = ro*vx*vx + pr;
    f[2] = ro*vx*vy;
    f[3] = ro*h*vx;
    g[0] = ro*vy;
    g[1] = f[2];
    g[2] = ro*vy*vy + pr;
    g[3] = ro*h*vy;
    return;
}

// U -> F,G 25
inline void u_to_fg(float* __restrict__ f, float* __restrict__ g, const float* __restrict__ u, float gamma){
    const float vx = u[1] / u[0];
    const float vy = u[2] / u[0];
    const float rovx2 = vx*u[1]; // ro*vx*vx
    const float rovy2 = vy*u[2]; // ro*vy*vy
    const float ke = 0.5*(rovx2 + rovy2); // Kinetic Energy
    f[0] = u[1];
    f[1] = (gamma - 1.0) * (u[3] - ke) + rovx2;
    f[2] = vx * u[2];
    f[3] = vx * (gamma*u[3] - (gamma-1.0)*ke);
    g[0] = u[2];
    g[1] = f[2];
    g[2] = (gamma - 1.0) * (u[3] - ke) + rovy2;
    g[3] = vy * (gamma*u[3] - (gamma-1.0)*ke);
    return;
}
