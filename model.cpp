
void bdfrex(float* u, int ix, int jx);
void bdfrey(float* u, int ix, int jx);

void bdperx(float* u, int ix, int jx);
void bdpery(float* u, int ix, int jx);

void DefaultValue(float* ro, float* pr, float* vx, float* vy,
                  int ix, int jx)
{
    int id = 0;
    for(int j=0; j<jx; j++){
        const float ro0 = 1.0, ro1 = 0.1;
        const float pr0 = 1.0, pr1 = 0.05;
        // const float ro0 = 1.0, ro1 = 0.1;
        // const float pr0 = 1.0, pr1 = pr0;
        const float vx0 = 0.0, vx1 = vx0;
        const float vy0 = 0.0, vy1 = vy0;
        // Left Side
        for(int i=0; i<(ix/2); i++, id++){
            ro[id] = ro0;
            pr[id] = pr0;
            vx[id] = vx0;
            vy[id] = vy0;
        }
        // Right Side
        for(int i=(ix/2); i<ix; i++, id++){
            ro[id] = ro1;
            pr[id] = pr1;
            vx[id] = vx1;
            vy[id] = vy1;
        }
    }
    return;
}

void bnd(float* ro, float* pr, float* vx, float* vy, int ix, int jx){
    // x-direction
    bdfrex(ro, ix, jx);
    bdfrex(pr, ix, jx);
    bdfrex(vx, ix, jx);
    bdfrex(vy, ix, jx);

    // y-direction
    bdpery(ro, ix, jx);
    bdpery(pr, ix, jx);
    bdpery(vx, ix, jx);
    bdpery(vy, ix, jx);
    return;
}
