#include <iostream>
#include <fstream>

void output(const float* ro, const float* pr,
            const float* vx, const float* vy,
            float dx, float dy, int OutNum, int ix, int jx)
{
    // FileOpen
    static int num = 0;
    std::ofstream fp;
    if(num){
        fp.open("bOutput", std::ios::app | std::ios::binary);
    }
    else{
        fp.open("bOutput", std::ios::out | std::ios::binary);
    }
    if( fp.fail() ){
        throw 1;
    }

    // FileWrite
    if(!num){
        fp.write( (char*)&OutNum, sizeof(int) );
        fp.write( (char*)&ix, sizeof(int) );
        fp.write( (char*)&jx, sizeof(int) );
        for(int j=0; j<jx; j++){
            for(int i=0; i<ix; i++){
                const float x = dx*(i-1);
                fp.write( (char*)&x, sizeof(float) );
            }
        }
    
        for(int j=0; j<jx; j++){
            for(int i=0; i<ix; i++){
                const float y = dx*(j-1);
                fp.write( (char*)&y, sizeof(float) );
            }
        }
    }
    const int ixjx = ix*jx;
    fp.write( (char*)ro, sizeof(float)*ixjx );
    fp.write( (char*)pr, sizeof(float)*ixjx );
    fp.write( (char*)vx, sizeof(float)*ixjx );
    fp.write( (char*)vy, sizeof(float)*ixjx );

    num++;
    return;
}
