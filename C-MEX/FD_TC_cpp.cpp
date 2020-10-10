// FD code
// version: 0.1
// Ivan Abakumov, Lei Li.
// Publication date: 10th October 2019

//
//////////////////////////////////////////////////////
//c++  libraries
#include "mex.h"
#include <math.h>
#include <algorithm>
#include <assert.h>  
#include <fstream>
#include <cstdlib>
#include <matrix.h>
#include <string.h>
#include <complex.h> 
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <numeric>
#include <iterator>
#include <vector>
#include <random>
#include <iostream>
#include <chrono>
#include <functional>


using namespace std;

#define info 0

class GridClass{
    public:
        int    nx;
        int    ny;
        int    nz;
        int    nt;
        double x0;
        double dx;
        double y0;
        double dy;
        double z0;
        double dz;
        double t0;
        double dt;
        double mx;
        double my;
        double mz;
        double mt;
        double xx;
        double yy;
        double zz;
        double tt;
};

typedef std::vector<double> Vector1D;
typedef std::vector<Vector1D> Vector2D;
typedef std::vector<Vector2D> Vector3D;

//double fdc[] = {1225.0/1024.0, -245.0/3072.0, 49.0/5120.0, -5.0/7168.0};
//int NN = 4;
double fdc[] = {1.211, -0.090, 0.014, -0.002, 0.0001};
int NN = 5;

///////////////////////////////////////////////////////////////////////////
///////////////        STOP WATCH                 /////////////////////////
///////////////////////////////////////////////////////////////////////////

class StopWatch {
public:
   typedef std::chrono::steady_clock clock;
   typedef std::chrono::microseconds microseconds;
   typedef std::chrono::milliseconds milliseconds;
   typedef std::chrono::seconds seconds;

   StopWatch();
   StopWatch(const StopWatch&);
   StopWatch& operator=(const StopWatch& rhs);

   uint64_t ElapsedUs() const;
   uint64_t ElapsedMs() const;
   uint64_t ElapsedSec() const;

   std::chrono::steady_clock::time_point Restart();

private:
   clock::time_point mStart;
};

StopWatch::StopWatch() : mStart(clock::now()) {
   static_assert(std::chrono::steady_clock::is_steady, "Serious OS/C++ library issues. Steady clock is not steady");
   // FYI:  This would fail  static_assert(std::chrono::high_resolution_clock::is_steady(), "High Resolution Clock is NOT steady on CentOS?!");
}

StopWatch::StopWatch(const StopWatch& other): mStart(other.mStart) { 
}

/// @return StopWatch::StopWatch&  - assignment operator.
StopWatch& StopWatch::operator=(const StopWatch& rhs) {
      mStart = rhs.mStart;
      return *this;
}

/// @return the elapsed microseconds since start
uint64_t StopWatch::ElapsedUs() const {
   return std::chrono::duration_cast<microseconds>(clock::now() - mStart).count();
}

/// @return the elapsed milliseconds since start
uint64_t StopWatch::ElapsedMs() const {
   return std::chrono::duration_cast<milliseconds>(clock::now() - mStart).count();
}

/// @return the elapsed seconds since start
uint64_t StopWatch::ElapsedSec() const {
   return std::chrono::duration_cast<seconds>(clock::now() - mStart).count();
}
/**
 * Resets the start point
 * @return the updated start point
 */
std::chrono::steady_clock::time_point StopWatch::Restart() {
   mStart = clock::now();
   return mStart;
}

///////////////////////////////////////////////////////////////////////////
///////////////        AUXILIARY ROUTINES         /////////////////////////
///////////////////////////////////////////////////////////////////////////


// Forward Derivatives

double Dxfm(Vector3D &Q, int i, int j, int k)
{
    double Dxf = fdc[0]*(Q[i+1][j][k] - Q[i  ][j][k]) +
                 fdc[1]*(Q[i+2][j][k] - Q[i-1][j][k]) +
                 fdc[2]*(Q[i+3][j][k] - Q[i-2][j][k]) +
                 fdc[3]*(Q[i+4][j][k] - Q[i-3][j][k]) +
                 fdc[4]*(Q[i+5][j][k] - Q[i-4][j][k]);
    return Dxf;
}

double Dyfm(Vector3D &Q, int i, int j, int k)
{
    double Dyf = fdc[0]*(Q[i][j+1][k] - Q[i][j  ][k]) +
                 fdc[1]*(Q[i][j+2][k] - Q[i][j-1][k]) +
                 fdc[2]*(Q[i][j+3][k] - Q[i][j-2][k]) +
                 fdc[3]*(Q[i][j+4][k] - Q[i][j-3][k]) +
                 fdc[4]*(Q[i][j+5][k] - Q[i][j-4][k]);
    return Dyf;
}

double Dzfm(Vector3D &Q, int i, int j, int k)
{
    double Dzf = fdc[0]*(Q[i][j][k+1] - Q[i][j][k  ]) +
                 fdc[1]*(Q[i][j][k+2] - Q[i][j][k-1]) +
                 fdc[2]*(Q[i][j][k+3] - Q[i][j][k-2]) +
                 fdc[3]*(Q[i][j][k+4] - Q[i][j][k-3]) +
                 fdc[4]*(Q[i][j][k+5] - Q[i][j][k-4]);
    return Dzf;
}

// Backward Derivatives

double Dxbm(Vector3D &Q, int i, int j, int k)
{
    double Dxb = fdc[0]*(Q[i  ][j][k] - Q[i-1][j][k]) +
                 fdc[1]*(Q[i+1][j][k] - Q[i-2][j][k]) +
                 fdc[2]*(Q[i+2][j][k] - Q[i-3][j][k]) +
                 fdc[3]*(Q[i+3][j][k] - Q[i-4][j][k]) +
                 fdc[4]*(Q[i+4][j][k] - Q[i-5][j][k]);
    return Dxb;
}

double Dybm(Vector3D &Q, int i, int j, int k)
{
    double Dyb = fdc[0]*(Q[i][j  ][k] - Q[i][j-1][k]) +
                 fdc[1]*(Q[i][j+1][k] - Q[i][j-2][k]) +
                 fdc[2]*(Q[i][j+2][k] - Q[i][j-3][k]) +
                 fdc[3]*(Q[i][j+3][k] - Q[i][j-4][k]) +
                 fdc[4]*(Q[i][j+4][k] - Q[i][j-5][k]);
    return Dyb;
}

double Dzbm(Vector3D &Q, int i, int j, int k)
{
    double Dzb = fdc[0]*(Q[i][j][k  ] - Q[i][j][k-1]) +
                 fdc[1]*(Q[i][j][k+1] - Q[i][j][k-2]) +
                 fdc[2]*(Q[i][j][k+2] - Q[i][j][k-3]) +
                 fdc[3]*(Q[i][j][k+3] - Q[i][j][k-4]) +
                 fdc[4]*(Q[i][j][k+4] - Q[i][j][k-5]);
    return Dzb;
}

///////////////////////////////////////////////////////////////////////////
///////////////           MAIN   CODE             /////////////////////////
///////////////////////////////////////////////////////////////////////////

void FD (GridClass G, double *c11, double *c22, double *c33, double *c44, double *c55, double *c66, double *c12, double *c13, double *c23,
                      double *bx,  double *by,  double *bz,  double *rx,  double *ry,  double *rz,  double *sx,  double *sy,  double *sz,
                      int Nrece, int Nshot, double *pmlx,double *pmly,double *pmlz,double *mt,  double *ft,  double *SX, double *SY, double *SZ) 
{
   
    StopWatch sw;
    
    int i,j,k,t,r;
    int gsx, gsy, gsz; 
    
    gsx = (int) (sx[0]-1);
    gsy = (int) (sy[0]-1); 
    gsz = (int) (sz[0]-1); 
    
    int grx[Nrece], gry[Nrece], grz[Nrece]; 
    
    for (r=0;r<Nrece;r++){
        grx[r] = (int) (rx[r]-1);
        gry[r] = (int) (ry[r]-1); 
        grz[r] = (int) (rz[r]-1); 
    }
    
    double MT[3][3];
    
    for (i=0;i<3;i++){
         for (j=0;j<3;j++){
             MT[i][j] = mt[i+j*3];
         }
    }
    
    double FT[G.nt];
    
    for (t=0;t<G.nt;t++){
        FT[t] = ft[t];
    }
    
    double dTxxdxf, dTxydyb, dTxzdzb;
    double dTxydxb, dTyydyf, dTyzdzb;
    double dTxzdxb, dTyzdyb, dTzzdzf;
    double dVxdxb, dVydyb, dVzdzb; 
    double dVydxf, dVxdyf, dVzdxf; 
    double dVxdzf, dVzdyf, dVydzf;
    
    Vector3D C11(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D C22(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D C33(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D C44(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D C55(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D C66(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D C12(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D C13(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D C23(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    
    Vector3D Vx  (G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Vx_x(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Vx_y(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Vx_z(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Vy  (G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Vy_x(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Vy_y(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Vy_z(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Vz  (G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Vz_x(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Vz_y(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Vz_z(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    
    Vector3D Txx  (G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Txx_x(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Txx_y(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Txx_z(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Tyy  (G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Tyy_x(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Tyy_y(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Tyy_z(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Tzz  (G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Tzz_x(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Tzz_y(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Tzz_z(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    
    Vector3D Txy  (G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Txy_x(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Txy_y(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Txz  (G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Txz_x(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Txz_z(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Tyz  (G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Tyz_y(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Tyz_z(G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    
    Vector3D Pxp  (G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Pyp  (G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Pzp  (G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Pxm  (G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Pym  (G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D Pzm  (G.nx,Vector2D(G.ny,Vector1D(G.nz)));
    
    Vector3D Bx  (G.nx-1,Vector2D(G.ny,Vector1D(G.nz)));
    Vector3D By  (G.nx,Vector2D(G.ny-1,Vector1D(G.nz)));
    Vector3D Bz  (G.nx,Vector2D(G.ny,Vector1D(G.nz-1)));
    
    // Assign values
    for (i=0;i<G.nx;i++){
        for (j=0;j<G.ny;j++){
            for (k=0;k<G.nz;k++){
                C11[i][j][k] = c11[i+j*G.nx+k*G.nx*G.ny];
                C22[i][j][k] = c22[i+j*G.nx+k*G.nx*G.ny];
                C33[i][j][k] = c33[i+j*G.nx+k*G.nx*G.ny];
                C44[i][j][k] = c44[i+j*G.nx+k*G.nx*G.ny];
                C55[i][j][k] = c55[i+j*G.nx+k*G.nx*G.ny];
                C66[i][j][k] = c66[i+j*G.nx+k*G.nx*G.ny];
                C12[i][j][k] = c12[i+j*G.nx+k*G.nx*G.ny];
                C13[i][j][k] = c13[i+j*G.nx+k*G.nx*G.ny];
                C23[i][j][k] = c23[i+j*G.nx+k*G.nx*G.ny];
                
                Vx[i][j][k] = 0.0;
                Vx_x[i][j][k] = 0.0;
                Vx_y[i][j][k] = 0.0;
                Vx_z[i][j][k] = 0.0;
                Vy[i][j][k] = 0.0;
                Vy_x[i][j][k] = 0.0;
                Vy_y[i][j][k] = 0.0;
                Vy_z[i][j][k] = 0.0;
                Vz[i][j][k] = 0.0;
                Vz_x[i][j][k] = 0.0;
                Vz_y[i][j][k] = 0.0;
                Vz_z[i][j][k] = 0.0;
                
                Txx[i][j][k] = 0.0;
                Txx_x[i][j][k] = 0.0;
                Txx_y[i][j][k] = 0.0;
                Txx_z[i][j][k] = 0.0;
                Tyy[i][j][k] = 0.0;
                Tyy_x[i][j][k] = 0.0;
                Tyy_y[i][j][k] = 0.0;
                Tyy_z[i][j][k] = 0.0;
                Tzz[i][j][k] = 0.0;
                Tzz_x[i][j][k] = 0.0;
                Tzz_y[i][j][k] = 0.0;
                Tzz_z[i][j][k] = 0.0;
                
                Txy[i][j][k] = 0.0;
                Txy_x[i][j][k] = 0.0;
                Txy_y[i][j][k] = 0.0;
                Tyz[i][j][k] = 0.0;
                Tyz_y[i][j][k] = 0.0;
                Tyz_z[i][j][k] = 0.0;
                Txz[i][j][k] = 0.0;
                Txz_x[i][j][k] = 0.0;
                Txz_z[i][j][k] = 0.0;
                
                Pxp[i][j][k] = 1.0 + 0.5*G.dt*pmlx[i+j*G.nx+k*G.nx*G.ny]; 
                Pyp[i][j][k] = 1.0 + 0.5*G.dt*pmly[i+j*G.nx+k*G.nx*G.ny]; 
                Pzp[i][j][k] = 1.0 + 0.5*G.dt*pmlz[i+j*G.nx+k*G.nx*G.ny]; 
                Pxm[i][j][k] = 1.0 - 0.5*G.dt*pmlx[i+j*G.nx+k*G.nx*G.ny]; 
                Pym[i][j][k] = 1.0 - 0.5*G.dt*pmly[i+j*G.nx+k*G.nx*G.ny]; 
                Pzm[i][j][k] = 1.0 - 0.5*G.dt*pmlz[i+j*G.nx+k*G.nx*G.ny]; 
            }
        }
    }
   
    for (i=0;i<G.nx-1;i++){
        for (j=0;j<G.ny;j++){
            for (k=0;k<G.nz;k++){
                Bx[i][j][k] = bx[i+j*(G.nx-1)+k*(G.nx-1)*G.ny];
            }
        }
    }
    for (i=0;i<G.nx;i++){
        for (j=0;j<G.ny-1;j++){
            for (k=0;k<G.nz;k++){
                By[i][j][k] = by[i+j*G.nx+k*G.nx*(G.ny-1)];
            }
        }
    }
    for (i=0;i<G.nx;i++){
        for (j=0;j<G.ny;j++){
            for (k=0;k<G.nz-1;k++){
                Bz[i][j][k] = bz[i+j*G.nx+k*G.nx*G.ny];
            }
        }
    }
    
    for (t=0;t<G.nt;t++){

        Txx_x[gsx][gsy][gsz] = Txx_x[gsx][gsy][gsz] + (-MT[0][0]/3)*FT[t];
        Txx_y[gsx][gsy][gsz] = Txx_y[gsx][gsy][gsz] + (-MT[0][0]/3)*FT[t];
        Txx_z[gsx][gsy][gsz] = Txx_z[gsx][gsy][gsz] + (-MT[0][0]/3)*FT[t];
        Tyy_x[gsx][gsy][gsz] = Tyy_x[gsx][gsy][gsz] + (-MT[1][1]/3)*FT[t];
        Tyy_y[gsx][gsy][gsz] = Tyy_y[gsx][gsy][gsz] + (-MT[1][1]/3)*FT[t];
        Tyy_z[gsx][gsy][gsz] = Tyy_z[gsx][gsy][gsz] + (-MT[1][1]/3)*FT[t];
        Tzz_x[gsx][gsy][gsz] = Tzz_x[gsx][gsy][gsz] + (-MT[2][2]/3)*FT[t];
        Tzz_y[gsx][gsy][gsz] = Tzz_y[gsx][gsy][gsz] + (-MT[2][2]/3)*FT[t];
        Tzz_z[gsx][gsy][gsz] = Tzz_z[gsx][gsy][gsz] + (-MT[2][2]/3)*FT[t];
        Txy_x[gsx][gsy][gsz] = Txy_x[gsx][gsy][gsz] + (-MT[0][1]/2)*FT[t];
        Txy_y[gsx][gsy][gsz] = Txy_y[gsx][gsy][gsz] + (-MT[0][1]/2)*FT[t];
        Txz_x[gsx][gsy][gsz] = Txz_x[gsx][gsy][gsz] + (-MT[0][2]/2)*FT[t];
        Txz_z[gsx][gsy][gsz] = Txz_z[gsx][gsy][gsz] + (-MT[0][2]/2)*FT[t];
        Tyz_y[gsx][gsy][gsz] = Tyz_y[gsx][gsy][gsz] + (-MT[1][2]/2)*FT[t];
        Tyz_z[gsx][gsy][gsz] = Tyz_z[gsx][gsy][gsz] + (-MT[1][2]/2)*FT[t];
        
        for (i=0;i<G.nx;i++){
            for (j=0;j<G.ny;j++){
                for (k=0;k<G.nz;k++){
        
                    Txx[i][j][k] = Txx_x[i][j][k] + Txx_y[i][j][k] + Txx_z[i][j][k];
                    Tyy[i][j][k] = Tyy_x[i][j][k] + Tyy_y[i][j][k] + Tyy_z[i][j][k];
                    Tzz[i][j][k] = Tzz_x[i][j][k] + Tzz_y[i][j][k] + Tzz_z[i][j][k];

                    Txy[i][j][k] = Txy_x[i][j][k] + Txy_y[i][j][k];
                    Txz[i][j][k] = Txz_x[i][j][k] + Txz_z[i][j][k];
                    Tyz[i][j][k] = Tyz_y[i][j][k] + Tyz_z[i][j][k];     
                }
            }
        }
        
        for (i=NN;i<(G.nx-NN);i++){
            for (j=NN;j<(G.ny-NN);j++){
                for (k=NN;k<(G.nz-NN);k++){
                    
                    /* 
                    dTxxdxf = Dxfm(Txx,i,j,k);
                    dTxxdxf = Dxfm(Txx,i,j,k);
                    dTxydyb = Dybm(Txy,i,j,k);
                    dTxzdzb = Dzbm(Txz,i,j,k);
                    dTxydxb = Dxbm(Txy,i,j,k);
                    dTyydyf = Dyfm(Tyy,i,j,k);
                    dTyzdzb = Dzbm(Tyz,i,j,k);
                    dTxzdxb = Dxbm(Txz,i,j,k);
                    dTyzdyb = Dybm(Tyz,i,j,k);
                    dTzzdzf = Dzfm(Tzz,i,j,k);
                    */
                    
                    dTxxdxf = fdc[0]*(Txx[i+1][j][k] - Txx[i  ][j][k]) +
                              fdc[1]*(Txx[i+2][j][k] - Txx[i-1][j][k]) +
                              fdc[2]*(Txx[i+3][j][k] - Txx[i-2][j][k]) +
                              fdc[3]*(Txx[i+4][j][k] - Txx[i-3][j][k]) +
                              fdc[4]*(Txx[i+5][j][k] - Txx[i-4][j][k]);
                    
                    dTxydyb = fdc[0]*(Txy[i][j  ][k] - Txy[i][j-1][k]) +
                              fdc[1]*(Txy[i][j+1][k] - Txy[i][j-2][k]) +
                              fdc[2]*(Txy[i][j+2][k] - Txy[i][j-3][k]) +
                              fdc[3]*(Txy[i][j+3][k] - Txy[i][j-4][k]) +
                              fdc[4]*(Txy[i][j+4][k] - Txy[i][j-5][k]);
                    
                    dTxzdzb = fdc[0]*(Txz[i][j][k  ] - Txz[i][j][k-1]) +
                              fdc[1]*(Txz[i][j][k+1] - Txz[i][j][k-2]) +
                              fdc[2]*(Txz[i][j][k+2] - Txz[i][j][k-3]) +
                              fdc[3]*(Txz[i][j][k+3] - Txz[i][j][k-4]) +
                              fdc[4]*(Txz[i][j][k+4] - Txz[i][j][k-5]);
                    
                    dTxydxb = fdc[0]*(Txy[i  ][j][k] - Txy[i-1][j][k]) +
                              fdc[1]*(Txy[i+1][j][k] - Txy[i-2][j][k]) +
                              fdc[2]*(Txy[i+2][j][k] - Txy[i-3][j][k]) +
                              fdc[3]*(Txy[i+3][j][k] - Txy[i-4][j][k]) +
                              fdc[4]*(Txy[i+4][j][k] - Txy[i-5][j][k]);
                    
                    dTyydyf = fdc[0]*(Tyy[i][j+1][k] - Tyy[i][j  ][k]) +
                              fdc[1]*(Tyy[i][j+2][k] - Tyy[i][j-1][k]) +
                              fdc[2]*(Tyy[i][j+3][k] - Tyy[i][j-2][k]) +
                              fdc[3]*(Tyy[i][j+4][k] - Tyy[i][j-3][k]) +
                              fdc[4]*(Tyy[i][j+5][k] - Tyy[i][j-4][k]);
                    
                    dTyzdzb = fdc[0]*(Tyz[i][j][k  ] - Tyz[i][j][k-1]) +
                              fdc[1]*(Tyz[i][j][k+1] - Tyz[i][j][k-2]) +
                              fdc[2]*(Tyz[i][j][k+2] - Tyz[i][j][k-3]) +
                              fdc[3]*(Tyz[i][j][k+3] - Tyz[i][j][k-4]) +
                              fdc[4]*(Tyz[i][j][k+4] - Tyz[i][j][k-5]);
                    
                    dTxzdxb = fdc[0]*(Txz[i  ][j][k] - Txz[i-1][j][k]) +
                              fdc[1]*(Txz[i+1][j][k] - Txz[i-2][j][k]) +
                              fdc[2]*(Txz[i+2][j][k] - Txz[i-3][j][k]) +
                              fdc[3]*(Txz[i+3][j][k] - Txz[i-4][j][k]) +
                              fdc[4]*(Txz[i+4][j][k] - Txz[i-5][j][k]);
                    
                    dTyzdyb = fdc[0]*(Tyz[i][j  ][k] - Tyz[i][j-1][k]) +
                              fdc[1]*(Tyz[i][j+1][k] - Tyz[i][j-2][k]) +
                              fdc[2]*(Tyz[i][j+2][k] - Tyz[i][j-3][k]) +
                              fdc[3]*(Tyz[i][j+3][k] - Tyz[i][j-4][k]) +
                              fdc[4]*(Tyz[i][j+4][k] - Tyz[i][j-5][k]);
                    
                    dTzzdzf = fdc[0]*(Tzz[i][j][k+1] - Tzz[i][j][k  ]) +
                              fdc[1]*(Tzz[i][j][k+2] - Tzz[i][j][k-1]) +
                              fdc[2]*(Tzz[i][j][k+3] - Tzz[i][j][k-2]) +
                              fdc[3]*(Tzz[i][j][k+4] - Tzz[i][j][k-3]) +
                              fdc[4]*(Tzz[i][j][k+5] - Tzz[i][j][k-4]);
                    
                    
                    Vx_x[i][j][k] = (Pxm[i][j][k]*Vx_x[i][j][k] + G.dt*Bx[i][j][k]*dTxxdxf/G.dx)/Pxp[i][j][k];  
                    Vx_y[i][j][k] = (Pym[i][j][k]*Vx_y[i][j][k] + G.dt*Bx[i][j][k]*dTxydyb/G.dy)/Pyp[i][j][k];
                    Vx_z[i][j][k] = (Pzm[i][j][k]*Vx_z[i][j][k] + G.dt*Bx[i][j][k]*dTxzdzb/G.dz)/Pzp[i][j][k];

                    Vy_x[i][j][k] = (Pxm[i][j][k]*Vy_x[i][j][k] + G.dt*By[i][j][k]*dTxydxb/G.dx)/Pxp[i][j][k];
                    Vy_y[i][j][k] = (Pym[i][j][k]*Vy_y[i][j][k] + G.dt*By[i][j][k]*dTyydyf/G.dy)/Pyp[i][j][k];
                    Vy_z[i][j][k] = (Pzm[i][j][k]*Vy_z[i][j][k] + G.dt*By[i][j][k]*dTyzdzb/G.dz)/Pzp[i][j][k];
   
                    Vz_x[i][j][k] = (Pxm[i][j][k]*Vz_x[i][j][k] + G.dt*Bz[i][j][k]*dTxzdxb/G.dx)/Pxp[i][j][k];
                    Vz_y[i][j][k] = (Pym[i][j][k]*Vz_y[i][j][k] + G.dt*Bz[i][j][k]*dTyzdyb/G.dy)/Pyp[i][j][k];
                    Vz_z[i][j][k] = (Pzm[i][j][k]*Vz_z[i][j][k] + G.dt*Bz[i][j][k]*dTzzdzf/G.dz)/Pzp[i][j][k];
                    
                    Vx[i][j][k] = Vx_x[i][j][k] + Vx_y[i][j][k] + Vx_z[i][j][k];
                    Vy[i][j][k] = Vy_x[i][j][k] + Vy_y[i][j][k] + Vy_z[i][j][k];
                    Vz[i][j][k] = Vz_x[i][j][k] + Vz_y[i][j][k] + Vz_z[i][j][k];
                    
                }
            }
        } 
        
        /*    % topFS with the assumption of weak anisotropy near the surface
        if strcmp(BCtype,'topFS')
            Vz(h-1,ii, jj) = Vz(h,ii, jj) +... 
                (lam(h,ii, jj)./lamu(h,ii, jj)).*(...
                Vx(h,ii, jj)  -  Vx(h,(ii)-1, jj)...
                +Vx(h,ii, jj)  -  Vx(h,ii, (jj)-1)   ); 
            %Vx
            %Vy
             Vx(h-1,ii, jj) = Vx(h,ii, jj) ... 
                +Vz(h-1,(ii)+1, jj)  -  Vz(h-1,ii, jj)...
                +Vz(h  ,(ii)+1, jj)  -  Vz(h  ,ii, jj)...
                +Vx(h+1,ii, jj)  -  Vx(h,ii, jj)  ;
            
            Vy(h-1,ii, jj) = Vy(h,ii, jj) ... 
                +Vz(h-1,(ii)+1, jj)  -  Vz(h-1,ii, jj)...
                +Vz(h  ,(ii)+1, jj)  -  Vz(h  ,ii, jj)...
                +Vx(h+1,ii, jj)  -  Vx(h,ii, jj)  ;
                
            Vz(h-2,ii, jj) = Vz(h-1,ii, jj) +... 
                (lam(h,ii, jj)./lamu(h,ii, jj)).*(...
                 Vx(h-1,ii, jj)  -  Vx(h-1,(ii)-1, jj)...
                +Vy(h-1,ii, jj)  -  Vy(h-1,ii, (jj)-1)   ); 
           
        end
         **/
        
        for (i=NN;i<(G.nx-NN);i++){
            for (j=NN;j<(G.ny-NN);j++){
                for (k=NN;k<(G.nz-NN);k++){
        
                    /*
                    dVxdxb = Dxbm(Vx,i,j,k);
                    dVydyb = Dybm(Vy,i,j,k);
                    dVzdzb = Dzbm(Vz,i,j,k);
                    dVydxf = Dxfm(Vy,i,j,k);
                    dVxdyf = Dyfm(Vx,i,j,k);
                    dVzdxf = Dxfm(Vz,i,j,k);
                    dVxdzf = Dzfm(Vx,i,j,k);
                    dVzdyf = Dyfm(Vz,i,j,k);
                    dVydzf = Dzfm(Vy,i,j,k);
                    */
                    
                    dVxdxb = fdc[0]*(Vx[i  ][j][k] - Vx[i-1][j][k]) +
                             fdc[1]*(Vx[i+1][j][k] - Vx[i-2][j][k]) +
                             fdc[2]*(Vx[i+2][j][k] - Vx[i-3][j][k]) +
                             fdc[3]*(Vx[i+3][j][k] - Vx[i-4][j][k]) +
                             fdc[4]*(Vx[i+4][j][k] - Vx[i-5][j][k]);
                    
                    dVydyb = fdc[0]*(Vy[i][j  ][k] - Vy[i][j-1][k]) +
                             fdc[1]*(Vy[i][j+1][k] - Vy[i][j-2][k]) +
                             fdc[2]*(Vy[i][j+2][k] - Vy[i][j-3][k]) +
                             fdc[3]*(Vy[i][j+3][k] - Vy[i][j-4][k]) +
                             fdc[4]*(Vy[i][j+4][k] - Vy[i][j-5][k]);
                        
                    dVzdzb = fdc[0]*(Vz[i][j][k  ] - Vz[i][j][k-1]) +
                             fdc[1]*(Vz[i][j][k+1] - Vz[i][j][k-2]) +
                             fdc[2]*(Vz[i][j][k+2] - Vz[i][j][k-3]) +
                             fdc[3]*(Vz[i][j][k+3] - Vz[i][j][k-4]) +
                             fdc[4]*(Vz[i][j][k+4] - Vz[i][j][k-5]);
                    
                    dVydxf = fdc[0]*(Vy[i+1][j][k] - Vy[i  ][j][k]) +
                             fdc[1]*(Vy[i+2][j][k] - Vy[i-1][j][k]) +
                             fdc[2]*(Vy[i+3][j][k] - Vy[i-2][j][k]) +
                             fdc[3]*(Vy[i+4][j][k] - Vy[i-3][j][k]) +
                             fdc[4]*(Vy[i+5][j][k] - Vy[i-4][j][k]);
                        
                    dVxdyf = fdc[0]*(Vx[i][j+1][k] - Vx[i][j  ][k]) +
                             fdc[1]*(Vx[i][j+2][k] - Vx[i][j-1][k]) +
                             fdc[2]*(Vx[i][j+3][k] - Vx[i][j-2][k]) +
                             fdc[3]*(Vx[i][j+4][k] - Vx[i][j-3][k]) +
                             fdc[4]*(Vx[i][j+5][k] - Vx[i][j-4][k]);
                        
                    dVzdxf = fdc[0]*(Vz[i+1][j][k] - Vz[i  ][j][k]) +
                             fdc[1]*(Vz[i+2][j][k] - Vz[i-1][j][k]) +
                             fdc[2]*(Vz[i+3][j][k] - Vz[i-2][j][k]) +
                             fdc[3]*(Vz[i+4][j][k] - Vz[i-3][j][k]) +
                             fdc[4]*(Vz[i+5][j][k] - Vz[i-4][j][k]);
                        
                    dVxdzf = fdc[0]*(Vx[i][j][k+1] - Vx[i][j][k  ]) +
                             fdc[1]*(Vx[i][j][k+2] - Vx[i][j][k-1]) +
                             fdc[2]*(Vx[i][j][k+3] - Vx[i][j][k-2]) +
                             fdc[3]*(Vx[i][j][k+4] - Vx[i][j][k-3]) +
                             fdc[4]*(Vx[i][j][k+5] - Vx[i][j][k-4]);
                        
                    dVzdyf = fdc[0]*(Vz[i][j+1][k] - Vz[i][j  ][k]) +
                             fdc[1]*(Vz[i][j+2][k] - Vz[i][j-1][k]) +
                             fdc[2]*(Vz[i][j+3][k] - Vz[i][j-2][k]) +
                             fdc[3]*(Vz[i][j+4][k] - Vz[i][j-3][k]) +
                             fdc[4]*(Vz[i][j+5][k] - Vz[i][j-4][k]);
                        
                    dVydzf = fdc[0]*(Vy[i][j][k+1] - Vy[i][j][k  ]) +
                             fdc[1]*(Vy[i][j][k+2] - Vy[i][j][k-1]) +
                             fdc[2]*(Vy[i][j][k+3] - Vy[i][j][k-2]) +
                             fdc[3]*(Vy[i][j][k+4] - Vy[i][j][k-3]) +
                             fdc[4]*(Vy[i][j][k+5] - Vy[i][j][k-4]);
                    
                    Txx_x[i][j][k] = (Pxm[i][j][k]*Txx_x[i][j][k] + G.dt*C11[i][j][k]*dVxdxb/G.dx)/Pxp[i][j][k];
                    Txx_y[i][j][k] = (Pym[i][j][k]*Txx_y[i][j][k] + G.dt*C12[i][j][k]*dVydyb/G.dy)/Pyp[i][j][k];
                    Txx_z[i][j][k] = (Pzm[i][j][k]*Txx_z[i][j][k] + G.dt*C13[i][j][k]*dVzdzb/G.dz)/Pzp[i][j][k];

                    Tyy_x[i][j][k] = (Pxm[i][j][k]*Tyy_x[i][j][k] + G.dt*C12[i][j][k]*dVxdxb/G.dx)/Pxp[i][j][k];
                    Tyy_y[i][j][k] = (Pym[i][j][k]*Tyy_y[i][j][k] + G.dt*C22[i][j][k]*dVydyb/G.dy)/Pyp[i][j][k];
                    Tyy_z[i][j][k] = (Pzm[i][j][k]*Tyy_z[i][j][k] + G.dt*C23[i][j][k]*dVzdzb/G.dz)/Pzp[i][j][k];

                    Tzz_x[i][j][k] = (Pxm[i][j][k]*Tzz_x[i][j][k] + G.dt*C13[i][j][k]*dVxdxb/G.dx)/Pxp[i][j][k];
                    Tzz_y[i][j][k] = (Pym[i][j][k]*Tzz_y[i][j][k] + G.dt*C23[i][j][k]*dVydyb/G.dy)/Pyp[i][j][k];
                    Tzz_z[i][j][k] = (Pzm[i][j][k]*Tzz_z[i][j][k] + G.dt*C33[i][j][k]*dVzdzb/G.dz)/Pzp[i][j][k];

                    Txy_x[i][j][k] = (Pxm[i][j][k]*Txy_x[i][j][k] + G.dt*C66[i][j][k]*dVydxf/G.dx)/Pxp[i][j][k];
                    Txy_y[i][j][k] = (Pym[i][j][k]*Txy_y[i][j][k] + G.dt*C66[i][j][k]*dVxdyf/G.dy)/Pyp[i][j][k];
                    Txz_x[i][j][k] = (Pxm[i][j][k]*Txz_x[i][j][k] + G.dt*C55[i][j][k]*dVzdxf/G.dx)/Pxp[i][j][k];
                    Txz_z[i][j][k] = (Pzm[i][j][k]*Txz_z[i][j][k] + G.dt*C55[i][j][k]*dVxdzf/G.dz)/Pzp[i][j][k];
                    Tyz_y[i][j][k] = (Pym[i][j][k]*Tyz_y[i][j][k] + G.dt*C44[i][j][k]*dVzdyf/G.dy)/Pyp[i][j][k];
                    Tyz_z[i][j][k] = (Pzm[i][j][k]*Tyz_z[i][j][k] + G.dt*C44[i][j][k]*dVydzf/G.dz)/Pzp[i][j][k];

                    Txx[i][j][k] = Txx_x[i][j][k] + Txx_y[i][j][k] + Txx_z[i][j][k];
                    Tyy[i][j][k] = Tyy_x[i][j][k] + Tyy_y[i][j][k] + Tyy_z[i][j][k];
                    Tzz[i][j][k] = Tzz_x[i][j][k] + Tzz_y[i][j][k] + Tzz_z[i][j][k];
                    Txy[i][j][k] = Txy_x[i][j][k] + Txy_y[i][j][k];
                    Txz[i][j][k] = Txz_x[i][j][k] + Txz_z[i][j][k];
                    Tyz[i][j][k] = Tyz_y[i][j][k] + Tyz_z[i][j][k]; 
                }
            }
        }
        
        /*    % topFS with the assumption of weak anisotropy near the surface
            %if strcmpi(BCtype,'topFS')
            %        Tzz(h,:,:)=0;
            %        Tzz(h-1,:,:)=-Tzz(h+1,:,:);
            %        Txz(h-1,:,:)=-Txz(h,:,:);
            %        Txz(h-2,:,:)=-Txz(h+1,:,:);
            %        Tyz(h-1,:,:)=-Tyz(h,:,:);
            %        Tyz(h-2,:,:)=-Tyz(h+1,:,:);
            %end
         */ 
        //sw.Restart();
        //cout << "CPU time: " << sw.ElapsedUs() << endl;

        for (r=0;r<Nrece;r++){
            SX[t+r*G.nt] = Vx[grx[r]][gry[r]][grz[r]];
            SY[t+r*G.nt] = Vy[grx[r]][gry[r]][grz[r]];
            SZ[t+r*G.nt] = Vz[grx[r]][gry[r]][grz[r]];
        }
    }
    return; 
}
    




///////////////////////////////////////////////////////////////////////////
//////////////////      Gateway routine      //////////////////////////////
///////////////////////////////////////////////////////////////////////////
// Transfer Data and results from MATLAB to C++ and back
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    
    const mwSize *dims;
    mwSize ndim;

    mwSize mrows, ncols; 
    int Nrece, Nshot;
    double *tmp;
    double *C11, *C22, *C33, *C44, *C55, *C66, *C12, *C13, *C23, *bx, *by, *bz;
    double *rx, *ry, *rz, *sx, *sy, *sz; 
    double *pmlx, *pmly, *pmlz; 
    double *MT, *FT;
    
    GridClass G;
    
    
    double *SX, *SY, *SZ;
     
     //check size of input and output arguments
    if (nrhs!=5) mexErrMsgTxt("5 input arguments");  
    if (nlhs!=3) mexErrMsgTxt("10 output argument");    
    
    //////////////////////////////////////////////////////////////////////
     
    // Get G - parameters of the velocity grid
    if (!mxIsClass(prhs[0], "GridClass")) mexErrMsgTxt("Input (1st arg.) must be an object of class GridClass");
    if (mxGetProperty(prhs[0],0,"nx")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'nx' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"nx"));
	G.nx = int(*tmp);
    if (mxGetProperty(prhs[0],0,"ny")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'ny' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"ny"));
	G.ny = int(*tmp);
    if (mxGetProperty(prhs[0],0,"nz")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'nz' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"nz"));
	G.nz = int(*tmp);
    if (mxGetProperty(prhs[0],0,"nt")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'nt' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"nt"));
	G.nt = int(*tmp);
    if (mxGetProperty(prhs[0],0,"x0")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'x0' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"x0"));
	G.x0 = double(*tmp);
    if (mxGetProperty(prhs[0],0,"y0")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'y0' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"y0"));
	G.y0 = double(*tmp);
    if (mxGetProperty(prhs[0],0,"z0")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'z0' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"z0"));
	G.z0 = double(*tmp);
    if (mxGetProperty(prhs[0],0,"t0")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 't0' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"t0"));
	G.t0 = double(*tmp);
    if (mxGetProperty(prhs[0],0,"dx")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'dx' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"dx"));
	G.dx = double(*tmp);
    if (mxGetProperty(prhs[0],0,"dy")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'dy' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"dy"));
	G.dy = double(*tmp);
    if (mxGetProperty(prhs[0],0,"dz")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'dz' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"dz"));
	G.dz = double(*tmp);
    if (mxGetProperty(prhs[0],0,"dt")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'dt' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"dt"));
	G.dt = double(*tmp);
    if (mxGetProperty(prhs[0],0,"mx")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'mx' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"mx"));
	G.mx = double(*tmp);
    if (mxGetProperty(prhs[0],0,"my")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'my' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"my"));
	G.my = double(*tmp);
    if (mxGetProperty(prhs[0],0,"mz")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'mz' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"mz"));
	G.mz = double(*tmp);
    if (mxGetProperty(prhs[0],0,"mt")==NULL)   mexErrMsgTxt("Input (1st arg): Required Property 'mt' is missing.");
	tmp = (double *)mxGetPr(mxGetProperty(prhs[0],0,"mt"));
	G.mt = double(*tmp);

    // Get MODEL 
    if (!(mxIsStruct(prhs[1]))) mexErrMsgTxt("Input (2d arg.) must be a 'Model' structure.");
    if (mxGetField(prhs[1],0,"C11")==NULL)   mexErrMsgTxt("Input (2d arg): Required field 'C11' is missing.");
    C11 = (double *)mxGetPr(mxGetField(prhs[1],0,"C11"));
    if (mxGetField(prhs[1],0,"C22")==NULL)   mexErrMsgTxt("Input (2d arg): Required field 'C22' is missing.");
    C22 = (double *)mxGetPr(mxGetField(prhs[1],0,"C22"));
    if (mxGetField(prhs[1],0,"C33")==NULL)   mexErrMsgTxt("Input (2d arg): Required field 'C33' is missing.");
    C33 = (double *)mxGetPr(mxGetField(prhs[1],0,"C33"));
    if (mxGetField(prhs[1],0,"C44")==NULL)   mexErrMsgTxt("Input (2d arg): Required field 'C44' is missing.");
    C44 = (double *)mxGetPr(mxGetField(prhs[1],0,"C44"));
    if (mxGetField(prhs[1],0,"C55")==NULL)   mexErrMsgTxt("Input (2d arg): Required field 'C55' is missing.");
    C55 = (double *)mxGetPr(mxGetField(prhs[1],0,"C55"));
    if (mxGetField(prhs[1],0,"C66")==NULL)   mexErrMsgTxt("Input (2d arg): Required field 'C66' is missing.");
    C66 = (double *)mxGetPr(mxGetField(prhs[1],0,"C66"));
    if (mxGetField(prhs[1],0,"C12")==NULL)   mexErrMsgTxt("Input (2d arg): Required field 'C12' is missing.");
    C12 = (double *)mxGetPr(mxGetField(prhs[1],0,"C12"));
    if (mxGetField(prhs[1],0,"C13")==NULL)   mexErrMsgTxt("Input (2d arg): Required field 'C13' is missing.");
    C13 = (double *)mxGetPr(mxGetField(prhs[1],0,"C13"));
    if (mxGetField(prhs[1],0,"C23")==NULL)   mexErrMsgTxt("Input (2d arg): Required field 'C23' is missing.");
    C23 = (double *)mxGetPr(mxGetField(prhs[1],0,"C23"));
    if (mxGetField(prhs[1],0,"bx")==NULL)   mexErrMsgTxt("Input (2d arg): Required field 'bx' is missing.");
    bx = (double *)mxGetPr(mxGetField(prhs[1],0,"bx"));    
    if (mxGetField(prhs[1],0,"by")==NULL)   mexErrMsgTxt("Input (2d arg): Required field 'by' is missing.");
    by = (double *)mxGetPr(mxGetField(prhs[1],0,"by"));    
    if (mxGetField(prhs[1],0,"bz")==NULL)   mexErrMsgTxt("Input (2d arg): Required field 'bz' is missing.");
    bz = (double *)mxGetPr(mxGetField(prhs[1],0,"bz"));
    
    // Get ACQ
    if (!(mxIsStruct(prhs[2]))) mexErrMsgTxt("Input (3rd arg.) must be a structure.");
    if (mxGetField(prhs[2],0,"rx")==NULL)   mexErrMsgTxt("Input (3rd arg): Required field 'rx' is missing.");
    if (mxGetField(prhs[2],0,"ry")==NULL)   mexErrMsgTxt("Input (3rd arg): Required field 'ry' is missing.");
    if (mxGetField(prhs[2],0,"rz")==NULL)   mexErrMsgTxt("Input (3rd arg): Required field 'rz' is missing.");
    if (mxGetField(prhs[2],0,"sx")==NULL)   mexErrMsgTxt("Input (3rd arg): Required field 'sx' is missing.");
    if (mxGetField(prhs[2],0,"sy")==NULL)   mexErrMsgTxt("Input (3rd arg): Required field 'sy' is missing.");
    if (mxGetField(prhs[2],0,"sz")==NULL)   mexErrMsgTxt("Input (3rd arg): Required field 'sz' is missing.");
    if (mxGetField(prhs[2],0,"nrece")==NULL)   mexErrMsgTxt("Input (3rd arg): Required field 'nrece' is missing.");
    if (mxGetField(prhs[2],0,"nshot")==NULL)   mexErrMsgTxt("Input (3rd arg): Required field 'nshot' is missing.");
    rx = (double *)mxGetPr(mxGetField(prhs[2],0,"rx"));    
    ry = (double *)mxGetPr(mxGetField(prhs[2],0,"ry"));    
    rz = (double *)mxGetPr(mxGetField(prhs[2],0,"rz"));    
    sx = (double *)mxGetPr(mxGetField(prhs[2],0,"sx"));    
    sy = (double *)mxGetPr(mxGetField(prhs[2],0,"sy"));    
    sz = (double *)mxGetPr(mxGetField(prhs[2],0,"sz"));  
    tmp = (double *)mxGetPr(mxGetField(prhs[2],0,"nrece")); 
    Nrece = int(*tmp);
    tmp = (double *)mxGetPr(mxGetField(prhs[2],0,"nshot"));
    Nshot = int(*tmp);
    // Get BC
    if (!(mxIsStruct(prhs[3]))) mexErrMsgTxt("Input (4th arg.) must be a structure.");
    if (mxGetField(prhs[3],0,"pmlx")==NULL)   mexErrMsgTxt("Input (4th arg): Required field 'pmlx' is missing.");
    if (mxGetField(prhs[3],0,"pmly")==NULL)   mexErrMsgTxt("Input (4th arg): Required field 'pmly' is missing.");
    if (mxGetField(prhs[3],0,"pmlz")==NULL)   mexErrMsgTxt("Input (4th arg): Required field 'pmlz' is missing.");
    pmlx = (double *)mxGetPr(mxGetField(prhs[3],0,"pmlx"));    
    pmly = (double *)mxGetPr(mxGetField(prhs[3],0,"pmly"));    
    pmlz = (double *)mxGetPr(mxGetField(prhs[3],0,"pmlz"));    
     
    // Get Source
    if (!(mxIsStruct(prhs[4]))) mexErrMsgTxt("Input (5th arg.) must be a structure.");
    if (mxGetField(prhs[4],0,"MT")==NULL)   mexErrMsgTxt("Input (5th arg): Required field 'MT' is missing.");
    if (mxGetField(prhs[4],0,"FT")==NULL)   mexErrMsgTxt("Input (5th arg): Required field 'FT' is missing.");
    MT = (double *)mxGetPr(mxGetField(prhs[4],0,"MT"));    
    FT = (double *)mxGetPr(mxGetField(prhs[4],0,"FT"));
     
    // OUTPUT
    mrows = (mwSize)G.nt;
    ncols = (mwSize)Nrece;
     
    plhs[0] = mxCreateDoubleMatrix(mrows,ncols, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(mrows,ncols, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(mrows,ncols, mxREAL);
    
    //associate output pointer
    SX = (double *)mxGetPr(plhs[0]);   
    SY = (double *)mxGetPr(plhs[1]);   
    SZ = (double *)mxGetPr(plhs[2]);   
     
    //////////////////////////////////////////////////////////////////////
    
    if(info!= 0) mexPrintf("Calling routine from Mex gateway...\n");
    FD(G,C11,C22,C33,C44,C55,C66,C12,C13,C23,bx,by,bz,rx,ry,rz,sx,sy,sz,Nrece,Nshot,pmlx,pmly,pmlz,MT,FT,SX,SY,SZ);
    if(info!= 0) mexPrintf("FD finished.\n");
     
    //////////////////////////////////////////////////////////////////////
    
    return;
}
