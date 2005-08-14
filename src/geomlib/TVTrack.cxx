//*************************************************************************
//* ===============
//*  TVTrack Class
//* ===============
//*
//* (Description)
//*   A class to implement a track object.
//* (Requires)
//*     TVCurve
//* (Provides)
//*     class TVTrack
//* (Update Recored)
//*   2003/10/24  K.Fujii       Original version.
//*
//*************************************************************************
//

#include <iostream>
#include "TVTrack.h"

using namespace std;

//_____________________________________________________________________
//  -----------------------------------
//  Track Class
//  -----------------------------------
//
//_____________________________________________________________________
//  --------------
//  Ctors and Dtor
//  --------------
//

TVTrack::TVTrack(Double_t dr,
                 Double_t phi0,
                 Double_t kappa,
                 Double_t dz,
                 Double_t tanl,
                 Double_t x0,
                 Double_t y0,
                 Double_t z0,
                 Double_t b)
             : fDrho(dr), fPhi0(phi0), fKappa(kappa), fDz(dz), fTanL(tanl),
               fX0(x0,y0,z0)
{
   SetMagField(b);
}

TVTrack::TVTrack(const TMatrixD &a, const TVector3 & x0, Double_t b)
             : fDrho(a(0,0)), fPhi0(a(1,0)), fKappa(a(2,0)), 
               fDz(a(3,0)), fTanL(a(4,0)),
               fX0(x0)
{
   SetMagField(b);
}