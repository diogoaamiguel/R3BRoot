#include "R3BCalifaPoint.h"

#include <iostream>

using std::cout;
using std::endl;
using std::flush;

R3BCalifaPoint::R3BCalifaPoint()
    : FairMCPoint()
    ,
    /*fX_out(0.), fY_out(0.), fZ_out(0.),
    fPx_out(0.), fPy_out(0.), fPz_out(0.),*/
    fCrystalId(-1)
    , fNf(0.)
    , fNs(0.)
{
}

R3BCalifaPoint::R3BCalifaPoint(Int_t trackID,
                               Int_t detID,
                               Int_t ident,
                               TVector3 posIn,
                               TVector3 momIn,
                               Double_t tof,
                               Double_t length,
                               Double_t eLoss,
                               Double_t Nf,
                               Double_t Ns,
                               UInt_t EventId)
    : FairMCPoint(trackID, detID, posIn, momIn, tof, length, eLoss, EventId)
    , fCrystalId(ident)
    , fNf(Nf)
    , fNs(Nf)
{
}

R3BCalifaPoint::R3BCalifaPoint(const R3BCalifaPoint& right)
    : FairMCPoint(right)
    , fCrystalId(right.fCrystalId)
    , fNf(right.fNf)
    , fNs(right.fNs)
{
}

R3BCalifaPoint::~R3BCalifaPoint() {}

void R3BCalifaPoint::Print(const Option_t* opt) const
{
    cout << "-I- R3BCalifaPoint: CALIFA Point for track " << fTrackID << " in crystal " << fCrystalId << endl;
    cout << "    Position (" << fX << ", " << fY << ", " << fZ << ") cm" << endl;
    cout << "    Momentum (" << fPx << ", " << fPy << ", " << fPz << ") GeV" << endl;
    cout << "    Time " << fTime << " ns,  Length " << fLength << " cm,  Energy loss " << fELoss * 1.0e06 << " keV"
         << endl;
    cout << "Nf " << fNf << ", Ns " << fNs << endl;
}

ClassImp(R3BCalifaPoint)
