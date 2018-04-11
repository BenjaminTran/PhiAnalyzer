#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/kaon.h"

kaon::kaon(TVector3 momentum, double eta, double phi, int charge, double dedx, double mass, bool isGen)
{
    //TVector3 momentum(px,py,pz);
    double energy = sqrt(TMath::Power(mass,2) + TMath::Power(momentum.Mag(),2));
    TLorentzVector PtEtaPhiE(momentum.Perp(),eta,phi,energy);
    momentum_ = momentum;
    PtEtaPhiE_ = PtEtaPhiE;
    dedx_ = dedx;
    charge_ = charge;
    isGen_ = isGen;
}

kaon::kaon(TVector3 momentum, double eta, double phi, int charge, double mass, bool isGen)
{
    //TVector3 momentum(px,py,pz);
    double energy = sqrt(TMath::Power(mass,2) + TMath::Power(momentum.Mag(),2));
    TLorentzVector PtEtaPhiE(momentum.Perp(),eta,phi,energy);
    momentum_ = momentum;
    PtEtaPhiE_ = PtEtaPhiE;
    charge_ = charge;
    dedx_ = 0;
    isGen_ = isGen;
}

double kaon::deltaR(TLorentzVector otherLV)
{
    return PtEtaPhiE_.DeltaR(otherLV);
}

double kaon::getRapidity()
{
    TLorentzVector LV(momentum_.X(), momentum_.Y(), momentum_.Z(), PtEtaPhiE_.E());

    return LV.Rapidity();
}

bool kaon::matched(kaon genKaon)
{
    double genPt = -999;
    double trkPt = -999;
    //check which kaon is gen level. If neither then throw exception
    if(genKaon.getIsGen())
    {
        genPt = genKaon.getPt();
        trkPt = momentum_.Perp();
    }
    else if(isGen_)
    {
        genPt = momentum_.Perp();
        trkPt = genKaon.getPt();
    }
    else
        throw std::invalid_argument("Neither of these is a gen level kaon");

    double dpt = genPt - trkPt;

    if(kaon::deltaR(genKaon.getLorentzVect()) < 0.1 && fabs(dpt/genPt) < 0.1)
        return true;
}

bool kaon::getIsGen()
{
    return isGen_;
}

double kaon::getP()
{
    return momentum_.Mag();
}

double kaon::getPt()
{
    return momentum_.Perp();
}

double kaon::getPx()
{
    return momentum_.X();
}

double kaon::getPy()
{
    return momentum_.Y();
}

double kaon::getPz()
{
    return momentum_.Z();
}

TVector3 kaon::getMomentumVect()
{
    return momentum_;
}

TLorentzVector kaon::getLorentzVect()
{
    return PtEtaPhiE_;
}

double getDedx()
{
    return dedx_;
}

double getEta()
{
    return PtEtaPhiE_.Eta();
}

double getPhi()
{
    return PtEtaPhiE_.Phi();
}

double getEnergy()
{
    return PtEtaPhiE_.Energy();
}

/*
 * Setters
 */
void setIsGen(bool isGen)
{
    isGen_ = isGen;
}

void setDedx(double dedx)
{
    dedx_ = dedx;
}

void setCharge(int charge)
{
    charge_ = charge;
}

void setMomentumVect(TVector3 momentum)
{
    momentum_ = momentum;
}

void setLorentzVect(TLorentzVector PtEtaPhiE)
{
    PtEtaPhiE_ = PtEtaPhiE;
}
