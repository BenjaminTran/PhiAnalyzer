#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/phi.h"

PhiMeson::PhiMeson()
{
}

PhiMeson::PhiMeson(double mass, double pt, double eta, double phi, double rapidity, double p, bool isMatched)
{
    mass_ = mass;
    pt_ = pt;
    eta_ = eta;
    rapidity_ = rapidity;
    phi_ = phi;
    p_ = p;
    isMatched_ = isMatched;
}

std::vector<PhiMeson> PhiMeson::EventCombinatorialPhi(std::vector<kaon> PKp_, std::vector<kaon> PKm_)
{
    std::vector<PhiMeson> phiCollection;
    for(kaon Pkp : PKp_)
    {
        for(kaon Pkm : PKm_)
        {
            PhiMeson pgf = BuildPhi(Pkp, Pkm);
            PhiMeson pgf(mass,pt,eta,phi,rapidity,p);

            phiCollection.push_back(pgf);
        }
    }

    return phiCollection;
}

PhiMeson PhiMeson::BuildPhi(kaon Pkp, kaon Pkm, bool isMatched)
{
    TVector3 dau1p(Pkp.getPx(), Pkp.getPy(), Pkp.getPz());
    TVector3 dau2p(Pkm.getPx(), Pkm.getPy(), Pkm.getPz());
    TVector3 dauPsum(dau1p + dau2p);
    double mass = sqrt(TMath::Power(Pkp.getEnergy()+Pkm.getEnergy(),2) - dauPsum.Mag2());

    double p = sqrt(dauPsum.Mag2());
    double pt = dauPsum.Perp();

    TLorentzVector phiLV(dauPsum,Pkp.getEnergy() + Pkm.getEnergy());

    double rapidity = phiLV.Rapidity();
    double eta = phiLV.Eta();
    double phi = phiLV.Phi();

    PhiMeson pgf(mass,pt,eta,phi,rapidity,p,isMatched);
    pgf.addKaonDau(Pkp);
    pgf.addKaonDau(Pkm);

    return pgf;
}

double PhiMeson::getMass()
{
    return mass_;
}

double PhiMeson::getEta()
{
    return eta_;
}

double PhiMeson::getRapidity()
{
    return rapidity_;
}

double PhiMeson::getPhi()
{
    return phi_;
}

double PhiMeson::getPt()
{
    return pt_;
}

double PhiMeson::getP()
{
    return p_;
}

bool PhiMeson::getIsGen()
{
    return isMatched_;
}

kaon PhiMeson::getKaonDau(int dauID)
{
    return KaonDau_.at(dauID);
}

void PhiMeson::setMass(double mass)
{
    mass_ = mass;
}

void PhiMeson::setEta(double eta)
{
    eta_ = eta;
}

void PhiMeson::setRapidity(double rapidity)
{
    rapidity_ = rapidity;
}

void PhiMeson::setPhi(double phi)
{
    phi_ = phi;
}

void PhiMeson::setPt(double pt)
{
    pt_ = pt;
}

void PhiMeson::setP(double p)
{
    p_ = p;
}

void PhiMeson::setIsGen(bool isMatched)
{
    isMatched_ = isMatched;
}

void addKaonDau(kaon dau)
{
    KaonDau_.push_back(dau);
}
