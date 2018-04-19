#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/phi.h"

PhiMeson::PhiMeson()
{
}

PhiMeson::PhiMeson(double mass, TVector3 momentum, TLorentzVector PtEtaPhiE, bool isGen)
{
    mass_ = mass;
    isGen_ = isGen;

    momentum_ = momentum;
    PtEtaPhiE_ = PtEtaPhiE;
}

std::vector<PhiMeson> PhiMeson::EventCombinatorialPhi(std::vector<kaon> PKp_, std::vector<kaon> PKm_)
{
    std::vector<PhiMeson> phiCollection;
    for(kaon Pkp : PKp_)
    {
        for(kaon Pkm : PKm_)
        {
            PhiMeson pgf = BuildPhi(Pkp, Pkm);

            phiCollection.push_back(pgf);
        }
    }

    return phiCollection;
}

PhiMeson PhiMeson::BuildPhi(kaon Pkp, kaon Pkm, bool isGen)
{
    TVector3 dau1p(Pkp.getPx(), Pkp.getPy(), Pkp.getPz());
    TVector3 dau2p(Pkm.getPx(), Pkm.getPy(), Pkm.getPz());
    TVector3 dauPsum(dau1p + dau2p);
    double mass = sqrt(TMath::Power(Pkp.getEnergy()+Pkm.getEnergy(),2) - dauPsum.Mag2());

    double p = sqrt(dauPsum.Mag2());
    double pt = dauPsum.Perp();

    TLorentzVector phiLV(dauPsum,Pkp.getEnergy() + Pkm.getEnergy());

    TLorentzVector PtEtaPhiE(pt,phiLV.Eta(),phiLV.Phi(),Pkp.getEnergy() + Pkm.getEnergy());

    PhiMeson pgf(mass,dauPsum,PtEtaPhiE,isGen);
    pgf.addKaonDau(Pkp*);
    pgf.addKaonDau(Pkm*);

    return pgf;
}

kaon PhiMeson::getKaonDau(int dauID)
{
    return KaonDau_.at(dauID);
}

void addKaonDau(kaon dau)
{
    KaonDau_.push_back(dau);
}
