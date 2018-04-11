#ifndef PHIANALYZER__KAON_H
#define PHIANALYZER__KAON_H

#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"

class kaon{
    private:
        //double p;
        //double pt;
        //double px;
        //double py;
        //double pz;
        TVector3 momentum_;
        TLorentzVector PtEtaPhiE_;
        double dedx_;
        int charge_;
        bool isGen_;

    public:
        kaon(TVector3 momentum, double eta, double phi, int charge = -9, double dedx = 0, bool isGen = false, double mass = utility::kaonMass);
        kaon(TVector3 momentum, double eta, double phi, int charge = -9, bool isGen = false, double mass = utility::kaonMass);

        double deltaR(TLorentzVector otherLV);
        bool matched(kaon genKaon);

        double getP();
        double getPt();
        double getPx();
        double getPy();
        double getPz();
        TVector3 getMomentumVect();
        TLorentzVector getLorentzVect();
        double getDedx();
        double getRapidity();
        double getEta();
        double getPhi();
        double getEnergy();
        bool getIsGen();

        void setMomentumVect(TVector3 momentum);
        void setLorentzVect(TLorentzVector PtEtaPhiE);
        void setDedx(double dedx);
        void setCharge(int charge);
        void setIsGen(bool isGen);
};
#endif
