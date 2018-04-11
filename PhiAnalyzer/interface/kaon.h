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
        kaon(TVector3 momentum, double eta, double phi, double dedx, int charge, double mass, bool isGen);
        kaon(TVector3 momentum, double eta, double phi, int charge, double mass, bool isGen);

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
        double getEta();
        double getPhi();
        double getEnergy();
        bool getIsGen();

        void setMomentumVect();
        void setLorentzVect();
        void setDedx();
        void setCharge();
        void setIsGen();
};
#endif
