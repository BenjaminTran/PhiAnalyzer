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
        double ptError_ = -999;
        double dz_ = -999;
        double dzError_ = -999;
        double dxy_ = -999;
        double dxyError_ = -999;
        int nhits_ = -999;
        double chi2_ = -999;
        double chi2norm_ = -999;
        double vx_ = -999;
        double vy_ = -999;
        double vz_ = -999;
        double ndof_ = -999;

    public:
        struct cutVariables{
            double ptError;
            double dz;
            double dzError;
            double dxy;
            double dxyError;
            int nhits;
            double chi2;
            double chi2norm;
            double vx;
            double vy;
            double vz;
            double ndof;
        };
        kaon(TVector3 momentum, double eta, double phi, int charge = -9, double dedx = 0, bool isGen = false, double mass = utility::kaonMass);
        kaon(TVector3 momentum, double eta, double phi, int charge = -9, bool isGen = false, double mass = utility::kaonMass);
        kaon(TVector3 momentum, double eta, double phi, cutVariables cutValues, int charge = -9, double dedx = 0, bool isGen = false, double mass = utility::kaonMass);

        double deltaR(TLorentzVector otherLV);
        bool matched(kaon genKaon);

        double getP();
        double getPt();
        double getPx();
        double getPy();
        double getPz();
        double getPtError();
        double getDz();
        double getDzError();
        double getDxy();
        double getDxyError();
        int getNhits();
        double getChi2();
        double getChi2norm();
        double getVx();
        double getVy();
        double getVz();
        double getNdof();
        TVector3 getMomentumVect();
        TLorentzVector getLorentzVect();
        double getDedx();
        double getRapidity();
        double getEta();
        double getPhi();
        double getEnergy();
        double getCharge();
        bool getIsGen();

        void setMomentumVect(TVector3 momentum);
        void setLorentzVect(TLorentzVector PtEtaPhiE);
        void setDedx(double dedx);
        void setCharge(int charge);
        void setIsGen(bool isGen);
};
#endif
