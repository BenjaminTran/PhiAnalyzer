#ifndef PHIANALYZER__PHIMESON_H
#define PHIANALYZER__PHIMESON_H
#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"

class PhiMeson{
    private:
        double mass;
        double eta;
        double rapidity;
        double phi;
        double pt;
        double p;

    public:
        PhiMeson();
        PhiMeson(double mass_, double pt_, double eta_, double phi_, double rapidity_, double p_ = 0);
        double getMass();
        double getEta();
        double getRapidity();
        double getPhi();
        double getPt();
        double getP();

        void setMass(double mass_);
        void setEta(double eta_);
        void setRapidity(double rapidity_);
        void setPhi(double phi_);
        void setPt(double pt_);
        void setP(double p_);

};

#endif
