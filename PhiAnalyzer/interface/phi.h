#ifndef PHIANALYZER__PHIMESON_H
#define PHIANALYZER__PHIMESON_H
#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"

class PhiMeson{
    private:
        double mass_;
        double eta_;
        double rapidity_;
        double phi_;
        double pt_;
        double p_;

    public:
        PhiMeson();
        PhiMeson(double mass, double pt, double eta, double phi, double rapidity, double p);
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
