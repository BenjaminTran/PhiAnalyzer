#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/Phi.h"

PhiMeson::PhiMeson()
{
}

PhiMeson::PhiMeson(double mass, double pt, double eta, double phi, double rapidity, double p = 0)
{
    mass_ = mass;
    pt_ = pt;
    eta_ = eta;
    rapidity_ = rapidity;
    phi_ = phi;
    p_ = p;
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
