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
    return mass;
}

double PhiMeson::getEta()
{
    return eta;
}

double PhiMeson::getRapidity()
{
    return rapidity;
}

double PhiMeson::getPhi()
{
    return phi;
}

double PhiMeson::getPt()
{
    return pt;
}

double PhiMeson::getP()
{
    return p;
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
