#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/Phi.h"

PhiMeson::PhiMeson()
{
}

PhiMeson::PhiMeson(double mass_, double pt_, double eta_, double phi_, double rapidity_, double p_ = 0)
{
    mass = mass_;
    pt = pt_;
    eta = eta_;
    rapidity = rapidity_;
    phi = phi_;
    eta = p_;
}

PhiMeson::getMass()
{
    return mass;
}

PhiMeson::getEta()
{
    return eta;
}

PhiMeson::getRapidity()
{
    return rapidity;
}

PhiMeson::getPhi()
{
    return phi;
}

PhiMeson::getPt()
{
    return pt;
}

PhiMeson::getP()
{
    return p;
}

PhiMeson::setMass(double mass_)
{
    mass = mass_;
}

PhiMeson::setEta(double eta_)
{
    eta = eta_;
}

PhiMeson::setRapidity(double rapidity_)
{
    rapidity = rapidity_;
}

PhiMeson::setPhi(double phi_)
{
    phi = phi_;
}

PhiMeson::setPt(double pt_)
{
    pt = pt_;
}

PhiMeson::setP(double p_)
{
    p = p_;
}
