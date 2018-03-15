#include "PhiAnalyzer/PhiAnalyzer/interface/PhiGenMatch.h"

using namespace std;

//Constructors and destructor
PhiGenMatch::PhiGenMatch(const edm::ParameterSet& iConfig)
{

}

PhiGenMatch::~PhiGenMatch()
{}

// Member functions

// ------------ method called for each event  ------------
void
PhiGenMatch::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    h_nEvt->Fill(1);

    edm::Handle<reco::GenParticleCollection> gencand;
    try{
        utility::GetCollection(iEvent,_genCollection,gencand);
    }
    catch(const std::invalid_argument& e){
        std::cerr << e.what();
        return;
    }

    for(reco::GenParticleCollection::const_iterator gncand = gencand->begin();
            gncand != gencand->end();
            ++gncand)
    {
        int id = gncand->pdgId();
        double rap = gncand->rapidity();
        double eta = gncand->eta();

        if(fabs(id) == 333 && fabs(rap) < 1.0)
            h_phi_yield_rap_1->Fill(gncand->mass());

        if(fabs(id) == 333 && fabs(eta) < 2.4)
            h_phi_yield_norap->Fill(gncand->mass());
    }

    //iEvent.getByToken(_genCollection,gencand);
    //if(!gencand.isValid())
    //{
        //cout << "GenCollection is invalid" << endl;
        //return;
    //}
}

void
PhiGenMatch::beginJob()
{
    TH1::SetDefaultSumw2();

    h_nEvt = fs->make<TH1D>("h_nEvt","Events",10,0,10);
    h_phi_yield_rap_1 = fs->make<TH1D>("h_phi_yield_rap_1","Gen Phi Yield rap < 1",50,1,1.05);
    h_phi_yield_norap = fs->make<TH1D>("h_phi_yield_norap","Gen Phi Yield no rap",50,1,1.05);

}

void
PhiGenMatch::endJob()
{
}

void
PhiGenMatch::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.addUntracked<edm::InputTag>("genCollection",edm::InputTag("genParticles"));
    descriptions.add("PhiGenMatch",desc);
}

DEFINE_FWK_MODULE(PhiGenMatch);
