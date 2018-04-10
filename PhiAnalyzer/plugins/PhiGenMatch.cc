#include "PhiAnalyzer/PhiAnalyzer/interface/PhiGenMatch.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/Phi.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"

using namespace std;

//Constructors and destructor
PhiGenMatch::PhiGenMatch(const edm::ParameterSet& iConfig)
{
    _genCollection = consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genCollection"));
    _trkSrc = consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trkSrc"));
    _vertexCollName = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertexCollName"));
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
    edm::Handle<reco::TrackCollection> trkSrc;
    edm::Handle<reco::VertexCollection> vertices;
    try{
        utility::GetCollection<reco::GenParticleCollection>(iEvent,_genCollection,gencand);
        utility::GetCollection<reco::TrackCollection>(iEvent,_trkSrc,trkSrc);
        utility::GetCollection<reco::VertexCollection>(iEvent,_vertexCollName,vertices);
    }
    catch(const std::invalid_argument& e){
        std::cerr << e.what();
        return;
    }

    //utility::myVertex vertex = utility::MyVertexBuild(vertices);

    for(reco::TrackCollection::const_iterator trk = trkSrc->begin();
            trk != trkSrc->end();
            ++trk)
    {
        //if(!utility::SelectionCut(trk,vertex,false,1.0,1.0,2.4,0,5))
    }

    for(reco::GenParticleCollection::const_iterator gncand = gencand->begin();
            gncand != gencand->end();
            ++gncand)
    {
        int id = gncand->pdgId();
        double rap = gncand->rapidity();
        double eta = gncand->eta();

        if(fabs(id) == 333 && fabs(rap) < 1.0)
        {
            h_phi_yield_rap_1->Fill(gncand->mass());
            const reco::Candidate *d1 = gncand->daughter(0);
            const reco::Candidate *d2 = gncand->daughter(1);

            if(fabs(d1->pdgId()) != fabs(d2->pdgId())) continue;

            h_phid1_mass->Fill(d1->mass());
            h_phid2_mass->Fill(d2->mass());
        }

        if(fabs(id) == 333 && fabs(eta) < 2.4)
            h_phi_yield_norap->Fill(gncand->mass());
    }
}

void
PhiGenMatch::beginJob()
{
    TH1::SetDefaultSumw2();

    h_nEvt = fs->make<TH1D>("h_nEvt","Events",10,0,10);
    h_phi_yield_rap_1 = fs->make<TH1D>("h_phi_yield_rap_1","Gen Phi Yield rap < 1",50,1,1.05);
    h_phi_yield_norap = fs->make<TH1D>("h_phi_yield_norap","Gen Phi Yield no rap",50,1,1.05);
    h_phid1_mass = fs->make<TH1D>("h_phid1_mass","Dau1 mass",100,0.45,0.55);
    h_phid2_mass = fs->make<TH1D>("h_phid2_mass","Dau2 mass",100,0.45,0.55);

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
    desc.addUntracked<edm::InputTag>("trkSrc",edm::InputTag("generalTracks"));
    desc.addUntracked<edm::InputTag>("vertexCollName",edm::InputTag("offlinePrimaryVertices"));
    descriptions.add("PhiGenMatch",desc);
}

DEFINE_FWK_MODULE(PhiGenMatch);
