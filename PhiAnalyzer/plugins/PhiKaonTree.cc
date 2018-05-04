#include "PhiAnalyzer/PhiAnalyzer/interface/PhiKaonTree.h"


using namespace std;

PhiKaonTree::PhiKaonTree(const edm::ParameterSet& iConfig)
{
    _trkSrc = consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trkSrc"));
    _vtxSrc = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vtxSrc"));
    _Dedx_Harmonic2 = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));
}

PhiKaonTree::~PhiKaonTree()
{
}


void
PhiKaonTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    h_nEvt->Fill(1);
}

void
PhiKaonTree::beginJob()
{
    TH1::SetDefaultSumw2();
    edm::Service<TFileService> fs;

    h_nEvt = fs->make<TH1D>("nEvt","",10,0,10);
    phiKaonTree = fs->make<TTree>("PhiKaonTree","PhiKaonTree");
}

void
PhiKaonTree::endJob()
{}

void
PhiKaonTree::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    edm::ParameterSetDescription desc;
    desc.addUntracked<edm::InputTag>("trkSrc",edm::InputTag("generalTracks"));
    desc.addUntracked<edm::InputTag>("vtxSrc",edm::InputTag("offlinePrimaryVertices"));
    descriptions.add("PhiKaonTree",desc);
}

DEFINE_FWK_MODULE(PhiKaonTree);
