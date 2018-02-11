// -*- C++ -*-
//
// Package:    PhiTree/PhiTree
// Class:      PhiTree
//
/**\class PhiTree PhiTree.cc PhiTree/PhiTree/plugins/PhiTree.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Benjamin Tran
//         Created:  Mon, 06 Feb 2018 19:09:22 GMT
//
//

#include "/afs/cern.ch/user/b/btran/work/CMSSW_8_0_24/src/PhiAnalyzer/PhiAnalyzer/interface/PhiTree.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/PhiSelector.h"

using namespace std;

//
// constructors and destructor
//
PhiTree::PhiTree(const edm::ParameterSet& iConfig)
{
    multMin_ = iConfig.getUntrackedParameter<int>("multMin");
    multMax_ = iConfig.getUntrackedParameter<int>("multMax");
    _trkSrc = consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trkSrc"));
    _vtxSrc = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vtxSrc"));
    _Dedx_Harmonic2 = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));
}


PhiTree::~PhiTree()
{
}


//
// member functions
//

// ------------ method called for each event  ------------
void
PhiTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    h_nEvt->Fill(1);

    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(_trkSrc,tracks);

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(_vtxSrc,vertices);

    edm::Handle<edm::ValueMap<reco::DeDxData> > DeDx_Harm;
    iEvent.getByToken(_Dedx_Harmonic2,DeDx_Harm);
    if(!DeDx_Harm.isValid())
    {
        cout << "Bad DeDx collection" << endl;
        return;
    }

    utility::myVertex vertex = utility::MyVertexBuild(vertices);

    for(reco::TrackCollection::const_iterator it = tracks->begin();
            it != tracks->end();
            ++it)
    {
        if(!it->quality(reco::TrackBase::highPurity)) continue;
        reco::TrackRef track_ref = reco::TrackRef(tracks,it - tracks->begin());
        utility::track_combo track_bundle(it,track_ref);

        if(!utility::AcceptTrackDeDx(track_bundle,DeDx_Harm,false)) continue;

        track_particle_.momentum      = it->p();
        track_particle_.px            = it->px();
        track_particle_.py            = it->py();
        track_particle_.pz            = it->pz();
        track_particle_.pt            = it->pt();
        track_particle_.ptError       = it->ptError();
        track_particle_.energy        = sqrt(TMath::Power(utility::kaonMass,2) + TMath::Power(it->p(),2));
        track_particle_.dedx          = utility::getDeDx(track_bundle,DeDx_Harm);
        track_particle_.charge        = it->charge();
        track_particle_.dz            = it->dz(vertex.bestvtx);
        track_particle_.dzError       = sqrt(TMath::Power(it->dzError(),2) + TMath::Power(vertex.bestvzError,2));
        track_particle_.dxy           = it->dxy(vertex.bestvtx);
        track_particle_.dxyError      = sqrt(TMath::Power(it->d0Error(),2) + vertex.bestvxError*vertex.bestvyError);
        track_particle_.eta           = it->eta();
        track_particle_.phi           = it->phi();
        track_particle_.ndof          = it->ndof();
        track_particle_.vx            = it->vx();
        track_particle_.vy            = it->vy();
        track_particle_.vz            = it->vz();
        track_particle_.vzFlip        = -track_particle_.vz;
        track_particle_.chi2          = it->chi2();
        track_particle_.chi2norm      = it->normalizedChi2();
        track_particle_.nhits         = it->numberOfValidHits();

        trackTree->Fill();
    }

}


// ------------ method called once each job just before starting event loop  ------------
void
PhiTree::beginJob()
{
    TH1::SetDefaultSumw2();
    edm::Service<TFileService> fs;

    trackTree = fs->make<TTree>("TrackTree","TrackTTree");

    h_nEvt = fs->make<TH1D>("Evt","",10,0,10);

    trackTree->Branch("momentum"      , &track_particle_.momentum      , "momentum/F");
    trackTree->Branch("px"            , &track_particle_.px            , "px/F");
    trackTree->Branch("py"            , &track_particle_.py            , "py/F");
    trackTree->Branch("pz"            , &track_particle_.pz            , "pz/F");
    trackTree->Branch("pt"            , &track_particle_.pt            , "pt/F");
    trackTree->Branch("ptError"       , &track_particle_.ptError       , "ptError/F");
    trackTree->Branch("energy"        , &track_particle_.energy        , "energy/F");
    trackTree->Branch("dedx"          , &track_particle_.dedx          , "dedx/F");
    trackTree->Branch("dz"            , &track_particle_.dz            , "dz/F");
    trackTree->Branch("dzError"       , &track_particle_.dzError       , "dzError/F");
    trackTree->Branch("dxy"           , &track_particle_.dxy           , "dxy/F");
    trackTree->Branch("dxyError"      , &track_particle_.dxyError      , "dxyError/F");
    trackTree->Branch("eta"           , &track_particle_.eta           , "eta/F");
    trackTree->Branch("phi"           , &track_particle_.phi           , "phi/F");
    trackTree->Branch("vx"            , &track_particle_.vx            , "vx/F");
    trackTree->Branch("vy"            , &track_particle_.vy            , "vy/F");
    trackTree->Branch("vz"            , &track_particle_.vz            , "vz/F");
    trackTree->Branch("vzFlip"        , &track_particle_.vzFlip        , "vzFlip/F");
    trackTree->Branch("chi2"          , &track_particle_.chi2          , "chi2/F");
    trackTree->Branch("chi2norm"      , &track_particle_.chi2norm      , "chi2norm/F");
    trackTree->Branch("ndof"          , &track_particle_.ndof          , "ndof/F");
    trackTree->Branch("nhits"         , &track_particle_.nhits         , "nhits/I");
    trackTree->Branch("charge"        , &track_particle_.charge        , "charge/I");
}

// ------------ method called once each job just after ending the event loop  ------------
void
PhiTree::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhiTree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.addUntracked<edm::InputTag>("trkSrc",edm::InputTag("generalTracks"));
  desc.addUntracked<edm::InputTag>("vtxSrc",edm::InputTag("offlinePrimaryVertices"));
  desc.addUntracked<int>("multMin",0);
  desc.addUntracked<int>("multMax",999);
  descriptions.add("PhiTree",desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhiTree);
