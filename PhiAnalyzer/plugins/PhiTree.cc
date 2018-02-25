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
    dedxConstraint_ = iConfig.getUntrackedParameter<string>("dedxConstraint");
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

    mult = utility::Multiplicity(tracks, vertices);
    multRaw = tracks->size();

    h_mult->Fill(mult);
    h_multRaw->Fill(multRaw);

    for(reco::TrackCollection::const_iterator it = tracks->begin();
            it != tracks->end();
            ++it)
    {
        if(!it->quality(reco::TrackBase::highPurity)) continue;
        reco::TrackRef track_ref = reco::TrackRef(tracks,it - tracks->begin());
        utility::track_combo track_bundle(it,track_ref);

        try
        {
            if(!utility::AcceptTrackDeDx(track_bundle,DeDx_Harm,dedxConstraint_)) continue;
        }
        catch(const std::invalid_argument& e)
        {
            std::cerr << e.what();
        }

        double energy = sqrt(TMath::Power(utility::kaonMass,2) + TMath::Power(it->p(),2));

        TLorentzVector tmpLV = new TLorentzVector(it->px(),it->py(),it->pz(),energy);

        double rap = tmpLV.rapidity();

        track_particle_.momentum.push_back(it->p());
        track_particle_.px.push_back(it->px());
        track_particle_.py.push_back(it->py());
        track_particle_.pz.push_back(it->pz());
        track_particle_.pt.push_back(it->pt());
        track_particle_.ptError.push_back(it->ptError());
        track_particle_.energy.push_back(energy);
        track_particle_.dedx.push_back(utility::getDeDx(track_bundle,DeDx_Harm));
        track_particle_.charge.push_back(it->charge());
        track_particle_.dz.push_back(it->dz(vertex.bestvtx));
        track_particle_.dzError.push_back(sqrt(TMath::Power(it->dzError(),2) + TMath::Power(vertex.bestvzError,2)));
        track_particle_.dxy.push_back(it->dxy(vertex.bestvtx));
        track_particle_.dxyError.push_back(sqrt(TMath::Power(it->d0Error(),2) + vertex.bestvxError*vertex.bestvyError));
        track_particle_.eta.push_back(it->eta());
        track_particle_.rapidity.push_back(rap);
        track_particle_.phi.push_back(it->phi());
        track_particle_.ndof.push_back(it->ndof());
        track_particle_.vx.push_back(it->vx());
        track_particle_.vy.push_back(it->vy());
        track_particle_.vz.push_back(it->vz());
        track_particle_.vzFlip.push_back(-(it->vz()));
        track_particle_.chi2.push_back(it->chi2());
        track_particle_.chi2norm.push_back(it->normalizedChi2());
        track_particle_.nhits.push_back(it->numberOfValidHits());
    }

    trackTree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void
PhiTree::beginJob()
{
    TH1::SetDefaultSumw2();
    edm::Service<TFileService> fs;

    trackTree = fs->make<TTree>("TrackTree","TrackTTree");

    h_nEvt = fs->make<TH1D>("Evt","",10,0,10);
    h_mult = fs->make<TH1D>("mult","",900,0,900);
    h_multRaw = fs->make<TH1D>("multRaw","",5000,0,5000);

    trackTree->Branch("momentum" , &track_particle_.momentum , "momentum/F");
    trackTree->Branch("px"       , &track_particle_.px       , "px/F");
    trackTree->Branch("py"       , &track_particle_.py       , "py/F");
    trackTree->Branch("pz"       , &track_particle_.pz       , "pz/F");
    trackTree->Branch("pt"       , &track_particle_.pt       , "pt/F");
    trackTree->Branch("ptError"  , &track_particle_.ptError  , "ptError/F");
    trackTree->Branch("energy"   , &track_particle_.energy   , "energy/F");
    trackTree->Branch("dedx"     , &track_particle_.dedx     , "dedx/F");
    trackTree->Branch("dz"       , &track_particle_.dz       , "dz/F");
    trackTree->Branch("dzError"  , &track_particle_.dzError  , "dzError/F");
    trackTree->Branch("dxy"      , &track_particle_.dxy      , "dxy/F");
    trackTree->Branch("dxyError" , &track_particle_.dxyError , "dxyError/F");
    trackTree->Branch("eta"      , &track_particle_.eta      , "eta/F");
    trackTree->Branch("rapidity" , &track_particle_.rapidity , "rapidity/F");
    trackTree->Branch("phi"      , &track_particle_.phi      , "phi/F");
    trackTree->Branch("vx"       , &track_particle_.vx       , "vx/F");
    trackTree->Branch("vy"       , &track_particle_.vy       , "vy/F");
    trackTree->Branch("vz"       , &track_particle_.vz       , "vz/F");
    trackTree->Branch("vzFlip"   , &track_particle_.vzFlip   , "vzFlip/F");
    trackTree->Branch("chi2"     , &track_particle_.chi2     , "chi2/F");
    trackTree->Branch("chi2norm" , &track_particle_.chi2norm , "chi2norm/F");
    trackTree->Branch("ndof"     , &track_particle_.ndof     , "ndof/F");
    trackTree->Branch("nhits"    , &track_particle_.nhits    , "nhits/F");
    trackTree->Branch("charge"   , &track_particle_.charge   , "charge/F");
    trackTree->Branch("mult"     , &mult                     , "mult/F");
    trackTree->Branch("multRaw"  , &multRaw                  , "multRaw/F");

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
  desc.addUntracked<string>("dedxConstraint","loose");
  descriptions.add("PhiTree",desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhiTree);
