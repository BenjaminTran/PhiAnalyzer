// -*- C++ -*-
//
// Package:    PhiSelector/PhiSelector
// Class:      PhiSelector
//
/**\class PhiSelector PhiSelector.cc PhiSelector/PhiSelector/plugins/PhiSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Benjamin Tran
//         Created:  Tue, 30 Jan 2018 19:09:22 GMT
//
//

#include "/afs/cern.ch/user/b/btran/work/CMSSW_8_0_24/src/PhiAnalyzer/PhiAnalyzer/interface/PhiSelector.h"
#include "PhiAnalyzer/PhiAnalyzer/interface/utility.h"

using namespace std;

//
// constructors and destructor
//
PhiSelector::PhiSelector(const edm::ParameterSet& iConfig)
{
    multMin_ = iConfig.getUntrackedParameter<int>("multMin");
    multMax_ = iConfig.getUntrackedParameter<int>("multMax");
    trackPtCut_ = iConfig.getUntrackedParameter<bool>("trackPtCut");
    _trkSrc = consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trkSrc"));
    _vtxSrc = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vtxSrc"));
    _Dedx_Harmonic2 = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));
}


PhiSelector::~PhiSelector()
{
}


//
// member functions
//

void
PhiSelector::DeDxFiller(utility::track_combo track_combo_, edm::Handle<edm::ValueMap<reco::DeDxData> > DeDxTrack, TH2D* dedx_p)
{
        double dedx     = -999.9;
        double momentum = track_combo_.track->p();
        dedx = getDeDx(track_combo_,DeDxTrack);
        if(AcceptTrack(track_combo_,DeDxTrack))
            dedx_p->Fill(momentum,dedx);

}

double
PhiSelector::getDeDx(utility::track_combo track_combo_, edm::Handle<edm::ValueMap<reco::DeDxData> > DeDxTrack)
{
    double dedx_ = -1;
    const edm::ValueMap<reco::DeDxData> dedxTrack = *DeDxTrack.product();
    dedx_ = dedxTrack[track_combo_.track_ref].dEdx();

    return dedx_;
}

void
PhiSelector::FillKaonContainer(utility::track_combo track_combo_, edm::Handle<edm::ValueMap<reco::DeDxData> > DeDxTrack, std::vector<kaon> &pkp, std::vector<kaon> &pkm)
{
    if(!AcceptTrack(track_combo_, DeDxTrack))
        return;

    double energy = sqrt(TMath::Power(kaonMass,2) + TMath::Power(track_combo_.track->p(),2));
    kaon pk(track_combo_.track->p(), getDeDx(track_combo_,DeDxTrack), energy, track_combo_.track->charge());

    //positive kaons
    if(track_combo_.track->charge() == 1)
        pkp.push_back(pk);

    //negative kaons
    if(track_combo_.track->charge() == -1)
        pkm.push_back(pk);
}

void
PhiSelector::CombinatorialMass(std::vector<PhiSelector::kaon> PKp, std::vector<PhiSelector::kaon> PKm, TH1D* h_mass_)
{
    for(PhiSelector::kaon Pkp : PKp)
    {
        for(PhiSelector::kaon Pkm : PKm)
        {
            double mass = sqrt(2*TMath::Power(kaonMass,2) + 2*(Pkp.energy*Pkm.energy - Pkp.p*Pkm.p));
            h_mass_->Fill(mass);
        }
    }
}

bool
PhiSelector::AcceptTrack(utility::track_combo track_combo_, edm::Handle<edm::ValueMap<reco::DeDxData> > DeDxTrack)
{
    double functionValueTop = -999;
    double functionValueBot = -999;
    double momentum = track_combo_.track->p();
    double dedx = getDeDx(track_combo_, DeDxTrack);
    int nhits = track_combo_.track->numberOfValidHits();
    functionValueTop = 0.55*(TMath::Power(1.6/momentum,2) - 2*TMath::Power(0.6/momentum,1)) + 3.3;
    functionValueBot = 0.55*(TMath::Power(1.15/momentum,2) - 2*TMath::Power(0.6/momentum,1)) + 3;
    if(dedx < functionValueTop && dedx > functionValueBot && nhits > 11)
        return true;
    else
        return false;
}

// ------------ method called for each event  ------------
void
PhiSelector::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //Vectors to hold kaons to perform combinatorial mass reconstruction. Following PEN naming scheme
   std::vector<PhiSelector::kaon> PKp_Harm; //Positive charged
   std::vector<PhiSelector::kaon> PKm_Harm; //Negative charged

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

   // Multiplicity selection
   int mult = utility::Multiplicity(tracks,vertices);
   h_mult->Fill(mult);

   for(reco::TrackCollection::const_iterator it = tracks->begin();
           it != tracks->end();
           ++it)
   {
       //Use only tracks good enough to pass utility::Multiplicity
       if(!utility::isTrackGood(it,vertex,trackPtCut_)) continue;

       reco::TrackRef track_ref = reco::TrackRef(tracks,it - tracks->begin());
       utility::track_combo track_bundle(it,track_ref);

       //Fill in the dEdx histograms
       DeDxFiller(track_bundle,DeDx_Harm,h_Dedx_p_Harm);

       // Make the vector of kaons to calculate invariant mass at the end
       FillKaonContainer(track_bundle,DeDx_Harm,PKp_Harm,PKm_Harm);
   }
   CombinatorialMass(PKp_Harm,PKm_Harm,h_mass_Harm);
}


// ------------ method called once each job just before starting event loop  ------------
void
PhiSelector::beginJob()
{
    TH1::SetDefaultSumw2();
    edm::Service<TFileService> fs;

    h_nEvt = fs->make<TH1D>("nEvt","",10,0,10);
    h_mult = fs->make<TH1D>("mult","",400,0,400);
    h_mass_Harm = fs->make<TH1D>("mass_harm",";GeV",800,1.010,1.030);
    h_Dedx_p_Harm = fs->make<TH2D>("Dedx_harm",";p;dE/dx",200,0,5,250,0,15);
}

// ------------ method called once each job just after ending the event loop  ------------
void
PhiSelector::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhiSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.addUntracked<edm::InputTag>("trkSrc",edm::InputTag("generalTracks"));
  desc.addUntracked<edm::InputTag>("vtxSrc",edm::InputTag("offlinePrimaryVertices"));
  desc.addUntracked<int>("multMin",0);
  desc.addUntracked<int>("multMax",999);
  desc.addUntracked<bool>("trackPtCut",true);
  descriptions.add("PhiSelector",desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhiSelector);
