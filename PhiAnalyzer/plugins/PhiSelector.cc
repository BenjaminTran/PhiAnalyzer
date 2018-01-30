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


//
// constructors and destructor
//
PhiSelector::PhiSelector(const edm::ParameterSet& iConfig)
 :
  _trkSrc = consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trkSrc"));
  _vtxSrc = consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vtxSrc"));
  _Dedx_Harmonic2 = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("Dedx_Harmonic2"));
  _Dedx_Trunc40 = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("Dedx_Trunc40"));

{
   //now do what ever initialization is needed
   usesResource("TFileService");

}


PhiSelector::~PhiSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PhiSelector::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   using reco::TrackCollection;

    Handle<TrackCollection> tracks;
    iEvent.getByLabel(trackTags_,tracks);
    for(TrackCollection::const_iterator itTrack = tracks->begin();
        itTrack != tracks->end();                      
        ++itTrack) {
       int charge = 0;
       charge = itTrack->charge();  
       histo->Fill(charge);
    }

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
PhiSelector::beginJob()
{
    TH1::SetDefaultSumw2();
    edm::Service<TFileService> fs;

    fs->make<TH1D>("Charges","",10,-5,5);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhiSelector::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhiSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  //edm::ParameterSetDescription desc;
  //desc.setUnknown();
  //descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  edm::ParameterSetDescription desc;
  desc.addUntracked<edm::InputTag>("trkSrc",edm::InputTag("generalTracks"));
  desc.addUntracked<edm::InputTag>("vtxSrc",edm::InputTag("offlinePrimaryVertices"));
  desc.addUntracked<edm::InputTag>("Dedx_Harmonic2",edm::InputTag("dedxHarmonic2"));
  desc.addUntracked<edm::InputTag>("Dedx_Trunc40",edm::InputTag("dedxTruncated40"));
  descriptions.add("PhiSelector",desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhiSelector);
