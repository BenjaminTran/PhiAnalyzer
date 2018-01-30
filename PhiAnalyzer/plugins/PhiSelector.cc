// -*- C++ -*-
//
// Package:    PhiAnalyzer/PhiAnalyzer
// Class:      PhiAnalyzer
// 
/**\class PhiAnalyzer PhiAnalyzer.cc PhiAnalyzer/PhiAnalyzer/plugins/PhiAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Benjamin Tran
//         Created:  Tue, 30 Jan 2018 19:09:22 GMT
//
//

#include $CMSSW_BASE/src/PhiAnalyzer/PhiAnalyzer/interface/PhiSelector.h


//
// constructors and destructor
//
PhiAnalyzer::PhiAnalyzer(const edm::ParameterSet& iConfig)
 :
  trackTags_(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))

{
   //now do what ever initialization is needed
   usesResource("TFileService");
@example_histo   edm::Service<TFileService> fs;
@example_histo   histo = fs->make<TH1D>("charge" , "Charges" , 200 , -2 , 2 );

}


PhiAnalyzer::~PhiAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PhiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
@example_histo       histo->Fill( charge );
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
PhiAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PhiAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PhiAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  //edm::ParameterSetDescription desc;
  //desc.setUnknown();
  //descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  ParameterSetDescription desc;
  desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  desc.addUntracked<bool>("test","False");
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PhiAnalyzer);
