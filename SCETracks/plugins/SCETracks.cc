// -*- C++ -*-
//
// Package:    tracks/SCETracks
// Class:      SCETracks
// 
/**\class SCETracks SCETracks.cc tracks/SCETracks/plugins/SCETracks.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Tue, 23 Feb 2016 15:12:55 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Utilities/interface/InputTag.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class SCETracks : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit SCETracks(const edm::ParameterSet&);
      ~SCETracks();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

  int indexEvent;
  edm::EDGetTokenT<edm::View<reco::Track> > generalTracksToken_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SCETracks::SCETracks(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");

   indexEvent = 0;
   generalTracksToken_ = consumes<edm::View<reco::Track> >(edm::InputTag("generalTracks"));

}


SCETracks::~SCETracks()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SCETracks::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;



   //#ifdef THIS_IS_AN_EVENT_EXAMPLE
   //   Handle<ExampleData> pIn;
   //   iEvent.getByLabel("example",pIn);
   //#endif
   
   //#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   //   ESHandle<SetupData> pSetup;
   //   iSetup.get<SetupRecord>().get(pSetup);
   //#endif

   std::cout<<"Event "<< indexEvent<< std::endl;

   Handle<edm::View<reco::Track> > trackHandle_;
   iEvent.getByToken(generalTracksToken_,trackHandle_);

   for (int j = 0 ; j < (int)trackHandle_->size(); j++){
     const reco::Track& track = trackHandle_->at(j);
     std::cout << "    Track " << j << " " << track.charge()*track.pt() << " " << track.phi()
               << " " << track.eta() << " " << track.dxy() << " " << track.dz() << std::endl;
   }
   indexEvent++;


}


// ------------ method called once each job just before starting event loop  ------------
void 
SCETracks::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SCETracks::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SCETracks::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SCETracks);
