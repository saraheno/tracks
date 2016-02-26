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
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

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
  //  edm::EDGetTokenT<edm::View<reco::GenParticle> genParticlesToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

  edm::Service<TFileService> fs;
  TH1F *h1;

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
   genParticlesToken_ = consumes<reco::GenParticleCollection>(edm::InputTag("genParticles"));

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

   std::cout<<"Event "<< indexEvent<< std::endl;

   // printing header
   std::cout<<" pT phi eta"<<std::endl;


   // tracks
   Handle<edm::View<reco::Track> > trackHandle_;
   iEvent.getByToken(generalTracksToken_,trackHandle_);


   h1->Fill((int)trackHandle_->size());
   for (int j = 0 ; j < (int)trackHandle_->size(); j++){
     const reco::Track& track = trackHandle_->at(j);
     std::cout << "    Track " << j << " " << track.pt() << " " << track.phi()
               << " " << track.eta() << " " << track.dxy() << " " << track.dz() <<std::endl;
   }


   // genparticles

   
   Handle<reco::GenParticleCollection > GenParticleHandle_;
   iEvent.getByToken(genParticlesToken_,GenParticleHandle_);

   for (int j = 0 ; j < (int)GenParticleHandle_->size(); j++){
          const reco::GenParticle& genparticle = GenParticleHandle_->at(j);
          std::cout << "    GenParticle " << j << " "<<genparticle.pdgId()<<" "<<genparticle.pt()<<" "<<genparticle.phi()<<
	    " "<<genparticle.eta()<<" "<<genparticle.vx()<<" "<<genparticle.vy()<<" "<<
	    sqrt(pow(genparticle.vx(),2)+pow(genparticle.vy(),2))<<" "<<genparticle.vz()
               <<  std::endl;
   }


   // for each generator particle, fine the nearest track
   std::vector<int> pttrk(trackHandle_->size());
   for (int i = 0 ; i < (int)trackHandle_->size(); i++){
     const reco::Track& track = trackHandle_->at(i);
     //     std::cout<<" track "<<i<<std::endl;
     pttrk[i]=0;
     float delR=99999.;
     for (int j = 0 ; j < (int)GenParticleHandle_->size(); j++){
       //       std::cout<<" particle "<<j<<std::endl;
       const reco::GenParticle& genparticle = GenParticleHandle_->at(j);
       if(genparticle.numberOfDaughters()==0) {  // if a final state particle (is this really the right flag?)
	 if(genparticle.charge()!=0) {
           float delphi = (track.phi()-genparticle.phi());
           if(delphi>3.14159) delphi=2.*3.14159-delphi;
           float deleta = (track.eta()-genparticle.eta());
           float DR = sqrt(pow(delphi,2)+pow(deleta,2));
           if(DR<delR) {
  	     delR=DR;
  	     pttrk[i]=j;
           }  // end dR test
	 }  // end zero daughters
       }  //end final state
     } //end loop over gen particles
     std::cout<<" track "<<i<<" matches with genparticle "<<pttrk[i]<<" with delR of "<<delR<<std::endl;
   }
   

   // next event
   indexEvent++;


}


// ------------ method called once each job just before starting event loop  ------------
void 
SCETracks::beginJob()
{
  h1 = fs->make<TH1F>("h1" , "snumber of tracks", 100, 0, 100.);
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
