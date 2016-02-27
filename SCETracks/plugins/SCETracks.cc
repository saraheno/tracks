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

// Root objects                                                                                                                              
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TFile.h"
#include "TProfile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TMath.h"



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
  TH1F *hntrk;
  TH1F *hnchst;
  TH2D *hgenrecopt;

  const int idbg = 1;

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

   if(idbg>0) std::cout<<"Event "<< indexEvent<< std::endl;

   // printing header
   if(idbg>0) std::cout<<" pT phi eta"<<std::endl;


   // tracks
   Handle<edm::View<reco::Track> > trackHandle_;
   iEvent.getByToken(generalTracksToken_,trackHandle_);


   hntrk->Fill((int)trackHandle_->size());
   for (int j = 0 ; j < (int)trackHandle_->size(); j++){
     const reco::Track& track = trackHandle_->at(j);
     if(idbg>0) std::cout << "    Track " << j << " " << track.pt() << " " << track.phi()
               << " " << track.eta() << " " << track.dxy() << " " << track.dz() <<std::endl;
   }


   // genparticles

   
   Handle<reco::GenParticleCollection > GenParticleHandle_;
   iEvent.getByToken(genParticlesToken_,GenParticleHandle_);

   int nchst=0; // count numbe of charged stable particles
   for (int j = 0 ; j < (int)GenParticleHandle_->size(); j++){
     const reco::GenParticle& genparticle = GenParticleHandle_->at(j);
     if(idbg>0) std::cout << "    GenParticle " << j << " "<<genparticle.pdgId()<<" "<<genparticle.pt()<<" "<<
            genparticle.phi()<<
	    " "<<genparticle.eta()<<" "<<genparticle.vx()<<" "<<genparticle.vy()<<" "<<
	    sqrt(pow(genparticle.vx(),2)+pow(genparticle.vy(),2))<<" "<<genparticle.vz()
            <<  std::endl;
     if(genparticle.numberOfDaughters()==0) {  // if a final state particle (is this really the right flag?)
       if(genparticle.charge()!=0) {
	 nchst+=1;
       }
     }
   }
   hnchst->Fill(nchst);
   if(idbg>0) std::cout<<" nmber of stable charged is "<<nchst<<std::endl;


   // for each generator particle, fine the nearest track
   std::vector<int> pttrk(GenParticleHandle_->size());
   for (int j = 0 ; j < (int)GenParticleHandle_->size(); j++){
     const reco::GenParticle& genparticle = GenParticleHandle_->at(j);
     pttrk[j]=-1;
     float delR=99999.;
     for (int i = 0 ; i < (int)trackHandle_->size(); i++){
     const reco::Track& track = trackHandle_->at(i);
       if(genparticle.numberOfDaughters()==0) {  // if a final state particle (is this really the right flag?)
	 if(genparticle.charge()!=0) {
           float delphi = (track.phi()-genparticle.phi());
           if(delphi>3.14159) delphi=2.*3.14159-delphi;
           float deleta = (track.eta()-genparticle.eta());
           float DR = sqrt(pow(delphi,2)+pow(deleta,2));
           if(DR<delR) {
  	     delR=DR;
  	     pttrk[j]=i;
           }  // end dR test
	 }  // end zero daughters
       }  //end final state
     } //end loop over gen particles
     if(idbg>0) std::cout<<" genparticle  "<<j<<" matches with track "<<pttrk[j]<<" with delR of "<<delR<<std::endl;
   }
   


   // Scatter plot some gen-sim quantities
   for (int j = 0 ; j < (int)GenParticleHandle_->size(); j++){
     if(idbg>0) std::cout<<" for gen particle "<<j<<" matching to track "<<pttrk[j]<<std::endl;
     const reco::GenParticle& genparticle = GenParticleHandle_->at(j);
     if(pttrk[j]>=0) {
     const reco::Track& track = trackHandle_->at(pttrk[j]);
       if(genparticle.numberOfDaughters()==0) {  // if a final state particle (is this really the right flag?)
	 if(genparticle.charge()!=0) {
	   hgenrecopt->Fill(genparticle.pt(),track.pt());
	 }
       }
     }
   }


   // next event
   indexEvent++;


}


// ------------ method called once each job just before starting event loop  ------------
void 
SCETracks::beginJob()
{
  hntrk = fs->make<TH1F>("hntrk" , "snumber of tracks", 100, 0, 100.);
  hnchst = fs->make<TH1F>("hnchst" , "snumber of gen stable charged tracks", 20, 0, 20.);
  hgenrecopt = fs->make<TH2D>("hgenrecopt" , "gen vs reco pt", 100, 0, 100.,100,0.,100.);
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
