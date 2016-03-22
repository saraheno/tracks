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
  TH1F *hntrk,*hnchst,*hnchst1,*hrgen,*hrgen2,*hrreco,*hdR,*hgenreco,*hptgen,*hptreco,*hnchgen,*hnchreco,*hegen,*hereco,*hlgen,*himp;
  TH2D *hgenrecopt,*hratvr,*hphiphi,*hphiphip,*hphiphim,*hphiphi2,*hphid1,*hphid2,*hii,*hii1,*hii2;


  const int idbg = 20;
  const float ptcut = 0.8;
  const float radcut = 40.;
  const float Bmag=3.8;
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
   if(idbg>0) std::cout<<" pT phi eta dxy  dz"<<std::endl;


   // tracks
   Handle<edm::View<reco::Track> > trackHandle_;
   iEvent.getByToken(generalTracksToken_,trackHandle_);


   hntrk->Fill((int)trackHandle_->size());
   for (int j = 0 ; j < (int)trackHandle_->size(); j++){
     const reco::Track& track = trackHandle_->at(j);
     if(idbg>0) std::cout << "    Track " << j << " " << track.pt() << " " << track.phi()
               << " " << track.eta() << " " << track.dxy() << " " << track.dz() <<std::endl;
     himp->Fill(track.dxy());
   }


   // genparticles

   
   Handle<reco::GenParticleCollection > GenParticleHandle_;
   iEvent.getByToken(genParticlesToken_,GenParticleHandle_);

   // energy of original dark pion
   const reco::GenParticle& DarkPion = GenParticleHandle_->at(0);
   float DarkPionE = DarkPion.energy();


   int nchst=0; // count numbe of charged stable particles
   int nchst1=0;  // count ones with pt> 1 GeV
   for (int j = 0 ; j < (int)GenParticleHandle_->size(); j++){
     const reco::GenParticle& genparticle = GenParticleHandle_->at(j);
     float rgen =sqrt(pow(genparticle.vx(),2)+pow(genparticle.vy(),2)); 
     float lgen = sqrt(pow(genparticle.vx(),2)+pow(genparticle.vy(),2)+pow(genparticle.vz(),2)) *genparticle.mass() / genparticle.energy(); 
         if(idbg>0) std::cout << "    GenParticle " << j << " "<<genparticle.pdgId()<<" "<<genparticle.pt()<<" "<<
            genparticle.phi()<<
	    " "<<genparticle.eta()<<" "<<genparticle.vx()<<" "<<genparticle.vy()<<" "<<
		      rgen<<" "<<genparticle.vz()<<" "<<lgen
            <<  std::endl;
     if(genparticle.numberOfDaughters()==0) {  // if a final state particle (is this really the right flag?)
       if(genparticle.charge()!=0) {
	 nchst+=1;
	 if(genparticle.pt()>ptcut) {
	   nchst1+=1;
	 }
       }
     }
   }
   hnchst->Fill(nchst);
   hnchst1->Fill(nchst1);
   if(idbg>0) std::cout<<" nmber of stable charged is "<<nchst<<std::endl;
   if(idbg>0) std::cout<<" nmber of stable charged pt>ptcut is is "<<nchst1<<std::endl;

   const reco::GenParticle& genparticle = GenParticleHandle_->at(1);
   float lgen = sqrt(pow(genparticle.vx()-DarkPion.vx(),2)+pow(genparticle.vy()-DarkPion.vy(),2)+pow(genparticle.vz()-DarkPion.vz(),2)) *DarkPion.mass() / DarkPion.energy(); 
   hlgen->Fill(lgen);

   for (int j = 0 ; j < (int)GenParticleHandle_->size(); j++){
     const reco::GenParticle& genparticle = GenParticleHandle_->at(j);
     float rgen =sqrt(pow(genparticle.vx(),2)+pow(genparticle.vy(),2)); 

     if(genparticle.numberOfDaughters()==0) {  // if a final state particle (is this really the right flag?)
       if(genparticle.charge()!=0) {
	 if(genparticle.pt()>ptcut) {
         hrgen->Fill(rgen);
         hrgen2->Fill(rgen);
	 if(rgen<radcut) {
  	   hptgen->Fill(genparticle.pt());
	   hnchgen->Fill(nchst);
	   hegen->Fill(DarkPionE);
	 }
	 }
       }
     }
   }

   std::cout<<"test"<<std::endl;


   // for each generator particle, fine the nearest track
   std::vector<int> pttrk(GenParticleHandle_->size());
   std::vector<float> dRgt(trackHandle_->size());
   std::vector<int> ptgen(trackHandle_->size());
   for (int i = 0 ; i < (int)trackHandle_->size(); i++){ ptgen[i]=-1;}

   for (int j = 0 ; j < (int)GenParticleHandle_->size(); j++){
     const reco::GenParticle& genparticle = GenParticleHandle_->at(j);
     pttrk[j]=-1;
     float delR=99999.;
     if(genparticle.numberOfDaughters()==0) {  // if a final state particle (is this really the right flag?)
       if(genparticle.charge()!=0) {
         if(genparticle.pt()>ptcut) {
           for (int i = 0 ; i < (int)trackHandle_->size(); i++){
             const reco::Track& track = trackHandle_->at(i);
	     if(genparticle.charge()==track.charge() ) {
             float delphi = (track.phi()-genparticle.phi());
             if(delphi>3.14159) delphi=2.*3.14159-delphi;
             float deleta = (track.eta()-genparticle.eta());
             float DR = sqrt(pow(delphi,2)+pow(deleta,2));
       // try pT matching
       DR = fabs(genparticle.pt()-track.pt())/genparticle.pt();
       // try 1/pT matching
       //       DR = fabs((1/genparticle.pt())-(1/track.pt()))/(1/genparticle.pt());
             if(idbg>9) std::cout<<" j "<<j<<" i "<< i<<" DR "<<DR<<std::endl;
             if(DR<delR) {
	       if(ptgen[i]<0) {
                 delR=DR;
  	         pttrk[j]=i;
               } else {
	         if(DR<dRgt[i]) {
	           pttrk[ptgen[i]]=-1;
	           delR=DR;
	           pttrk[j]=i;
	         }
               } // end test on pointer to gen particle
             }  // end dR test
	     }  // end of charge test
           }  //end loop over tracks
           if(idbg>0) std::cout<<" genparticle  "<<j<<" matches with track "<<pttrk[j]<<" with delR of "<<delR<<std::endl;
         }  // tests on gen part
       }// tests of gen part
     } // tests on the gen part

     if(pttrk[j]>=0) {
       ptgen[pttrk[j]]=j;
       dRgt[pttrk[j]]=delR;
       if(delR>0.6) {
         ptgen[pttrk[j]]=-1;
         pttrk[j]=-1;
       }
     }
       hdR->Fill(delR);
   } //end loop over gen particles
      


   // Scatter plot some gen-sim quantities
   for (int j = 0 ; j < (int)GenParticleHandle_->size(); j++){
     const reco::GenParticle& genparticle = GenParticleHandle_->at(j);
     float rgen =sqrt(pow(genparticle.vx(),2)+pow(genparticle.vy(),2)); 
     if(genparticle.numberOfDaughters()==0) {  // if a final state particle (is this really the right flag?)
       if(genparticle.charge()!=0) {
	 if(genparticle.pt()>ptcut) {
           if(idbg>0) std::cout<<" for gen particle "<<j<<" matching to track "<<pttrk[j]<<std::endl;
           if(pttrk[j]>=0) {
             const reco::Track& track = trackHandle_->at(pttrk[j]);
             hgenrecopt->Fill(genparticle.pt(),track.pt());
             hgenreco->Fill(track.pt()/genparticle.pt());
             hratvr->Fill(rgen,track.pt()/genparticle.pt());


             if(fabs(1-(track.pt()/genparticle.pt()))>0.2) {
	        std::cout<<"danger danger will robinson bad match between gen particle "<<j<<" and track "<<pttrk[j]<<std::endl;
             } else {
               hrreco->Fill(rgen);
	       // compare impact parameter from generator and reconstructed vertex and phi
	       if(idbg>8) std::cout<<" gen charge phi vx vy are "<<genparticle.charge()<<" "<<genparticle.phi()<<" "<<genparticle.vx()<<" "<<genparticle.vy()<<std::endl;
	       if(idbg>8) std::cout<<" reco charge phi vx vy are "<<track.charge()<<" "<<track.phi()<<" "<<track.vx()<<" "<<track.vy()<<std::endl;
	       // calculate radius of curvature in cm
	       float Rg=100*genparticle.pt()/0.3/Bmag;
	       float Rr=100*track.pt()/0.3/Bmag;
	       if(idbg>8) std::cout<<" Rg and Rr are "<<Rg<<" "<<Rr<<std::endl;
	       // calculate x and y at center of circle 

	       float xcenterg = -(genparticle.vx()+genparticle.charge()*Rg*sin(genparticle.phi()));
	       float ycenterg = -(genparticle.vy()-genparticle.charge()*Rg*cos(genparticle.phi()));

	       float xcenterr = -(track.vx()+track.charge()*Rr*sin(track.phi()));
	       float ycenterr = -(track.vy()-track.charge()*Rr*cos(track.phi()));

	       if(idbg>8) std::cout<<" center gen is "<<xcenterg<<" "<<ycenterg<<" rad gen is "<<rgen<<std::endl;
	       if(idbg>8) std::cout<<" center reco is "<<xcenterr<<" "<<ycenterr<<std::endl;

	       // calculate phi at distance of closest approach
	       float phig=0.;
	       if(genparticle.charge()<0) {
		 phig = atan2(xcenterg,-ycenterg);
	       } else {
		 phig = atan2(-xcenterg,ycenterg);
	       }
	       if(idbg>8) std::cout<<" phi at creation, phi at dca, phi reco are "<<genparticle.phi()<<" "<<phig<<" "<<track.phi()<<std::endl;
	       hphiphi->Fill(phig,track.phi());
	       hphid1->Fill(rgen,phig-track.phi());
	       if((idbg>8)&&fabs(phig-track.phi())>2) {
		 std::cout<<" danger danger will robinson"<<std::endl;
	       }
	       hphiphi2->Fill(genparticle.phi(),track.phi());
	       hphid2->Fill(rgen,genparticle.phi()-track.phi());
	       if(genparticle.charge()>0) hphiphip->Fill(phig,track.phi());
	       if(genparticle.charge()<0) hphiphim->Fill(phig,track.phi());

	       // calculate impact parameter
	       float dxyg=sqrt(pow(xcenterg,2)+pow(ycenterg,2))-Rg;
	       if(genparticle.charge()>0) dxyg=-1.*dxyg;
	       hii->Fill(dxyg,track.dxy());
	       hii1->Fill(rgen,track.dxy());
	       hii2->Fill(rgen,dxyg);
	       if(idbg>8) std::cout<<" rgen dxyg tdxy are "<<rgen<<" "<<dxyg<<" "<<track.dxy()<<std::endl;
	       if(rgen>20) std::cout<<"danger will"<<std::endl;

               if(rgen<radcut) {
                 hptreco->Fill(genparticle.pt());
    	         hnchreco->Fill(nchst);
	         hereco->Fill(DarkPionE);
               }  // within r cut cm


             } // momentum measurement easonable	 
           }
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
  hntrk = fs->make<TH1F>("hntrk" , "snumber of tracks", 20, 0, 20.);
  hnchst = fs->make<TH1F>("hnchst" , "snumber of gen stable charged tracks", 20, 0, 20.);
  hnchst1 = fs->make<TH1F>("hnchst1" , "snumber of gen stable charged tracks pt>ptcut", 20, 0, 20.);
  hrgen = fs->make<TH1F>("hrgen" , "radius for track creatin", 100, 0, 500.);
  hlgen = fs->make<TH1F>("hlgen" , "reduced length for track creatin", 100, 0, 50.);
  hrgen2 = fs->make<TH1F>("hrgen2" , "radius for track creatin", 100, 0, 100.);
  hptgen = fs->make<TH1F>("hptgen" , "pt all", 60, 0, 15.);
  hptreco = fs->make<TH1F>("hptreco" , "pt reconstructed", 60, 0, 15.);
  hrreco = fs->make<TH1F>("hrreco" , "radius of dark pion for reconstructed matched tracks", 100, 0, 100.);
  hdR = fs->make<TH1F>("hdR" , "del R between track and gen", 100, 0,1.);
  hgenrecopt = fs->make<TH2D>("hgenrecopt" , "gen vs reco pt", 100, 0, 10.,100,0.,10.);
  hratvr = fs->make<TH2D>("hratvr" , "reco/gen pt vs r", 100, 0, 100.,100,0.2,2.0);
  hgenreco =   fs->make<TH1F>("hgenreco"," ratio reco to gen pt",100,0.5,1.5);
  hnchgen =   fs->make<TH1F>("hnchgen"," number of charged tracks in decay",20,0.,20.);
  hnchreco =   fs->make<TH1F>("hnchreco"," number of charged tracks in decay",20,0.,20.);
  hegen =   fs->make<TH1F>("hegen","energy of dark pion",30,0.,30.);
  hereco =   fs->make<TH1F>("hereco","energy of dark pion",30,0.,30.);
  hphiphi = fs->make<TH2D>("hphiphi" , "phi phi", 100, -3.2, 3.2,100,-3.2,3.2);
  hphiphip = fs->make<TH2D>("hphiphip" , "phi phi", 100, -3.2, 3.2,100,-3.2,3.2);
  hphiphim = fs->make<TH2D>("hphiphim" , "phi phi", 100, -3.2, 3.2,100,-3.2,3.2);
  hphiphi2 = fs->make<TH2D>("hphiphi2" , "phi phi", 100, -3.2, 3.2,100,-3.2,3.2);
  hphid1 = fs->make<TH2D>("hphid1" , "del phi v rgen", 100, 0.,40.,100,-1,1);
  hphid2 = fs->make<TH2D>("hphid2" , "del phi v rgen", 100, 0.,40.,100,-1,1);
  hii = fs->make<TH2D>("hii" , "impact parameters", 100, -20.,20.,100,-20.,20.);
  hii1 = fs->make<TH2D>("hii1" , "impact parameter versus rgen", 100, 0.,30.,100,-20.,20.);
  hii2 = fs->make<TH2D>("hii2" , "impact parameter versus rgen", 100, 0.,30.,100,-20.,20.);
  himp =   fs->make<TH1F>("himp","track impact parameters",30,-30.,30.);
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
