// -*- C++ -*-
//
// Package:    JPsiKaon
// Class:      JPsiKaon
// 

//=================================================
// Original author:  Jhovanny Andres Mejia        |
//         created:  Fryday Sep 23                |
//         <jhovanny.andres.mejia.guisao@cern.ch> | 
//=================================================

// system include files
#include <memory>


// user include files
#include "myAnalyzers/bph-hlt-tools/src/JPsiKaon.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
//#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"

//For kinematic fit:
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"


#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TH2F.h"

//
// constants, enums and typedefs
//

typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//

//
// constructors and destructor
//

JPsiKaon::JPsiKaon(const edm::ParameterSet& iConfig)
  :
  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  //trakCollection_label(consumes<std::vector<pat::GenericParticle>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  //trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  //trakCollection_label(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("Trak"))),
  //trakCollection_label_lowpt(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak_lowpt"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  triggerCollection_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("TriggerInput"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  algTok_(consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<edm::InputTag>("algInputTag"))),
  algInputTag_(consumes<GlobalAlgBlkBxCollection>(iConfig.getParameter<edm::InputTag>("algInputTag"))),
  l1MuonsToken_(consumes<BXVector<l1t::Muon>>(iConfig.getParameter<edm::InputTag>("l1Muons"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),
  builderToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
  gtUtil_( new l1t::L1TGlobalUtil( iConfig, consumesCollector(), *this, iConfig.getParameter<edm::InputTag>("algInputTag"), iConfig.getParameter<edm::InputTag>("algInputTag"), l1t::UseEventSetupIn::RunAndEvent  )),
  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
  HLTPaths_(iConfig.getParameter<std::vector<std::string>>("HLTPaths")),
  L1Seeds_(iConfig.getParameter<std::vector<std::string>>("L1Seeds")),
  debug_(iConfig.getParameter<bool>("debug")),
  tree_(0), 

  mumC2(0), mumNHits(0), mumNPHits(0),
  mupC2(0), mupNHits(0), mupNPHits(0),
  mumdxy(0), mupdxy(0), mumdz(0), mupdz(0),
  muon_dca(0),

  L10er1p5dR(0), L10er1p4dR(0), L14dR(0), L14p5dR(0),
  L10er2p0dEta1p5(0), L10er2p0dEta1p6(0),
  L10er1p4dEta1p6(0), L13er2p0dR(0),
  L14p5ups(0),

  L1mu_pt(0), L1mu_eta(0), L1mu_phi(0), L1mu_etaAtVtx(0), L1mu_phiAtVtx(0), L1mu_charge(0), L1mu_quality(0),
  hltsVector(HLTPaths_.size()),
  l1sVector(L1Seeds_.size()),

  HLTLowMassDisplaced(0), HLTLowMassInclusive(0),
  HLTMuMuTrkDisplaced(0), HLTBsMuMu(0), HLTUpsilon(0),

  tri_LowMassInclusive(0), tri_LowMassDisplaced(0),
  tri_Dim25(0), tri_Dim20(0), tri_JpsiTk(0),
 
  mu1HLTmatched(0), mu2HLTmatched(0),
 
  mu1pt(-1), mu1eta(-1), mu1phi(-1), mu1charge(-1), mu1dRProbe(-1), mu1dRTag(-1),
  mu2pt(-1), mu2eta(-1), mu2phi(-1), mu2charge(-1), mu2dRProbe(-1), mu2dRTag(-1),

  mu1soft(0), mu2soft(0), mu1tight(0), mu2tight(0), 
  mu1PF(0), mu2PF(0), mu1loose(0), mu2loose(0),
 
  // *******************************************************
 
  nB(0), nMu(0),
  //B_mass(0), B_px(0), B_py(0), B_pz(0), B_charge(0),
  //B_k_px(0), B_k_py(0), B_k_pz(0), B_k_charge1(0),
  //B_k_px_track(0), B_k_py_track(0), B_k_pz_track(0),
  //B_J_mass(0), B_J_px(0), B_J_py(0), B_J_pz(0),

  //B_J_pt1(0), B_J_px1(0), B_J_py1(0), B_J_pz1(0),
  //B_J_pt2(0), B_J_px2(0), B_J_py2(0), B_J_pz2(0), 
  //B_J_charge1(0), B_J_charge2(0),

  J_mass(0), J_px(0), J_py(0), J_pz(0), J_charge(0),
  J_pt1(0), J_eta1(0), J_phi1(0), J_px1(0), J_py1(0), J_pz1(0),
  J_pt2(0), J_eta2(0), J_phi2(0), J_px2(0), J_py2(0), J_pz2(0), 
  J_charge1(0), J_charge2(0),

  // Primary Vertex (PV)
  nVtx(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priVtxXYE(0), priVtxXZE(0), priVtxYZE(0),
  
  // ************************ ****************************************************

  //B_chi2(0), B_J_chi2(0),
  //B_Prob(0), B_J_Prob(0), 
 
  //B_DecayVtxX(0),     B_DecayVtxY(0),     B_DecayVtxZ(0),
  //B_DecayVtxXE(0),    B_DecayVtxYE(0),    B_DecayVtxZE(0),
  //B_DecayVtxXYE(0),   B_DecayVtxXZE(0),   B_DecayVtxYZE(0),

  J_chi2(0),
  J_Prob(0), 
  
  J_DecayVtxX(0), J_DecayVtxY(0), J_DecayVtxZ(0),
  J_DecayVtxXE(0), J_DecayVtxYE(0), J_DecayVtxZE(0),
  J_DecayVtxXYE(0), J_DecayVtxXZE(0), J_DecayVtxYZE(0),
 
  run(0), event(0),
  lumiblock(0)

{
   //now do what ever initialization is needed
}

JPsiKaon::~JPsiKaon()
{

}


// ------------ method called to for each event  ------------
void JPsiKaon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if (debug_) std::cout << "\n\n_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n"<< std::endl;
  if (debug_) std::cout << "Begin analyze!!" << std::endl;
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  //*********************************
  // Get event content information
  //*********************************  
 

  /*if (debug_){
    for(auto &&hlt_path: HLTPaths_)
      std::cout << hlt_path << std::endl;
  }*/
  


  // Kinematic fit
  //edm::ESHandle<TransientTrackBuilder> theB; 
  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 
  auto const &theB = iSetup.getData(builderToken_);

  //edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;
  //iEvent.getByToken(trakCollection_label,thePATTrackHandle);

  edm::Handle< View<pat::Muon> > thePATMuonHandle;
  iEvent.getByToken(dimuon_Label,thePATMuonHandle);

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerResults_Label, triggerBits);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerCollection;
  iEvent.getByToken(triggerCollection_, triggerCollection);
 
  //*********************************
  //Now we get the primary vertex 
  //*********************************
  reco::Vertex bestVtx;
  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  // get primary vertex
  bestVtx = *(primaryVertices_handle->begin());

  priVtxX = bestVtx.x();
  priVtxY = bestVtx.y();
  priVtxZ = bestVtx.z();
  priVtxXE = bestVtx.covariance(0, 0);
  priVtxYE = bestVtx.covariance(1, 1);
  priVtxZE = bestVtx.covariance(2, 2);
  priVtxXYE = bestVtx.covariance(0, 1);
  priVtxXZE = bestVtx.covariance(0, 2);
  priVtxYZE = bestVtx.covariance(1, 2);



  priVtxCL = ChiSquaredProbability((double)(bestVtx.chi2()),(double)(bestVtx.ndof())); 
  nVtx = primaryVertices_handle->size(); 
 
  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event();


  // Here we are going to save the information of the HLT_paths and L1seeds
  hltsVector.resize( HLTPaths_.size(), 0);
  l1sVector.resize( L1Seeds_.size(), 0);

  //Unpack trigger info
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  const pat::TriggerObjectStandAloneCollection unPackedCollection;

  for (unsigned int i = 0; i  < triggerBits->size(); i++) {
    std::string iName = names.triggerName(i);
    if (triggerBits->accept(i)) {
      if (iName.find("HLT_DoubleMu4_LowMass_Displaced_v") != std::string::npos) HLTLowMassDisplaced = 1;
      if (iName.find("HLT_DoubleMu4_3_LowMass_v") != std::string::npos) HLTLowMassInclusive = 1;
      if (iName.find("HLT_DoubleMu4_MuMuTrk_Displaced_v") != std::string::npos) HLTMuMuTrkDisplaced = 1;
      if (iName.find("HLT_DoubleMu4_3_Bs_v") != std::string::npos) HLTBsMuMu = 1;
      if (iName.find("HLT_Dimuon10_y1p4_v") != std::string::npos) HLTUpsilon = 1;
    
      for(std::size_t i = 0; i < HLTPaths_.size(); ++i) {
        if (iName.find(HLTPaths_[i]) != std::string::npos) hltsVector[i] = 1;
      }

    }


  } 

  /*if (debug_){
    for(std::size_t i = 0; i < HLTPaths_.size(); ++i)
      std::cout<<HLTPaths_[i] << "\t--> " << hltsVector[i] << std::endl; 
  }*/

  if (debug_) std::cout << "Len Triggers : " << triggerCollection->size() << std::endl;
  for (pat::TriggerObjectStandAlone trig : *triggerCollection) {
      trig.unpackPathNames(names);
      trig.unpackFilterLabels(iEvent, *triggerBits);
  }

  gtUtil_->retrieveL1(iEvent, iSetup, algInputTag_);
  const vector<pair<string, bool> > finalDecisions = gtUtil_->decisionsFinal();
  for (size_t i_l1t = 0; i_l1t < finalDecisions.size(); i_l1t++){
    string l1tName = (finalDecisions.at(i_l1t)).first;
    if ((finalDecisions.at(i_l1t)).second == 1) {
      if( l1tName.find("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4") != string::npos )    L10er1p5dR = 1;
      if( l1tName.find("L1_DoubleMu4_SQ_OS_dR_Max1p2") != string::npos )         L14dR = 1;
      if( l1tName.find("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4") != string::npos )    L10er1p4dR = 1;
      if( l1tName.find("L1_DoubleMu4p5_SQ_OS_dR_Max1p2") != string::npos )       L14p5dR = 1;
      if( l1tName.find("L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6") != string::npos )  L10er2p0dEta1p6 = 1;
      if( l1tName.find("L1_DoubleMu0er1p4_SQ_OS_dEta_Max1p6") != string::npos )  L10er1p4dEta1p6 = 1;
      if( l1tName.find("L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5") != string::npos )  L10er2p0dEta1p5 = 1;
      if( l1tName.find("L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4") != string::npos )    L13er2p0dR = 1;
      if( l1tName.find("L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18") != string::npos ) L14p5ups = 1;

      for(std::size_t i = 0; i < L1Seeds_.size(); ++i) {
        if (l1tName.find(L1Seeds_[i]) != std::string::npos) l1sVector[i] = 1;
      }

    }
  }

  edm::Handle<BXVector<l1t::Muon> > gmuons;
  iEvent.getByToken(l1MuonsToken_, gmuons);

//  for (int ibx = gmuons->getFirstBX(); ibx <= gmuons->getLastBX(); ++ibx) {

  for (auto itr = gmuons->begin(0); itr != gmuons->end(0); ++itr) {
    L1mu_pt.push_back(itr->pt());
    L1mu_eta.push_back(itr->eta());
    L1mu_phi.push_back(itr->phi());
    L1mu_etaAtVtx.push_back(itr->etaAtVtx());
    L1mu_phiAtVtx.push_back(itr->phiAtVtx());
    L1mu_quality.push_back(itr->hwQual());
    L1mu_charge.push_back(itr->charge());
  }

//  }

  std::string hltMuColl = "hltIterL3MuonCandidates";

  //*****************************************
  //Let's begin by looking for J/psi+K^+

  unsigned int nMu_tmp = thePATMuonHandle->size();
  //nMu = nMu_tmp;

  for(View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) {

    if((iMuon1->track()).isNull()) continue;
    if(iMuon1->track()->pt()<2.0) continue;

    if (debug_) std::cout << "Muon1 track()->eta() =  " << iMuon1->track()->eta() << std::endl;
    //if (debug_) std::cout << "Muon1 eta()          =  " << iMuon1->eta() << std::endl;
    if (debug_) std::cout << "Muon1 track()->phi() =  " << iMuon1->track()->phi() << std::endl;
    //if (debug_) std::cout << "Muon1 phi()          =  " << iMuon1->phi() << std::endl;
    if (debug_) std::cout << "Muon1 track()->pt()  =  " << iMuon1->track()->pt() << std::endl;
    //if (debug_) std::cout << "Muon1 pt()           =  " << iMuon1->pt() << std::endl;
    if (debug_) std::cout << "\n\n" ;

    mu1HLTmatched = false;
    int i = -1;
    for (pat::TriggerObjectStandAlone obj : *triggerCollection) {
      i++;
      if( obj.collection().find(hltMuColl) != std::string::npos ) {
        /* float etaTrig = obj.eta();
        float phiTrig = obj.phi();
        float etaReco = iMuon1->track()->eta();
        float phiReco = iMuon1->track()->phi();
        float de = etaTrig - etaReco;
        auto dp = std::abs(phiTrig - phiReco);
        if (dp > float(M_PI)) dp -= float(2*M_PI); */
        //if (debug_) std::cout << " -> mu1 V " << sqrt(de*de + dp*dp) << std::endl;
        float dr = reco::deltaR(obj.eta(), 
                                  obj.phi(), 
                                  iMuon1->track()->eta(), 
                                  iMuon1->track()->phi());

        if (debug_) std::cout << "\t HLT Object find(hltIterL3MuonCandidates)  " << i << "     dr = " << dr <<std::endl;
        if (debug_) std::cout << "\t\t pt, eta, phi     = ";
        if (debug_) std::cout << obj.pt()  << " " ;
        if (debug_) std::cout << obj.eta() << " " ;
        if (debug_) std::cout << obj.phi() << "\n" ;

        if(debug_){
          obj.unpackPathNames(names);
          std::vector<std::string> paths_matched = obj.pathNames(true, true);
          std::cout << "\t\t Paths : " ;
          for (auto tmp_path: paths_matched){
            std::cout << tmp_path << " ";
          }
          std::cout<<"\n";
          
          obj.unpackFilterLabels(iEvent, *triggerBits);
          std::vector<std::string> conditions = obj.conditionNames();
          std::cout << "\t\t Conditions : " ;
          for (auto tmp_cond: conditions){
            std::cout << tmp_cond << " ";
          }
          std::cout << "\n";
          
          std::vector<std::string> algos_names = obj.algorithmNames(true);
          std::cout << "\t\t Algos : " ;
          for (auto tmp_algo: algos_names){
            std::cout << tmp_algo << " ";
          }
          
          std::cout << "\n\n";
        }

        //if (debug_) std::cout << " -> mu1 dR " << dr << std::endl;
        //if (sqrt(de*de + dp*dp) < 0.01) mu1HLTmatched = true;
        if (dr < 0.01) mu1HLTmatched = true;
      }
    }







    if (debug_) std::cout << "# triggerObjectMatches = " << iMuon1->triggerObjectMatches().size() << std::endl;
    
    std::string TagPath = "HLT_DoubleMu4_3_LowMass_v*";
    std::string ProbePath = "HLT_Mu8_v*";
    
    //char TagPath [ ("HLT_Mu8").size() + 1];
    //strcpy( TagPath, ("HLT_Mu8").c_str() ); 


    // Iterate over the triggerObjectMatches and check if any object matched the path we are looking for
    if(iMuon1->triggerObjectMatches().size()!=0 and false){
      for(size_t i=0; i<iMuon1->triggerObjectMatches().size();i++){
        // hasPathName :
        // // https://cmssdt.cern.ch/lxr/source/DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h#0153
        // // (const std::string &name, bool pathLastFilterAccepted,bool pathL3FilterAccepted))
        //if (iMuon1->triggerObjectMatch(i)!=0 && iMuon1->triggerObjectMatch(i)->hasPathName(TagPath,true,true) ){
          float dr = reco::deltaR(iMuon1->triggerObjectMatch(i)->eta(), 
                                  iMuon1->triggerObjectMatch(i)->phi(), 
                                  iMuon1->eta(), 
                                  iMuon1->phi());
          
          if (debug_) std::cout << "\t HLT Matched Object " << i << "     dr = " << dr <<std::endl;
          if (debug_) std::cout << "\t\t pt, eta, phi     = ";
          if (debug_) std::cout << iMuon1->triggerObjectMatch(i)->pt()  << " " ;
          if (debug_) std::cout << iMuon1->triggerObjectMatch(i)->eta() << " " ;
          if (debug_) std::cout << iMuon1->triggerObjectMatch(i)->phi() << "\n" ;

        //}


        std::vector<std::string> paths_matched = iMuon1->triggerObjectMatch(i)->pathNames(true, true);
        std::cout << "\t\t Paths : " ;
        for (auto tmp_path: paths_matched){
          std::cout << tmp_path << " ";
        }
        std::cout << "\n";

        std::vector<std::string> algos_names = iMuon1->triggerObjectMatch(i)->algorithmNames(true);
        std::cout << "\t\t Algos : " ;
        for (auto tmp_algo: algos_names){
          std::cout << tmp_algo << " ";
        }
        /*std::cout << "\n";

        iMuon1->triggerObjectMatch(i)->unpackFilterLabels(iEvent, *triggerBits);
        std::vector<std::string> conditions = iMuon1->triggerObjectMatch(i)->conditionNames();
        std::cout << "\t\t Conditions : " ;
        for (auto tmp_cond: conditions){
          std::cout << tmp_cond << " ";
        } */
        std::cout << "\n\n";

         

      }
    }
    //const pat::Muon& muon1_obj = *iMuon1;
    float muon1DrTag = matchedDr(iMuon1, TagPath, debug_);
    //float muon1DrProbe = matchedDr(iMuon1, ProbePath, debug_);
    
    //if (debug_) std::cout << "\ndRtag : "<< muon1DrTag << std::endl;






        
    for(View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2)	{

      if(iMuon1==iMuon2) continue;
      if((iMuon2->track()).isNull()) continue;
      if(iMuon2->track()->pt()<2.0) continue;
      
      //opposite charge 
      if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue;

            mu2HLTmatched = false;
            for (pat::TriggerObjectStandAlone obj : *triggerCollection) {
              if( obj.collection().find(hltMuColl) != std::string::npos ) {
              float etaTrig = obj.eta();
              float phiTrig = obj.phi();
              float etaReco = iMuon2->track()->eta();
              float phiReco = iMuon2->track()->phi();
              float de = etaTrig - etaReco;
              auto dp = std::abs(phiTrig - phiReco);
              if (dp > float(M_PI)) dp -= float(2*M_PI);
              //std::cout << "mu2 " << sqrt(de*de + dp*dp) << std::endl;
              if (sqrt(de*de + dp*dp) < 0.01) mu2HLTmatched = true;
              }
            }

      TrackRef glbTrackP;	  
      TrackRef glbTrackM;	  
      
      if(iMuon1->charge() == 1){ glbTrackP = iMuon1->track();}
      if(iMuon1->charge() == -1){ glbTrackM = iMuon1->track();}
      
      if(iMuon2->charge() == 1) { glbTrackP = iMuon2->track();}
      if(iMuon2->charge() == -1){ glbTrackM = iMuon2->track();}
      
      if( glbTrackP.isNull() || glbTrackM.isNull() ) 
        {
          if (debug_) std::cout << "continue due to no track ref" << endl;
          continue;
        }

      if(iMuon1->track()->pt()<2.0) continue;
      if(iMuon2->track()->pt()<2.0) continue;

      if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
      if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;	 
      
      reco::TransientTrack muon1TT(theB.build(glbTrackP));
      reco::TransientTrack muon2TT(theB.build(glbTrackM));

    // *****  Trajectory states to calculate DCA for the 2 muons *********************
      FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
      FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();

      if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() ) continue;

      // Measure distance between tracks at their closest approach
      ClosestApproachInRPhi cApp;
      cApp.calculate(mu1State, mu2State);
      if( !cApp.status() ) continue;
      float dca = fabs( cApp.distance() );	  
      //if (dca < 0. || dca > 0.5) continue;
      //cout<<" closest approach  "<<dca<<endl;

      // *****  end DCA for the 2 muons *********************

      //Let's check the vertex and mass

      //The mass of a muon and the insignificant mass sigma 
      //to avoid singularities in the covariance matrix.
      ParticleMass muon_mass = 0.10565837; //pdg mass
      ParticleMass psi_mass = 3.096916;
      float muon_sigma = muon_mass*1.e-6;
      
      //Creating a KinematicParticleFactory
      KinematicParticleFactoryFromTransientTrack pFactory;
      
      //initial chi2 and ndf before kinematic fits.
      float chi = 0.;
      float ndf = 0.;
      vector<RefCountedKinematicParticle> muonParticles;
      try {
        muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
        muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
      }
      catch(...) { 
        if (debug_) std::cout<<" Exception caught ... continuing 1 "<<std::endl; 
        continue;
      }
      
      KinematicParticleVertexFitter fitter;   
      
      RefCountedKinematicTree psiVertexFitTree;
      try {
        psiVertexFitTree = fitter.fit(muonParticles); 
      }
      catch (...) { 
        std::cout<<" Exception caught ... continuing 2 "<<std::endl; 
        continue;
      }
      
      if (!psiVertexFitTree->isValid()) 
        {
          //std::cout << "caught an exception in the psi vertex fit" << std::endl;
          continue; 
        }
      
      psiVertexFitTree->movePointerToTheTop();
      
      RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();
      RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();
      
      if( psi_vFit_vertex_noMC->chiSquared() < 0 ) continue;
      
      //some loose cuts go here
      
      if(psi_vFit_vertex_noMC->chiSquared()>50.) continue;
      if(psi_vFit_noMC->currentState().mass()<0.05) continue;

      double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
      if(J_Prob_tmp<0.01) continue;

      psiVertexFitTree->movePointerToTheFirstChild();
      RefCountedKinematicParticle mu1PsiCand = psiVertexFitTree->currentParticle();
                
      psiVertexFitTree->movePointerToTheNextChild();
      RefCountedKinematicParticle mu2PsiCand = psiVertexFitTree->currentParticle();

      J_mass = psi_vFit_noMC->currentState().mass();
      J_px = psi_vFit_noMC->currentState().globalMomentum().x();
      J_py = psi_vFit_noMC->currentState().globalMomentum().y();
      J_pz = psi_vFit_noMC->currentState().globalMomentum().z();
      J_charge = psi_vFit_noMC->currentState().particleCharge();

      J_pt1 = mu1PsiCand->currentState().globalMomentum().perp();
      J_eta1 = mu1PsiCand->currentState().globalMomentum().eta();
      J_phi1 = mu1PsiCand->currentState().globalMomentum().phi();
      J_px1 = mu1PsiCand->currentState().globalMomentum().x();
      J_py1 = mu1PsiCand->currentState().globalMomentum().y();
      J_pz1 = mu1PsiCand->currentState().globalMomentum().z();
      J_charge1 = mu1PsiCand->currentState().particleCharge();

      J_pt2 = mu2PsiCand->currentState().globalMomentum().perp();
      J_eta2 = mu2PsiCand->currentState().globalMomentum().eta();
      J_phi2 = mu2PsiCand->currentState().globalMomentum().phi();
      J_px2 = mu2PsiCand->currentState().globalMomentum().x();
      J_py2 = mu2PsiCand->currentState().globalMomentum().y();
      J_pz2 = mu2PsiCand->currentState().globalMomentum().z();
      J_charge2 = mu2PsiCand->currentState().particleCharge();

      J_chi2 = psi_vFit_vertex_noMC->chiSquared();
      J_Prob = J_Prob_tmp;

      J_DecayVtxX = (*psi_vFit_vertex_noMC).position().x();    
      J_DecayVtxY = (*psi_vFit_vertex_noMC).position().y();
      J_DecayVtxZ = (*psi_vFit_vertex_noMC).position().z();
      J_DecayVtxXE = (psi_vFit_vertex_noMC->error().cxx());   
      J_DecayVtxYE = (psi_vFit_vertex_noMC->error().cyy());   
      J_DecayVtxZE = (psi_vFit_vertex_noMC->error().czz());
      J_DecayVtxXYE = (psi_vFit_vertex_noMC->error().cyx());
      J_DecayVtxXZE = (psi_vFit_vertex_noMC->error().czx());
      J_DecayVtxYZE = (psi_vFit_vertex_noMC->error().czy());

      //Now that we have a J/psi candidate, we look for K^+ candidates
      /*	  
      for(View<pat::PackedCandidate>::const_iterator iTrack1 = thePATTrackHandle->begin(); iTrack1 != thePATTrackHandle->end(); ++iTrack1 ) {
        if(iTrack1->charge()==0) continue;
        if(fabs(iTrack1->pdgId())!=211) continue;
        if(iTrack1->pt()<1.3) continue;
        //if(iTrack1->pt()<0.95) continue;
        if(!(iTrack1->trackHighPurity())) continue;
        
        if ( IsTheSame(*iTrack1,*iMuon1) || IsTheSame(*iTrack1,*iMuon2) ) continue;
                    
        reco::TransientTrack kaonTT((*theB).build(iTrack1->pseudoTrack()));

        ParticleMass kaon_mass = 0.493677;
        float kaon_sigma = kaon_mass*1.e-6;

        float chi = 0.;
        float ndf = 0.;		 

        // ***************************
        // JpsiKaon invariant mass (before kinematic vertex fit)
        // ***************************
        TLorentzVector kaon14V, Jpsi4V; 
        kaon14V.SetXYZM(iTrack1->px(),iTrack1->py(),iTrack1->pz(),kaon_mass);
        
        Jpsi4V.SetXYZM(psi_vFit_noMC->currentState().globalMomentum().x(),psi_vFit_noMC->currentState().globalMomentum().y(),psi_vFit_noMC->currentState().globalMomentum().z(),psi_vFit_noMC->currentState().mass());

        if ( (kaon14V + Jpsi4V).M()<4.2 || (kaon14V + Jpsi4V).M()>6.8 ) continue;
      
        //Now we are ready to combine!
        // JPsi mass constraint is applied in the final Bplus fit,

        vector<RefCountedKinematicParticle> vFitMCParticles;
        vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
        vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
        vFitMCParticles.push_back(pFactory.particle(kaonTT,kaon_mass ,chi,ndf,kaon_sigma));
              
        MultiTrackKinematicConstraint *  j_psi_c = new  TwoTrackMassKinematicConstraint(psi_mass);
        KinematicConstrainedVertexFitter kcvFitter;
        RefCountedKinematicTree vertexFitTree = kcvFitter.fit(vFitMCParticles, j_psi_c);
        if (!vertexFitTree->isValid()) {
          //std::cout << "caught an exception in the B vertex fit with MC" << std::endl;
          continue;
        }
        vertexFitTree->movePointerToTheTop();

      
        RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
        RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
        if (!bDecayVertexMC->vertexIsValid()){
          // cout << "B MC fit vertex is not valid" << endl;
          continue;
        }
          
        if ( (bCandMC->currentState().mass() < 5.0) || (bCandMC->currentState().mass() > 6.0) ) {
          continue;
        }
          
        if ( bDecayVertexMC->chiSquared()<0 || bDecayVertexMC->chiSquared()>50) {
          //if ( bDecayVertexMC->chiSquared()<0 ) cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
          continue;
        }
          
        double B_Prob_tmp  = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
        if(B_Prob_tmp<0.01) {
          continue;
              }
                  
              // get children from final B fit

        vertexFitTree->movePointerToTheFirstChild();
        RefCountedKinematicParticle mu1CandMC = vertexFitTree->currentParticle();
          
        vertexFitTree->movePointerToTheNextChild();
        RefCountedKinematicParticle mu2CandMC = vertexFitTree->currentParticle();
        
        vertexFitTree->movePointerToTheNextChild();
        RefCountedKinematicParticle kCandMC = vertexFitTree->currentParticle();		  

        KinematicParameters psiMu1KP = mu1CandMC->currentState().kinematicParameters();
        KinematicParameters psiMu2KP = mu2CandMC->currentState().kinematicParameters();
        KinematicParameters psiMupKP;
        KinematicParameters psiMumKP;
          
              if ( mu1CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu1KP;
        if ( mu1CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu1KP;
        if ( mu2CandMC->currentState().particleCharge() > 0 ) psiMupKP = psiMu2KP;
        if ( mu2CandMC->currentState().particleCharge() < 0 ) psiMumKP = psiMu2KP;	 

        GlobalVector Jp1vec(mu1CandMC->currentState().globalMomentum().x(),
                mu1CandMC->currentState().globalMomentum().y(),
          mu1CandMC->currentState().globalMomentum().z());


              GlobalVector Jp2vec(mu2CandMC->currentState().globalMomentum().x(),
                mu2CandMC->currentState().globalMomentum().y(),
          mu2CandMC->currentState().globalMomentum().z());
        
              KinematicParameters VCandKP = kCandMC->currentState().kinematicParameters();
                
        // ************ fill candidate variables now
        
        // Only save the first time
        if(nB==0){	    
          nMu  = nMu_tmp;
          // cout<< "*Number of Muons : " << nMu_tmp << endl;
        } // end nB==0

        B_mass->push_back(bCandMC->currentState().mass());
        B_px->push_back(bCandMC->currentState().globalMomentum().x());
        B_py->push_back(bCandMC->currentState().globalMomentum().y());
        B_pz->push_back(bCandMC->currentState().globalMomentum().z());
        B_charge->push_back(bCandMC->currentState().particleCharge());

        // You can get the momentum components (for muons and kaon) from the final B childrens or of the original Tracks. Here, a example for the kaon:
        B_k_px->push_back(VCandKP.momentum().x() );
        B_k_py->push_back(VCandKP.momentum().y() );
        B_k_pz->push_back(VCandKP.momentum().z() );
        B_k_px_track->push_back(iTrack1->px() );
        B_k_py_track->push_back(iTrack1->py() );
        B_k_pz_track->push_back(iTrack1->pz() );
        B_k_charge1->push_back(kCandMC->currentState().particleCharge());
        
        B_J_mass->push_back( psi_vFit_noMC->currentState().mass() );
        B_J_px->push_back( psi_vFit_noMC->currentState().globalMomentum().x() );
        B_J_py->push_back( psi_vFit_noMC->currentState().globalMomentum().y() );
        B_J_pz->push_back( psi_vFit_noMC->currentState().globalMomentum().z() );
          B_J_pt1->push_back(Jp1vec.perp());
              B_J_px1->push_back(psiMu1KP.momentum().x());
              B_J_py1->push_back(psiMu1KP.momentum().y());
        B_J_pz1->push_back(psiMu1KP.momentum().z());
        B_J_charge1->push_back(mu1CandMC->currentState().particleCharge());

        B_J_pt2->push_back(Jp2vec.perp());
        B_J_px2->push_back(psiMu2KP.momentum().x());
        B_J_py2->push_back(psiMu2KP.momentum().y());
        B_J_pz2->push_back(psiMu2KP.momentum().z());
        B_J_charge2->push_back(mu2CandMC->currentState().particleCharge());
        
        B_J_chi2->push_back(psi_vFit_vertex_noMC->chiSquared());
        B_chi2->push_back(bDecayVertexMC->chiSquared());
              
        //double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
        //double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
        B_Prob    ->push_back(B_Prob_tmp);
        B_J_Prob  ->push_back(J_Prob_tmp);

        B_DecayVtxX ->push_back((*bDecayVertexMC).position().x());    
        B_DecayVtxY ->push_back((*bDecayVertexMC).position().y());
        B_DecayVtxZ ->push_back((*bDecayVertexMC).position().z());
        B_DecayVtxXE ->push_back(bDecayVertexMC->error().cxx());   
        B_DecayVtxYE ->push_back(bDecayVertexMC->error().cyy());   
        B_DecayVtxZE ->push_back(bDecayVertexMC->error().czz());
        B_DecayVtxXYE ->push_back(bDecayVertexMC->error().cyx());
        B_DecayVtxXZE ->push_back(bDecayVertexMC->error().czx());
        B_DecayVtxYZE ->push_back(bDecayVertexMC->error().czy());
      
  */
        mu1pt     = iMuon1->pt();
        mu1eta    = iMuon1->eta();
        mu1phi    = iMuon1->phi();
        mu1charge = iMuon1->charge();
        mu1dRProbe= matchedDr(iMuon1, TagPath, false);
        mu1dRTag  = matchedDr(iMuon1, ProbePath, false);

        mu2pt     = iMuon2->pt();
        mu2eta    = iMuon2->eta();
        mu2phi    = iMuon2->phi();
        mu2charge = iMuon2->charge();
        mu2dRProbe= matchedDr(iMuon2, TagPath, false);
        mu2dRTag  = matchedDr(iMuon2, ProbePath, false);

        mu1soft = iMuon1->isSoftMuon(bestVtx);
        mu2soft = iMuon2->isSoftMuon(bestVtx);
        mu1tight = iMuon1->isTightMuon(bestVtx);
        mu2tight = iMuon2->isTightMuon(bestVtx);
        mu1PF = iMuon1->isPFMuon();
        mu2PF = iMuon2->isPFMuon();
        mu1loose = muon::isLooseMuon(*iMuon1);
        mu2loose = muon::isLooseMuon(*iMuon2);
        mumC2 = glbTrackP->normalizedChi2();
        mumNHits = glbTrackP->numberOfValidHits();
        mumNPHits = glbTrackP->hitPattern().numberOfValidPixelHits();	       
        mupC2 = glbTrackM->normalizedChi2();
        mupNHits = glbTrackM->numberOfValidHits();
        mupNPHits = glbTrackM->hitPattern().numberOfValidPixelHits();
        mumdxy = glbTrackP->dxy(bestVtx.position());// 
        mupdxy = glbTrackM->dxy(bestVtx.position());// 
        mumdz = glbTrackP->dz(bestVtx.position());
        mupdz = glbTrackM->dz(bestVtx.position());
        muon_dca = dca;

        nB++;	       
        muonParticles.clear();
        // vFitMCParticles.clear();
        if (debug_) std::cout << "Before tree_->Fill()   " << nB <<  std::endl;
        if (nB > 0) tree_->Fill();


        // Initialize the objects
        J_charge = 0;
        J_mass = 0;    J_px = 0;    J_py = 0;    J_pz = 0;

        J_pt1 = 0;  J_eta1 = 0.; J_phi1 = 0.; J_px1 = 0;  J_py1 = 0;  J_pz1 = 0; J_charge1 = 0;
        J_pt2 = 0;  J_eta2 = 0.; J_phi2 = 0.; J_px2 = 0;  J_py2 = 0;  J_pz2 = 0; J_charge2 = 0;

        J_chi2 = 0;
        J_Prob = 0;

        J_DecayVtxX = 0;     J_DecayVtxY = 0;     J_DecayVtxZ = 0;
        J_DecayVtxXE = 0;    J_DecayVtxYE = 0;    J_DecayVtxZE = 0;
        J_DecayVtxXYE = 0;   J_DecayVtxXZE = 0;   J_DecayVtxYZE = 0;

        nVtx = 0;
        priVtxX = 0;     priVtxY = 0;     priVtxZ = 0;
        priVtxXE = 0;    priVtxYE = 0;    priVtxZE = 0; priVtxCL = 0;
        priVtxXYE = 0;   priVtxXZE = 0;   priVtxYZE = 0;

        mumC2 = 0;
        mumNHits = 0; mumNPHits = 0;
        mupC2 = 0;
        mupNHits = 0; mupNPHits = 0;
        mumdxy = 0; mupdxy = 0; mumdz = 0; mupdz = 0; muon_dca = 0;

        mu1pt = -1, mu1eta = -1, mu1phi = -1, mu1charge = -1, mu1dRProbe = -1, mu1dRTag = -1;
        mu2pt = -1, mu2eta = -1, mu2phi = -1, mu2charge = -1, mu2dRProbe = -1, mu2dRTag = -1;

        mu1soft = 0; mu2soft = 0; mu1tight = 0; mu2tight = 0;
        mu1PF = 0; mu2PF = 0; mu1loose = 0; mu2loose = 0;

        L10er1p5dR = 0; L10er1p4dR = 0; L14dR = 0; L14p5dR = 0;
        L10er2p0dEta1p5 = 0; L10er2p0dEta1p6 = 0;
        L10er1p4dEta1p6 = 0; L13er2p0dR = 0;
        L14p5ups = 0;
          

        if (debug_) std::cout << "Cleaning vectors 1 --> " << std::endl;
        L1mu_pt.clear(); L1mu_eta.clear(); L1mu_phi.clear();
        L1mu_etaAtVtx.clear(); L1mu_phiAtVtx.clear();
        L1mu_quality.clear(); L1mu_charge.clear();
        
        for(std::size_t i = 0; i < HLTPaths_.size(); ++i){
          hltsVector[i] = 0;
        }

        for(std::size_t i = 0; i < L1Seeds_.size(); ++i){
          l1sVector[i] = 0;
        }

        HLTLowMassDisplaced = 0; HLTLowMassInclusive = 0;
        HLTMuMuTrkDisplaced = 0; HLTBsMuMu = 0; HLTUpsilon = 0;


        //   tri_LowMassInclusive = 0; tri_LowMassDisplaced = 0;
        mu1HLTmatched = 0; mu2HLTmatched = 0;

  //	    }//track
	  }//mu2
  }//mu1
 
   
   nB = 0; nMu = 0;
   tri_LowMassInclusive = 0; tri_LowMassDisplaced = 0;

   L10er1p5dR = 0; L10er1p4dR = 0; L14dR = 0; L14p5dR = 0;
   L10er2p0dEta1p5 = 0; L10er2p0dEta1p6 = 0;
   L10er1p4dEta1p6 = 0; L13er2p0dR = 0;
   L14p5ups = 0;

   if (debug_) std::cout << "Cleaning vectors 2--> " << std::endl;
   L1mu_pt.clear(); L1mu_eta.clear(); L1mu_phi.clear();
   L1mu_etaAtVtx.clear(); L1mu_phiAtVtx.clear();
   L1mu_quality.clear(); L1mu_charge.clear();

   for(std::size_t i = 0; i < HLTPaths_.size(); ++i){
     hltsVector[i] = 0;
   }

   for(std::size_t i = 0; i < L1Seeds_.size(); ++i){
     l1sVector[i] = 0;
   }

   HLTLowMassDisplaced = 0; HLTLowMassInclusive = 0;
   HLTMuMuTrkDisplaced = 0; HLTBsMuMu = 0; HLTUpsilon = 0;


  /*
   B_charge->clear();
   B_mass->clear();    B_px->clear();    B_py->clear();    B_pz->clear(); 
   B_k_px->clear(); B_k_py->clear(); B_k_pz->clear();  B_k_charge1->clear();
   B_k_px_track->clear(); B_k_py_track->clear(); B_k_pz_track->clear();

   B_J_mass->clear();  B_J_px->clear();  B_J_py->clear();  B_J_pz->clear();
 
   B_J_pt1->clear();  B_J_px1->clear();  B_J_py1->clear();  B_J_pz1->clear(), B_J_charge1->clear();
   B_J_pt2->clear();  B_J_px2->clear();  B_J_py2->clear();  B_J_pz2->clear(), B_J_charge2->clear();

   B_chi2->clear(); B_J_chi2->clear(); 
   B_Prob->clear(); B_J_Prob->clear();

   B_DecayVtxX->clear();     B_DecayVtxY->clear();     B_DecayVtxZ->clear();
   B_DecayVtxXE->clear();    B_DecayVtxYE->clear();    B_DecayVtxZE->clear();
   B_DecayVtxXYE->clear();   B_DecayVtxXZE->clear();   B_DecayVtxYZE->clear();
  */
  /*
   J_charge->clear();
   J_mass->clear();    J_px->clear();    J_py->clear();    J_pz->clear(); 
 
   J_pt1->clear();  J_px1->clear();  J_py1->clear();  J_pz1->clear(), J_charge1->clear();
   J_pt2->clear();  J_px2->clear();  J_py2->clear();  J_pz2->clear(), J_charge2->clear();

   J_chi2->clear(); 
   J_Prob->clear();

   J_DecayVtxX->clear();     J_DecayVtxY->clear();     J_DecayVtxZ->clear();
   J_DecayVtxXE->clear();    J_DecayVtxYE->clear();    J_DecayVtxZE->clear();
   J_DecayVtxXYE->clear();   J_DecayVtxXZE->clear();   J_DecayVtxYZE->clear();

   nVtx = 0;
   priVtxX = 0;     priVtxY = 0;     priVtxZ = 0; 
   priVtxXE = 0;    priVtxYE = 0;    priVtxZE = 0; priVtxCL = 0;
   priVtxXYE = 0;   priVtxXZE = 0;   priVtxYZE = 0;    

   mumC2->clear();
   mumNHits->clear(); mumNPHits->clear();
   mupC2->clear();
   mupNHits->clear(); mupNPHits->clear();
   mumdxy->clear(); mupdxy->clear(); mumdz->clear(); mupdz->clear(); muon_dca->clear();

   tri_Dim25->clear(); tri_Dim20->clear(); tri_JpsiTk->clear();
 
   mu1soft->clear(); mu2soft->clear(); mu1tight->clear(); mu2tight->clear();
   mu1PF->clear(); mu2PF->clear(); mu1loose->clear(); mu2loose->clear(); 
  */ 
  if (debug_) std::cout << "End Analyze --> " << std::endl;
}





// // Checks if a given muon has a match with the HLTpath, evaluates its dR and returns it
// // If no match is found it returns -1
// // //The match is done with all `triggerObjectMatches` from the muon and must pass:
// // // // hasPathName(HLTpath,true,true)
float JPsiKaon::matchedDr(edm::View<pat::Muon>::const_iterator muon, const std::string HLTpath, bool debug_) {

    if (debug_) std::cout << "Hello from JPsiKaon::matchedDr !\n";
    // Returned dr
    float dr=-1;
    float dr_tmp;

    if (debug_) std::cout << "# triggerObjectMatches = " << muon->triggerObjectMatches().size() << std::endl;
    //std::string HLTpath = "HLT_Mu8";

    // Iterate over the triggerObjectMatches and check if any object matched the path we are looking for
    if(muon->triggerObjectMatches().size()!=0){
      for(size_t i=0; i<muon->triggerObjectMatches().size();i++){
        // hasPathName :
        // // https://cmssdt.cern.ch/lxr/source/DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h#0153
        // // (const std::string &name, bool pathLastFilterAccepted,bool pathL3FilterAccepted))
        if (muon->triggerObjectMatch(i)!=0 && muon->triggerObjectMatch(i)->hasPathName(HLTpath,true,true) ){
          dr_tmp = reco::deltaR(muon->triggerObjectMatch(i)->eta(), 
                                  muon->triggerObjectMatch(i)->phi(), 
                                  muon->eta(), 
                                  muon->phi());
          //Update in case there are multiple matches!
          if(dr<0 || dr_tmp<dr){
             if( dr>0 && debug_ ) std::cout << "  x- x- x- x- x- x- Updated dr! " << dr << "  :  "  << dr_tmp << std::endl;
            dr = dr_tmp;
          }
          if (debug_) std::cout << "\t HLT Matched Object " << i << "     dr = " << dr <<std::endl;
          if (debug_) std::cout << "\t\t pt, eta, phi     = ";
          if (debug_) std::cout << muon->triggerObjectMatch(i)->pt()  << " " ;
          if (debug_) std::cout << muon->triggerObjectMatch(i)->eta() << " " ;
          if (debug_) std::cout << muon->triggerObjectMatch(i)->phi() << "\n" ;

          // pathNames
          if (debug_){
            std::vector<std::string> paths_matched = muon->triggerObjectMatch(i)->pathNames(true, true);
            std::cout << "\t\t Paths : " ;
            for (auto tmp_path: paths_matched){
              std::cout << tmp_path << " ";
            }
            std::cout << "\n";

            std::vector<std::string> algos_names = muon->triggerObjectMatch(i)->algorithmNames(true);
            std::cout << "\t\t Algos : " ;
            for (auto tmp_algo: algos_names){
              std::cout << tmp_algo << " ";
            }
          }
        }

      }
    }
    return dr;
}

bool JPsiKaon::IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

// ------------ method called once each job just before starting event loop  ------------

void 
JPsiKaon::beginJob()
{

  std::cout << "Beginning analyzer job with value of isMC= " << isMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","Dimuons ntuple");

  tree_->Branch("nB",&nB,"nB/i");
  tree_->Branch("nMu",&nMu,"nMu/i");
  /*
  tree_->Branch("B_charge", &B_charge);
  tree_->Branch("B_mass", &B_mass);
  tree_->Branch("B_px", &B_px);
  tree_->Branch("B_py", &B_py);
  tree_->Branch("B_pz", &B_pz);

  tree_->Branch("B_k_charge1", &B_k_charge1);
  tree_->Branch("B_k_px", &B_k_px);
  tree_->Branch("B_k_py", &B_k_py);
  tree_->Branch("B_k_pz", &B_k_pz);
  tree_->Branch("B_k_px_track", &B_k_px_track);
  tree_->Branch("B_k_py_track", &B_k_py_track);
  tree_->Branch("B_k_pz_track", &B_k_pz_track);

  tree_->Branch("B_J_mass", &B_J_mass);
  tree_->Branch("B_J_px", &B_J_px);
  tree_->Branch("B_J_py", &B_J_py);
  tree_->Branch("B_J_pz", &B_J_pz);

  tree_->Branch("B_J_pt1", &B_J_pt1);
  tree_->Branch("B_J_px1", &B_J_px1);
  tree_->Branch("B_J_py1", &B_J_py1);
  tree_->Branch("B_J_pz1", &B_J_pz1);
  tree_->Branch("B_J_charge1", &B_J_charge1);

  tree_->Branch("B_J_pt2", &B_J_pt2);
  tree_->Branch("B_J_px2", &B_J_px2);
  tree_->Branch("B_J_py2", &B_J_py2);
  tree_->Branch("B_J_pz2", &B_J_pz2);
  tree_->Branch("B_J_charge2", &B_J_charge2);

  tree_->Branch("B_chi2",    &B_chi2);
  tree_->Branch("B_J_chi2",  &B_J_chi2);

  tree_->Branch("B_Prob",    &B_Prob);
  tree_->Branch("B_J_Prob",  &B_J_Prob);
       
  tree_->Branch("B_DecayVtxX",     &B_DecayVtxX);
  tree_->Branch("B_DecayVtxY",     &B_DecayVtxY);
  tree_->Branch("B_DecayVtxZ",     &B_DecayVtxZ);
  tree_->Branch("B_DecayVtxXE",    &B_DecayVtxXE);
  tree_->Branch("B_DecayVtxYE",    &B_DecayVtxYE);
  tree_->Branch("B_DecayVtxZE",    &B_DecayVtxZE);
  tree_->Branch("B_DecayVtxXYE",   &B_DecayVtxXYE);
  tree_->Branch("B_DecayVtxXZE",   &B_DecayVtxXZE);
  tree_->Branch("B_DecayVtxYZE",   &B_DecayVtxYZE);
  */

  tree_->Branch("J_charge", &J_charge);
  tree_->Branch("J_mass", &J_mass);
  tree_->Branch("J_px", &J_px);
  tree_->Branch("J_py", &J_py);
  tree_->Branch("J_pz", &J_pz);

  tree_->Branch("J_pt1", &J_pt1);
  tree_->Branch("J_eta1", &J_eta1);
  tree_->Branch("J_phi1", &J_phi1);
  tree_->Branch("J_px1", &J_px1);
  tree_->Branch("J_py1", &J_py1);
  tree_->Branch("J_pz1", &J_pz1);
  tree_->Branch("J_charge1", &J_charge1);

  tree_->Branch("J_pt2", &J_pt2);
  tree_->Branch("J_eta2", &J_eta2);
  tree_->Branch("J_phi2", &J_phi2);
  tree_->Branch("J_px2", &J_px2);
  tree_->Branch("J_py2", &J_py2);
  tree_->Branch("J_pz2", &J_pz2);
  tree_->Branch("J_charge2", &J_charge2);

  tree_->Branch("J_chi2",  &J_chi2);
  tree_->Branch("J_Prob",  &J_Prob);
       
  tree_->Branch("J_DecayVtxX",     &J_DecayVtxX);
  tree_->Branch("J_DecayVtxY",     &J_DecayVtxY);
  tree_->Branch("J_DecayVtxZ",     &J_DecayVtxZ);
  tree_->Branch("J_DecayVtxXE",    &J_DecayVtxXE);
  tree_->Branch("J_DecayVtxYE",    &J_DecayVtxYE);
  tree_->Branch("J_DecayVtxZE",    &J_DecayVtxZE);
  tree_->Branch("J_DecayVtxXYE",   &J_DecayVtxXYE);
  tree_->Branch("J_DecayVtxXZE",   &J_DecayVtxXZE);
  tree_->Branch("J_DecayVtxYZE",   &J_DecayVtxYZE);

  tree_->Branch("priVtxX",&priVtxX, "priVtxX/f");
  tree_->Branch("priVtxY",&priVtxY, "priVtxY/f");
  tree_->Branch("priVtxZ",&priVtxZ, "priVtxZ/f");
  tree_->Branch("priVtxXE",&priVtxXE, "priVtxXE/f");
  tree_->Branch("priVtxYE",&priVtxYE, "priVtxYE/f");
  tree_->Branch("priVtxZE",&priVtxZE, "priVtxZE/f");
  tree_->Branch("priVtxXYE",&priVtxXYE, "priVtxXYE/f");
  tree_->Branch("priVtxXZE",&priVtxXZE, "priVtxXZE/f");
  tree_->Branch("priVtxYZE",&priVtxYZE, "priVtxYZE/f");
  tree_->Branch("priVtxCL",&priVtxCL, "priVtxCL/f");

  tree_->Branch("nVtx",       &nVtx);
  tree_->Branch("run",        &run,       "run/I");
  tree_->Branch("event",        &event,     "event/I");
  tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");
    
  // *************************
 
  tree_->Branch("mumC2",&mumC2);  
  tree_->Branch("mumNHits",&mumNHits);
  tree_->Branch("mumNPHits",&mumNPHits);
  tree_->Branch("mupC2",&mupC2);  
  tree_->Branch("mupNHits",&mupNHits);
  tree_->Branch("mupNPHits",&mupNPHits);
  tree_->Branch("mumdxy",&mumdxy);
  tree_->Branch("mupdxy",&mupdxy);
  tree_->Branch("mumdz",&mumdz);
  tree_->Branch("mupdz",&mupdz);
  tree_->Branch("muon_dca",&muon_dca);

  tree_->Branch("L10er1p5dR", &L10er1p5dR);
  tree_->Branch("L14dR", &L14dR);
  tree_->Branch("L10er1p4dR", &L10er1p4dR);
  tree_->Branch("L14p5dR", &L14p5dR);
  tree_->Branch("L10er2p0dEta1p5", &L10er2p0dEta1p5);
  tree_->Branch("L10er2p0dEta1p6", &L10er2p0dEta1p6);
  tree_->Branch("L10er1p4dEta1p6", &L10er1p4dEta1p6);
  tree_->Branch("L13er2p0dR", &L13er2p0dR);
  tree_->Branch("L14p5ups", &L14p5ups);

  tree_->Branch("L1mu_pt", &L1mu_pt); 
  tree_->Branch("L1mu_eta", &L1mu_eta);
  tree_->Branch("L1mu_phi", &L1mu_phi);
  tree_->Branch("L1mu_etaAtVtx", &L1mu_etaAtVtx);
  tree_->Branch("L1mu_phiAtVtx", &L1mu_phiAtVtx);
  tree_->Branch("L1mu_charge", &L1mu_charge);
  tree_->Branch("L1mu_quality", &L1mu_quality);

  tree_->Branch("HLTLowMassInclusive", &HLTLowMassInclusive);
  tree_->Branch("HLTLowMassDisplaced", &HLTLowMassDisplaced);
  tree_->Branch("HLTUpsilon", &HLTUpsilon);
  tree_->Branch("HLTMuMuTrkDisplaced", &HLTMuMuTrkDisplaced);
  tree_->Branch("HLTBsMuMu", &HLTBsMuMu);

  if (debug_) std::cout << "Filling Branches" << std::endl;
  for(std::size_t i = 0; i < HLTPaths_.size(); ++i){ 
    tree_->Branch(HLTPaths_[i].c_str(), &hltsVector[i]);
  }
  for(std::size_t i = 0; i < L1Seeds_.size(); ++i){ 
    tree_->Branch(L1Seeds_[i].c_str(), &l1sVector[i]);
  }

  if (debug_) std::cout << "Branches filled :)" << std::endl;


              

  tree_->Branch("tri_LowMassInclusive", &tri_LowMassInclusive);
  tree_->Branch("tri_LowMassDisplaced", &tri_LowMassDisplaced);

  tree_->Branch("mu1HLTmatched", &mu1HLTmatched);
  tree_->Branch("mu2HLTmatched", &mu2HLTmatched);

  tree_->Branch("tri_Dim25",&tri_Dim25);
  tree_->Branch("tri_Dim20",&tri_Dim20);
  tree_->Branch("tri_JpsiTk",&tri_JpsiTk);

  tree_->Branch("mu1pt", &mu1pt);
  tree_->Branch("mu1eta", &mu1eta);
  tree_->Branch("mu1phi", &mu1phi);
  tree_->Branch("mu1charge", &mu1charge);
  tree_->Branch("mu1dRProbe", &mu1dRProbe);
  tree_->Branch("mu1dRTag", &mu1dRTag);

  tree_->Branch("mu2pt", &mu2pt);
  tree_->Branch("mu2eta", &mu2eta);
  tree_->Branch("mu2phi", &mu2phi);
  tree_->Branch("mu2charge", &mu2charge);
  tree_->Branch("mu2dRProbe", &mu2dRProbe);
  tree_->Branch("mu2dRTag", &mu2dRTag);

  tree_->Branch("mu1soft",&mu1soft);
  tree_->Branch("mu2soft",&mu2soft);
  tree_->Branch("mu1tight",&mu1tight);
  tree_->Branch("mu2tight",&mu2tight);
  tree_->Branch("mu1PF",&mu1PF);
  tree_->Branch("mu2PF",&mu2PF);
  tree_->Branch("mu1loose",&mu1loose);
  tree_->Branch("mu2loose",&mu2loose);

}

// ------------ method called once each job just after ending the event loop  ------------
void JPsiKaon::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiKaon);
