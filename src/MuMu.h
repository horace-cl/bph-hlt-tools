#ifndef _MuMu_h
#define _MuMu_h

// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

//For kinematic fit:
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
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


#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"





#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
//
// class decleration
//

//class MuMu : public edm::EDAnalyzer {
class MuMu : public edm::one::EDAnalyzer<> {  
public:
  explicit MuMu(const edm::ParameterSet&);
  //~MuMu();
  ~MuMu() override = default;
  void fillPsi(const reco::Candidate& genpsi);
  void fillV0(const reco::Candidate& genv0);
  //int const getMuCat(reco::Muon const& muon) const;
  bool IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu);
  bool   isAncestor(const reco::Candidate*, const reco::Candidate*);
  double GetLifetime(TLorentzVector, TVector3, TVector3);
  //bool checkDeltaR(double phi1, double eta1, const std::vector<float>& phiVec, const std::vector<float>& etaVec, double threshold)
  //int MatchedObjs(TLorentzVector, TVector3, TVector3);
  
private:

  HLTPrescaleProvider hltPrescale_;
  //HLTConfigProvider   hltConfig_;


  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void printout(const RefCountedKinematicVertex& myVertex) const;
  void printout(const RefCountedKinematicParticle& myParticle) const;
  void printout(const RefCountedKinematicTree& myTree) const;



  // ----------member data ---------------------------
  //const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttrkToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttrkToken_;
  edm::EDGetTokenT<edm::View<pat::Muon>> muon_Label;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trakCollection_label;

  edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenToken_;

  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  edm::EDGetTokenT<reco::BeamSpot> BSLabel_;

  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerCollection_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
  edm::EDGetToken algTok_;
  edm::EDGetTokenT<GlobalAlgBlkBxCollection> algInputTag_;
  edm::EDGetTokenT<BXVector<l1t::Muon>> l1MuonsToken_;
  //edm::EDGetTokenT<GlobalAlgBlkBxCollection> algInputTag_;
  std::vector<std::string> HLTPaths_;
  std::vector<std::string> L1Seeds_;
  l1t::L1TGlobalUtil* gtUtil_;

  bool OnlyBest_;
  bool isMC_;
  bool OnlyGen_;
  bool mumuMassConstraint_;
  std::vector<double> mumuMasscut_;
  double Trkmass_;
  std::vector<double> BarebMasscut_;
  std::vector<double> bMasscut_;
  bool debug_;

  TTree*      tree_;
  TTree*      tree_muons;
  TTree*      tree_gen_muons;
  TTree*      tree_L1muons;
  TTree*      tree_L2muons;
  TTree*      tree_L3muons;
  
  int         mu1_charge, mu2_charge;
  int         mu1_GEN_match, mu2_GEN_match;
  int         mu1_L1_match, mu2_L1_match;
  int         mu1_L2_match, mu2_L2_match;
  int         mu1_L3_match, mu2_L3_match;

  Double_t    mu1C2;
  int         mu1NHits, mu1NPHits; 
  Double_t    mu2C2;
  int         mu2NHits, mu2NPHits;
  Double_t    mu1dxy, mu2dxy, mu1dz, mu2dz;
  Double_t    muon_dca;

  std::vector<float> L1mu_pt, L1mu_eta, L1mu_phi, L1mu_etaAtVtx, L1mu_phiAtVtx, L1mu_charge, L1mu_quality;
  std::vector<float> mu_pt, mu_eta, mu_phi, mu_charge;
  std::vector<float> L2mu_pt, L2mu_eta, L2mu_phi;
  std::vector<float> L3mu_pt, L3mu_eta, L3mu_phi;
  std::vector<int> hltsVector, l1sVector;
  std::vector<int> mu1_hltsVector, mu2_hltsVector;


  int         HLT_Dim25, HLT_JpsiTrk_Bc, HLT_JpsiTk; 
  int         HLT_DMu4_3_LM;//HLT_DoubleMu4_3_LowMass
  int         HLT_DMu4_LM_Displaced;//HLT_DoubleMu4_LowMass_Displaced

  bool       mu1soft, mu2soft, mu1tight, mu2tight;  
  bool       mu1PF, mu2PF, mu1loose, mu2loose;  
  bool       mu1Tracker, mu2Tracker, mu1Global, mu2Global;  
 
  // *************************************
  unsigned int    nB;
  unsigned int    nMu;
    
  Double_t DiMu_mass,DiMu_mass_err;
  Double_t DiMu_pt, DiMu_eta, DiMu_phi;
  Double_t DiMu_mu1_pt, DiMu_mu1_eta, DiMu_mu1_phi;
  Double_t DiMu_mu2_pt, DiMu_mu2_eta, DiMu_mu2_phi;
  //Double_t DiMu_mu1_charge, DiMu_mu2_charge;

  // Double_t       B_mass, B_px, B_py, B_pz, B_charge;
  // Double_t       B_k_px, B_k_py, B_k_pz,  B_k_charge1; 
  // Double_t       B_k_px_track, B_k_py_track, B_k_pz_track;

  // Double_t       B_J_mass, B_J_massErr, B_J_px, B_J_py, B_J_pz;

  // Double_t       B_J_pt1, B_J_px1, B_J_py1, B_J_pz1;
  // Double_t       B_J_pt2, B_J_px2, B_J_py2, B_J_pz2;
  // int            B_J_charge1, B_J_charge2;

  // Primary Vertex (PV)
  UInt_t          nVtx;
  Double_t        priVtxX, priVtxY, priVtxZ, priVtxXE, priVtxYE, priVtxZE, priVtxCL;
  Double_t        priVtxXYE, priVtxXZE, priVtxYZE;
  
  // ********************************** ************************************************************************
 
  Double_t      DiMu_chi2, DiMu_Prob;

  Double_t      DiMu_DecayVtxX,  DiMu_DecayVtxY,  DiMu_DecayVtxZ;
  Double_t      DiMu_DecayVtxXE, DiMu_DecayVtxYE, DiMu_DecayVtxZE;
  Double_t      DiMu_DecayVtxXYE, DiMu_DecayVtxXZE, DiMu_DecayVtxYZE;
  Double_t      lxy, lxyerr;
  Double_t      lxy_pv, lxy_pv_err;

  UInt_t  run;
  ULong64_t event;
  UInt_t lumiblock;

  TLorentzVector gen_bc_p4,gen_jpsi_p4,gen_pion3_p4,gen_muon1_p4,gen_muon2_p4;
  TVector3       gen_bc_vtx,gen_jpsi_vtx;
  Double_t       gen_bc_ct;
  std::vector<float> GENmu_pt, GENmu_eta, GENmu_phi;
  std::vector<int> GENmu_charge, GENmu_status;
  std::vector<int> GENmu_mother, GENmu_grandmother; 


};
#endif
