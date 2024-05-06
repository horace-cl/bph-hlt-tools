// -*- C++ -*-
//
// Package:    MuMu
// Class:      MuMu
// 

//=================================================
// Original author:  Jhovanny Andres Mejia        |
//         created:  October of 2021              |
//         <jhovanny.andres.mejia.guisao@cern.ch> | 
//=================================================

// user include files
#include "myAnalyzers/bph-hlt-tools/src/MuMu.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TH2F.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"


#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

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

MuMu::MuMu(const edm::ParameterSet& iConfig)
  :
  //ttrkToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
  hltPrescale_(iConfig, consumesCollector(), *this),
  ttrkToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
  muon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  genCands_(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"))), 
  packedGenToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter <edm::InputTag> ("packedGenParticles"))), 
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),
  triggerCollection_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("TriggerInput"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  algTok_(consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<edm::InputTag>("algInputTag"))),
  algInputTag_(consumes<GlobalAlgBlkBxCollection>(iConfig.getParameter<edm::InputTag>("algInputTag"))),
  l1MuonsToken_(consumes<BXVector<l1t::Muon>>(iConfig.getParameter<edm::InputTag>("l1Muons"))),
  HLTPaths_(iConfig.getParameter<std::vector<std::string>>("HLTPaths")),
  L1Seeds_(iConfig.getParameter<std::vector<std::string>>("L1Seeds")),    
  gtUtil_( new l1t::L1TGlobalUtil( iConfig, consumesCollector(), *this, iConfig.getParameter<edm::InputTag>("algInputTag"), iConfig.getParameter<edm::InputTag>("algInputTag"), l1t::UseEventSetupIn::RunAndEvent  )),
  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
  mumuMassConstraint_(iConfig.getParameter<bool>("mumuMassConstraint")),
  mumuMasscut_(iConfig.getParameter<std::vector<double> >("mumuMasscut")),
  Trkmass_(iConfig.getParameter<double>("Trkmass")),
  BarebMasscut_(iConfig.getParameter<std::vector<double> >("BarebMasscut")),
  bMasscut_(iConfig.getParameter<std::vector<double> >("bMasscut")),
  debug_(iConfig.getParameter<bool>("debug")),
  tree_(0), tree_muons(0), tree_gen_muons(0), tree_L1muons(0), tree_L2muons(0), tree_L3muons(0),

  mu1_charge(0), mu2_charge(0),
  mu1_L1_match(0), mu2_L1_match(0),
  mu1_L2_match(0), mu2_L2_match(0),
  mu1_L3_match(0), mu2_L3_match(0),
  mu1C2(0), mu1NHits(0), mu1NPHits(0),
  mu2C2(0), mu2NHits(0), mu2NPHits(0),
  mu1dxy(0), mu2dxy(0), mu1dz(0), mu2dz(0),
  muon_dca(0),

  
  L1mu_pt(0), L1mu_eta(0), L1mu_phi(0), L1mu_etaAtVtx(0), L1mu_phiAtVtx(0), L1mu_charge(0), L1mu_quality(0),
  mu_pt(0), mu_eta(0), mu_phi(0), mu_charge(0),
  L2mu_pt(0), L2mu_eta(0), L2mu_phi(0),
  L3mu_pt(0), L3mu_eta(0), L3mu_phi(0),

  hltsVector(HLTPaths_.size()),
  l1sVector(L1Seeds_.size()),
  mu1_hltsVector(HLTPaths_.size()),
  mu2_hltsVector(HLTPaths_.size()),

  HLT_Dim25(0), HLT_JpsiTrk_Bc(0), HLT_JpsiTk(0),
  HLT_DMu4_3_LM(0), HLT_DMu4_LM_Displaced(0),
  
  mu1soft(0), mu2soft(0), mu1tight(0), mu2tight(0), 
  mu1PF(0), mu2PF(0), mu1loose(0), mu2loose(0),
  mu1Tracker(0), mu2Tracker(0), mu1Global(0), mu2Global(0),  
 
  // *******************************************************
 
  nB(0), nMu(0),
  //B_mass(0), B_px(0), B_py(0), B_pz(0), B_charge(0),
  //B_k_px(0), B_k_py(0), B_k_pz(0), B_k_charge1(0),
  //B_k_px_track(0), B_k_py_track(0), B_k_pz_track(0),
  
  DiMu_mass(0), DiMu_mass_err(0), 
  DiMu_pt(0), DiMu_eta(0), DiMu_phi(0),
  DiMu_mu1_pt(0), DiMu_mu1_eta(0), DiMu_mu1_phi(0), 
  DiMu_mu2_pt(0), DiMu_mu2_eta(0), DiMu_mu2_phi(0), 

  // Primary Vertex (PV)
  nVtx(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priVtxXYE(0), priVtxXZE(0), priVtxYZE(0),
  
  // ************************ ****************************************************

  DiMu_chi2(0), DiMu_Prob(0), 
 
  DiMu_DecayVtxX(0),     DiMu_DecayVtxY(0),     DiMu_DecayVtxZ(0),
  DiMu_DecayVtxXE(0),    DiMu_DecayVtxYE(0),    DiMu_DecayVtxZE(0),
  DiMu_DecayVtxXYE(0),   DiMu_DecayVtxXZE(0),   DiMu_DecayVtxYZE(0),

  lxy(0), lxyerr(0), lxy_pv(0), lxy_pv_err(0),

  run(0), event(0),
  lumiblock(0),
  GENmu_pt(0), GENmu_eta(0), GENmu_phi(0),
  GENmu_charge(0), GENmu_status(0), GENmu_mother(0), GENmu_grandmother(0)
{
   //now do what ever initialization is needed
}

//MuMu::~MuMu(){}
void MuMu::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
//--------------------------------------------------------------------------
{
    //using namespace std;
    //using namespace edm;

    bool changed(true);
    //hltConfig_.init(iRun,iSetup,"HLT",changed);
    hltPrescale_.init(iRun, iSetup, "HLT", changed);
    //hltConfig_ = hltPrescale_.hltConfigProvider();
    
    if (changed)
    {
        std::cout<<"HLTConfig has changed . . . "<<std::endl;
        
    }
    
}//------------------- beginRun()

void MuMu::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;  
  desc.add<unsigned int>("stageL1Trigger", 1);
  descriptions.add("isoTrigDefault", desc);
}

// ------------ method called to for each event  ------------
void MuMu::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //bool debug_ = true;
  if (debug_) std::cout << " - Welcome to the Analyzer " << std::endl;
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  //*********************************
  // Get event content information
  //*********************************  
 
  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB = iSetup.getHandle(ttrkToken_);
  //edm::ESHandle<TransientTrackBuilder> theB; 
  //iSetup.getHandle<TransientTrackRecord>().get(ttrkToken_,theB);//
  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);// old one 
  //auto const& theB = iSetup.getHandle(transientTrackRecordToken_);
  //const edm::ESHandle<TransientTrackBuilder> TTbuilder = setup.getHandle(ttrkToken_);

  edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);

  edm::Handle< View<pat::Muon> > thePATMuonHandle;
  iEvent.getByToken(muon_Label,thePATMuonHandle);

  edm::Handle<reco::GenParticleCollection> pruned;
  //edm::Handle<pat::PackedGenParticle> pruned; 
  iEvent.getByToken(genCands_, pruned);
  
  edm::Handle<pat::PackedGenParticleCollection> packed;
  iEvent.getByToken(packedGenToken_,packed);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerCollection;
  iEvent.getByToken(triggerCollection_, triggerCollection);

  reco::BeamSpot vertexBeamSpot;
  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByToken(BSLabel_, beamSpot);
  vertexBeamSpot = *beamSpot;


  lumiblock = iEvent.id().luminosityBlock();
  run       = iEvent.id().run();
  event     = iEvent.id().event();



  //*********************************
  // Get gen level information
  //*********************************
  
  gen_bc_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_pion3_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_bc_vtx.SetXYZ(0.,0.,0.);
  gen_jpsi_vtx.SetXYZ(0.,0.,0.);
  gen_bc_ct = -9999.;

  // if(debug_) std::cout << " pruned.isValid() = "  << pruned.isValid() << std::endl;
  // if(debug_) std::cout << " isMC_            = "  << isMC_    << std::endl;
  // if(debug_) std::cout << " OnlyGen_         = "  << OnlyGen_ << std::endl;

  if ( (isMC_ || OnlyGen_) && pruned.isValid() ) {
    int foundit = 0;
	  if(debug_) std::cout<< "---> isMC" << std::endl;
    if(debug_) std::cout<< " Pruned size: " << pruned->size() << std::endl;
    for (size_t i=0; i<pruned->size(); i++) {
      foundit = 0;
      const reco::Candidate *dau = &(*pruned)[i];
      
      unsigned int nMom  = dau->numberOfMothers(); 
      // if (debug_) std::cout<< "      Particle   ID: " << dau->pdgId() << std::endl;
      // if (debug_) std::cout<< "             Status: " << dau->status() << std::endl;
      // if (debug_) std::cout<< "          n MoThers: " << nMom << std::endl;
      // if (nMom>0) {
      //   const reco::Candidate *mom_ = dau->mother(0);  
      //   if (debug_) std::cout<< "      1st Mother ID: " << mom_->pdgId() << std::endl;
      // }
      if ( (abs(dau->pdgId()) == 521) ) { //&& (dau->status() == 2) ) {
	foundit++;
	gen_bc_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
	gen_bc_vtx.SetXYZ(dau->vx(),dau->vy(),dau->vz());
	for (size_t k=0; k<dau->numberOfDaughters(); k++) {
	  const reco::Candidate *gdau = dau->daughter(k);
	  if (gdau->pdgId()==443 ) { //&& gdau->status()==2) {
	    foundit++;
	    gen_jpsi_vtx.SetXYZ(gdau->vx(),gdau->vy(),gdau->vz());
	    gen_bc_ct = GetLifetime(gen_bc_p4,gen_bc_vtx,gen_jpsi_vtx);
	    int nm=0;
	    for (size_t l=0; l<gdau->numberOfDaughters(); l++) {
	      const reco::Candidate *mm = gdau->daughter(l);
	      if (mm->pdgId()==13) { foundit++;
		if (mm->status()!=1) {
		  for (size_t m=0; m<mm->numberOfDaughters(); m++) {
		    const reco::Candidate *mu = mm->daughter(m);
		    if (mu->pdgId()==13 ) { //&& mu->status()==1) {
		      nm++;
		      gen_muon1_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
		      break;
		    }
		  }
		} else {
		  gen_muon1_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
		  nm++;
		}
	      }
	      if (mm->pdgId()==-13) { foundit++;
		if (mm->status()!=1) {
		  for (size_t m=0; m<mm->numberOfDaughters(); m++) {
		    const reco::Candidate *mu = mm->daughter(m);
		    if (mu->pdgId()==-13 ) { //&& mu->status()==1) {
		      nm++;
		      gen_muon2_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
		      break;
		    }
		  }
		} else {
		  gen_muon2_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
		  nm++;
		}
	      }
	    }
	    if (nm==2) gen_jpsi_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
	    else foundit-=nm;
	  }
	} // for (size_t k

	for (size_t lk=0; lk<packed->size(); lk++) {
	  const reco::Candidate * dauInPrunedColl = (*packed)[lk].mother(0);
	  int stable_id = (*packed)[lk].pdgId();
	  
	  if (dauInPrunedColl != nullptr && isAncestor(dau,dauInPrunedColl)) {
	    if( abs(stable_id) == 321 || abs(stable_id) ==211 ) {foundit++;
	      gen_pion3_p4.SetPtEtaPhiM((*packed)[lk].pt(),(*packed)[lk].eta(),(*packed)[lk].phi(),(*packed)[lk].mass());
	    }
	  }
	}
	//if ((abs(gdau->pdgId())==211 || abs(gdau->pdgId())==321) && gdau->status()==1) { foundit++;
	//gen_pion3_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
	//}
	//} // for (size_t k
      }   // if (abs(dau->pdgId())==521 )
      if (foundit>=5) break;
    } // for i
    if (foundit!=5) {
      gen_bc_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
      gen_bc_vtx.SetXYZ(0.,0.,0.);
      gen_jpsi_vtx.SetXYZ(0.,0.,0.);
      gen_bc_ct = -9999.;
      //std::cout << "Does not found the given decay " << run << "," << event << " foundit=" << foundit << std::endl; // sanity check
    }
  }
 
 
 
 
 
 
 
 
 
  //*********************************
  //Now we get the primary vertex 
  //*********************************
  // if(debug_) std::cout<< "---> Primary Vertex" << std::endl;

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
 


  // Here we are going to save the information of the HLT_paths and L1seeds
  // if(debug_) std::cout<< "---> HLT and L1s" << std::endl;
  hltsVector.resize( HLTPaths_.size(), 0);
  mu1_hltsVector.resize( HLTPaths_.size(), 0);
  mu2_hltsVector.resize( HLTPaths_.size(), 0);
  l1sVector.resize( L1Seeds_.size(), 0);

  //Unpack trigger info
  edm::Handle<edm::TriggerResults> triggerBits;
  //hltConfig_.init(iRun,iSetup,processName_,changed);

  iEvent.getByToken(triggerResults_Label, triggerBits);
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  const pat::TriggerObjectStandAloneCollection unPackedCollection;

  // if(debug_) std::cout<< "   ----> HLT only" << std::endl;
    for (unsigned int i = 0; i  < triggerBits->size(); i++) {
      std::string iName = names.triggerName(i);
      if (triggerBits->accept(i)) {
        for(std::size_t i = 0; i < HLTPaths_.size(); ++i) {
          if (iName.find(HLTPaths_[i]) != std::string::npos) hltsVector[i] = 1;
        }

      }
  }



  // if(debug_) std::cout<< "   ----> L1 only" << std::endl;
  // Initialize the L1TGlobalUtil object
  gtUtil_->retrieveL1(iEvent, iSetup, algInputTag_);
  
  // Some attributes of the L1TGlobalUtil object:
  // // L1TGlobalUtil:  Utility class for parsing the L1 Trigger Menu
  // // L1Trigger/L1TGlobal/L1TGlobal/interface/L1TGlobalUtil
  const vector<pair<string, bool> > decisionsFinal = gtUtil_->decisionsFinal();
  const vector<pair<string, bool> > decisionsInitial = gtUtil_->decisionsInitial();
  const vector<pair<string, bool> > decisionsInterm = gtUtil_->decisionsInterm();
  const vector<pair<string, double> > prescales = gtUtil_->prescales();
  const vector<pair<string, vector<int> > > masks = gtUtil_->masks();

  // if (debug_) std::cout << "decisionsInitial (size) = " << decisionsInitial.size() << std::endl;
  // if (debug_) std::cout << "decisionsInterm  (size) = " << decisionsInterm.size() << std::endl;
  // if (debug_) std::cout << "decisionsFinal   (size) = " << decisionsFinal.size() << std::endl;
  // if (debug_) std::cout << "prescales        (size) = " << prescales.size() << std::endl;
  // if (debug_) std::cout << "masks            (size) = " << masks.size() << std::endl;

  const string gtTriggerMenuName    = gtUtil_->gtTriggerMenuName();
  const string gtTriggerMenuVersion = gtUtil_->gtTriggerMenuVersion();
  const string gtTriggerMenuComment = gtUtil_->gtTriggerMenuComment();

  // if (debug_) std::cout << "gtTriggerMenuName    = " << gtTriggerMenuName << std::endl;
  // if (debug_) std::cout << "gtTriggerMenuVersion = " << gtTriggerMenuVersion << std::endl;
  // if (debug_) std::cout << "gtTriggerMenuComment = " << gtTriggerMenuComment << std::endl;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // EXAMPLE: L1 and HLT prescale values via (L1) EventSetup
  // Current (default) prescale set index - to be taken from L1GtUtil via Event.
  // Try to get the L1 and HLT prescale values that were actually used 
  // for this event.
  // This example needs the conditions stored in the Global Tag,
  // which is some sort of snapshot of the by-then current detector
  // conditions.  They need to be extracted from the database
  // and for that the "Frontier Conditions" lines need to
  // be added in the python configuration file along with the 
  // name for the global tag.
  // This can make the job very slow at the very begining....
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //const std::pair<int,int> prescalesHLT(hltConfig_.prescaleValue<double>(iSetup,HLTPaths_[0]));
  //const std::pair<int,int> prescales(hltConfig_.prescaleValues(iEvent,iSetup,HLTPaths_[0]));
  auto const prescalesHLT = hltPrescale_.prescaleValuesInDetail<double>(iEvent, iSetup, HLTPaths_[0]);
  std::cout << "analyzeSimplePrescales: path "
      << HLTPaths_[0] << " [ - ] "
      << "prescales L1T,HLT: " 
      << prescalesHLT.first.at(0).first  << " - " 
      << prescalesHLT.first.at(0).second << "," 
      << prescalesHLT.second
      << std::endl;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for (size_t i_l1t = 0; i_l1t < decisionsFinal.size(); i_l1t++){
    string l1tName = (decisionsFinal.at(i_l1t)).first;
    if (debug_ && false){
      for(std::size_t i_input = 0; i_input < L1Seeds_.size(); ++i_input) {
        if (l1tName.find(L1Seeds_[i_input]) != std::string::npos){
          std::cout << "\t - Name = " << l1tName << std::endl;
          std::cout << "\t - - decisionsInitial = " << decisionsInitial.at(i_l1t).second << std::endl;
          std::cout << "\t - - decisionsInterm  = " << decisionsInterm.at(i_l1t).second << std::endl;
          std::cout << "\t - - decisionsFinal   = " << decisionsFinal.at(i_l1t).second << std::endl;
          std::cout << "\t - - prescales        = " << prescales.at(i_l1t).second << std::endl;
          std::cout << "\t - - masks  (size)    = " << masks.at(i_l1t).second.size() << std::endl;
        }
      }
    }

    if ((decisionsFinal.at(i_l1t)).second == 1) {
      for(std::size_t i = 0; i < L1Seeds_.size(); ++i) {
        if (l1tName.find(L1Seeds_[i]) != std::string::npos) l1sVector[i] = 1;
      }

    }
  }




  edm::Handle<BXVector<l1t::Muon> > gmuons;
  iEvent.getByToken(l1MuonsToken_, gmuons);


  if (debug_) std::cout << "Len Triggers : " << triggerCollection->size() << std::endl;
  for (pat::TriggerObjectStandAlone trig : *triggerCollection) {
      trig.unpackPathNames(names);
      trig.unpackFilterLabels(iEvent, *triggerBits);
  }




  //*****************************************
  //Let's begin by looking for J/psi+K^+

  unsigned int nMu_tmp = thePATMuonHandle->size();
  nMu = nMu_tmp;

  for(View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1)  {     
    for(View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2)  {
      if(iMuon1==iMuon2) continue;

      //opposite charge 
      //if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue;

      TrackRef glbTrack1;	  
      TrackRef glbTrack2;	  
      
      glbTrack1 = iMuon1->track();
      glbTrack2 = iMuon2->track();

      if( glbTrack1.isNull() || glbTrack2.isNull() ) 
         {
         //std::cout << "continue due to no track ref" << endl;
         continue;
         }

      if(iMuon1->track()->pt()<1.5) continue;
      if(iMuon2->track()->pt()<1.5) continue;

      if(!(glbTrack2->quality(reco::TrackBase::highPurity))) continue;
      if(!(glbTrack1->quality(reco::TrackBase::highPurity))) continue;	 

      reco::TransientTrack muon1TT((*theB).build(glbTrack1));
      reco::TransientTrack muon2TT((*theB).build(glbTrack2));

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
         std::cout<<" Exception caught ... continuing 1 "<<std::endl; 
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

      if (!psiVertexFitTree->isValid())  continue; 

      psiVertexFitTree->movePointerToTheTop();

      RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();
      RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();

      if( psi_vFit_vertex_noMC->chiSquared() < 0 )  continue;

      //some loose cuts go here	  
      if(psi_vFit_vertex_noMC->chiSquared()>50.) continue;
      // ******************************************** JPsi mass  ********************************************
      //if(psi_vFit_noMC->currentState().mass()<2.9 || psi_vFit_noMC->currentState().mass()>3.3) continue;
      // ******************************************** mu1u mass  ********************************************
      //if(psi_vFit_noMC->currentState().mass()<0.2 || psi_vFit_noMC->currentState().mass()>4.9) continue;

      if(psi_vFit_noMC->currentState().mass()<mumuMasscut_[0] || psi_vFit_noMC->currentState().mass()>mumuMasscut_[1]) continue;


      double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());

      if(J_Prob_tmp<0.001) continue;

      psiVertexFitTree->movePointerToTheFirstChild();
      RefCountedKinematicParticle mu1Cand = psiVertexFitTree->currentParticle();

      psiVertexFitTree->movePointerToTheNextChild();
      RefCountedKinematicParticle mu2Cand = psiVertexFitTree->currentParticle();

      // KinematicParameters psiMu1KP = mu1Cand->currentState().kinematicParameters();
      // KinematicParameters psiMu2KP = mu2Cand->currentState().kinematicParameters();
      KinematicState mu1cand_state = mu1Cand->currentState();
      KinematicState mu2cand_state = mu2Cand->currentState();

         
      DiMu_mass     =  psi_vFit_noMC->currentState().mass() ;
      DiMu_mass_err = sqrt(psi_vFit_noMC->currentState().kinematicParametersError().matrix()(6,6));
         
      DiMu_pt  = psi_vFit_noMC->currentState().globalMomentum().perp() ;
      DiMu_eta = psi_vFit_noMC->currentState().globalMomentum().eta() ;
      DiMu_phi = psi_vFit_noMC->currentState().globalMomentum().phi() ;

      DiMu_mu1_pt  = mu1cand_state.globalMomentum().perp();
      DiMu_mu1_eta = mu1cand_state.globalMomentum().eta();
      DiMu_mu1_phi = mu1cand_state.globalMomentum().phi();

      DiMu_mu2_pt  = mu2cand_state.globalMomentum().perp();
      DiMu_mu2_eta = mu2cand_state.globalMomentum().eta();
      DiMu_mu2_phi = mu2cand_state.globalMomentum().phi();

      /*B_J_pt2 = Jp2vec.perp();
      B_J_px2 = psiMu2KP.momentum().x();
      B_J_py2 = psiMu2KP.momentum().y();
      B_J_pz2 = psiMu2KP.momentum().z();
      B_J_charge2 = mu2Cand->currentState().particleCharge();
      */
      DiMu_chi2 = psi_vFit_vertex_noMC->chiSquared();
      //double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
      //double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
      DiMu_Prob   = J_Prob_tmp;

      DiMu_DecayVtxX  = (*psi_vFit_vertex_noMC).position().x();    
      DiMu_DecayVtxY  = (*psi_vFit_vertex_noMC).position().y();
      DiMu_DecayVtxZ  = (*psi_vFit_vertex_noMC).position().z();
      DiMu_DecayVtxXE  = psi_vFit_vertex_noMC->error().cxx();   
      DiMu_DecayVtxYE  = psi_vFit_vertex_noMC->error().cyy();   
      DiMu_DecayVtxZE  = psi_vFit_vertex_noMC->error().czz();
      DiMu_DecayVtxXYE  = psi_vFit_vertex_noMC->error().cyx();
      DiMu_DecayVtxXZE  = psi_vFit_vertex_noMC->error().czx();
      DiMu_DecayVtxYZE  = psi_vFit_vertex_noMC->error().czy();
      

      /// lxy HLT Way
      GlobalPoint secondaryVertex(DiMu_DecayVtxX, DiMu_DecayVtxY, DiMu_DecayVtxZ);
      GlobalError verr = psi_vFit_vertex_noMC->error();
      GlobalPoint displacementFromBeamspot(-1 * ((vertexBeamSpot.x0() - secondaryVertex.x()) +
                                               (secondaryVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dxdz()),
                                         -1 * ((vertexBeamSpot.y0() - secondaryVertex.y()) +
                                               (secondaryVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dydz()),
                                         0);
      lxy = displacementFromBeamspot.perp();
      lxyerr = sqrt(verr.rerr(displacementFromBeamspot));


      /// /// lxy aa way
      float dx_lxy = (DiMu_DecayVtxX-bestVtx.x());
      float dy_lxy = (DiMu_DecayVtxY-bestVtx.y());      
      //float dz_lxy = (DiMu_DecayVtxZ-bestVtx.z());
      math::XYZVector vPS( dx_lxy, dy_lxy, 0);
      float lxy_xhat  = (dx_lxy*dx_lxy)/vPS.perp2() ;
      float lxy_xyhat  = (dx_lxy*dy_lxy)/vPS.perp2();
      float lxy_yhat = (dy_lxy*dy_lxy)/vPS.perp2();

      float sp_xx = psi_vFit_vertex_noMC->error().cxx()+bestVtx.covariance(0,0);
      float sp_xy = psi_vFit_vertex_noMC->error().cyx()+bestVtx.covariance(0,1);
      float sp_yy = psi_vFit_vertex_noMC->error().cyy()+bestVtx.covariance(1,1);

      lxy_pv  = sqrt(vPS.perp2());
      lxy_pv_err = sqrt(lxy_xhat*sp_xx +2*lxy_xyhat*sp_xy + lxy_yhat*sp_yy);






      if (debug_) std::cout << " Trigger Matching " << std::endl;
      // ********************* muon-trigger-machint ****************

      const pat::Muon* muon1 = &(*iMuon1);
      const pat::Muon* muon2 = &(*iMuon2);

      for(std::size_t i = 0; i < HLTPaths_.size(); ++i) {
        if (muon1->triggerObjectMatchByPath((HLTPaths_[i]+"*").c_str())!=nullptr) mu1_hltsVector[i] = 1;
        if (muon2->triggerObjectMatchByPath((HLTPaths_[i]+"*").c_str())!=nullptr) mu2_hltsVector[i] = 1;
      }

      int HLT_Dim25_tmp = 0, HLT_JpsiTk_tmp = 0,  HLT_JpsiTrk_Bc_tmp = 0, HLT_DMu4_3_LM_tmp = 0, HLT_DMu4_LM_Displaced_tmp = 0; 

      if(muon1->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr) HLT_Dim25_tmp = 1;
      if(muon1->triggerObjectMatchByPath("HLT_DoubleMu4_MuMuTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_MuMuTrk_Displaced_v*")!=nullptr) HLT_JpsiTk_tmp = 1;
      if(muon1->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Bc_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Bc_v*")!=nullptr) HLT_JpsiTrk_Bc_tmp = 1;
      if(muon1->triggerObjectMatchByPath("HLT_DoubleMu4_3_LowMass_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_3_LowMass_v*")!=nullptr) HLT_DMu4_3_LM_tmp = 1;
      if(muon1->triggerObjectMatchByPath("HLT_DoubleMu4_LowMass_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_LowMass_Displaced_v*")!=nullptr) HLT_DMu4_LM_Displaced_tmp = 1;


      HLT_Dim25 =  HLT_Dim25_tmp ;	       
      HLT_JpsiTk =  HLT_JpsiTk_tmp ;
      HLT_JpsiTrk_Bc =  HLT_JpsiTrk_Bc_tmp ;
      HLT_DMu4_3_LM =  HLT_DMu4_3_LM_tmp ;
      HLT_DMu4_LM_Displaced =  HLT_DMu4_LM_Displaced_tmp ;

      if (debug_) std::cout << " L1  Matching " << std::endl;
      // ************ l1, l2, l3 ************ 
      float dR2_muon1, dR2_muon2; 
      double dR2_threshold;
      dR2_threshold = 0.5 * 0.5;
      mu1_L1_match = 0;
      mu2_L1_match = 0;
      for (auto itr = gmuons->begin(0); itr != gmuons->end(0); ++itr) {
        dR2_muon1 = reco::deltaR2(iMuon1->eta(),iMuon1->phi(), 
                                  itr->eta(), itr->phi());                                
        if (dR2_muon1 < dR2_threshold) mu1_L1_match=1; 

        dR2_muon2 = reco::deltaR2(iMuon2->eta(),iMuon2->phi(), 
                                  itr->eta(), itr->phi());
        if (dR2_muon2 < dR2_threshold) mu2_L1_match=1;

        if (mu1_L1_match==1 && mu2_L1_match==1) break;
      }
      

      if (debug_) std::cout << " L2, L3  Matching " << std::endl;
      mu1_L2_match = 0;
      mu2_L2_match = 0;
      mu1_L3_match = 0;
      mu2_L3_match = 0;

      dR2_threshold = 0.1 * 0.1;
      std::string hltMuColl_L2 = "hltL2MuonCandidates";
      std::string hltMuColl_L3 = "hltIterL3MuonCandidates";
      for (pat::TriggerObjectStandAlone obj : *triggerCollection) {

        if( obj.hasCollection(hltMuColl_L2) ) {
          dR2_muon1 = reco::deltaR2(iMuon1->eta(),iMuon1->phi(), 
                                    obj.eta(), obj.phi());                                
          if (dR2_muon1 < dR2_threshold) mu1_L2_match=1; 

          dR2_muon2 = reco::deltaR2(iMuon2->eta(),iMuon2->phi(), 
                                  obj.eta(), obj.phi());
          if (dR2_muon1 < dR2_threshold) mu2_L2_match=1; 
        }


        if( obj.hasCollection(hltMuColl_L3) ) {
          dR2_muon1 = reco::deltaR2(iMuon1->eta(),iMuon1->phi(), 
                                    obj.eta(), obj.phi());                                
          if (dR2_muon1 < dR2_threshold) mu1_L3_match=1; 

          dR2_muon2 = reco::deltaR2(iMuon2->eta(),iMuon2->phi(), 
                                  obj.eta(), obj.phi());
          if (dR2_muon1 < dR2_threshold) mu2_L3_match=1; 
        }
        
      }


      if (debug_) std::cout << " GEN  Matching " << std::endl;
      mu1_GEN_match = 0;
      mu2_GEN_match = 0;  
      if ( (isMC_ || OnlyGen_) && pruned.isValid() ){
        for (size_t i=0; i<pruned->size(); i++) {
          const reco::Candidate *gen_particle = &(*pruned)[i];    
          if (abs(gen_particle->pdgId())==13){
            dR2_muon1 = reco::deltaR2(iMuon1->eta(),iMuon1->phi(), 
                                      gen_particle->eta(), gen_particle->phi());                                
            if (dR2_muon1 < dR2_threshold) mu1_GEN_match=1; 

            dR2_muon2 = reco::deltaR2(iMuon2->eta(),iMuon2->phi(), 
                                      gen_particle->eta(), gen_particle->phi());
            if (dR2_muon2 < dR2_threshold) mu2_GEN_match=1;

            if (mu1_GEN_match==1 && mu2_GEN_match==1) break;
          
          }
        }
        

      }

      if (debug_) std::cout << " Muons IDs and properties " << std::endl;
      // ************ Different muons Id, and other properties  ****************
      mu1_charge = iMuon1->charge() ;
      mu2_charge = iMuon2->charge() ; 
      mu1soft    = iMuon1->isSoftMuon(bestVtx) ;
      mu2soft    = iMuon2->isSoftMuon(bestVtx) ;
      mu1tight   = iMuon1->isTightMuon(bestVtx) ;
      mu2tight   = iMuon2->isTightMuon(bestVtx) ;
      mu1PF = iMuon1->isPFMuon();
      mu2PF = iMuon2->isPFMuon();
      mu1Tracker = iMuon1->isTrackerMuon();
      mu2Tracker = iMuon2->isTrackerMuon();
      mu1Global = iMuon1->isGlobalMuon();
      mu2Global = iMuon2->isGlobalMuon();
      mu1loose = muon::isLooseMuon(*iMuon1);
      mu2loose = muon::isLooseMuon(*iMuon2);

      mu1C2 =  glbTrack1->normalizedChi2() ;
      mu1NHits =  glbTrack1->numberOfValidHits() ;
      mu1NPHits =  glbTrack1->hitPattern().numberOfValidPixelHits() ;	       
      mu2C2 =  glbTrack2->normalizedChi2() ;
      mu2NHits =  glbTrack2->numberOfValidHits() ;
      mu2NPHits =  glbTrack2->hitPattern().numberOfValidPixelHits() ;
      mu1dxy = glbTrack1->dxy(bestVtx.position()) ;// 
      mu2dxy = glbTrack2->dxy(bestVtx.position()) ;// 
      mu1dz = glbTrack1->dz(bestVtx.position()) ;
      mu2dz = glbTrack2->dz(bestVtx.position()) ;
      muon_dca = dca;

      //fill the tree
      tree_->Fill();

      nB++;	       
      muonParticles.clear();
      //vFitMCParticles.clear();

	    
	  }
  }
 

 
  if (nB>0) {

    if (debug_) std::cout << " ---> L1 Tree " << std::endl;
    for (auto itr = gmuons->begin(0); itr != gmuons->end(0); ++itr) {
      L1mu_pt.push_back(itr->pt());
      L1mu_eta.push_back(itr->eta());
      L1mu_phi.push_back(itr->phi());
      L1mu_etaAtVtx.push_back(itr->etaAtVtx());
      L1mu_phiAtVtx.push_back(itr->phiAtVtx());
      L1mu_quality.push_back(itr->hwQual());
      L1mu_charge.push_back(itr->charge());
    }
    tree_L1muons->Fill();


    if (debug_) std::cout << " ---> L3 Tree " << std::endl;
    std::string hltMuColl = "hltIterL3MuonCandidates";
    for (pat::TriggerObjectStandAlone obj : *triggerCollection) {
      if( ! obj.hasCollection(hltMuColl) ) continue;
      L3mu_eta.push_back(obj.eta());
      L3mu_phi.push_back(obj.phi());
      L3mu_pt.push_back(obj.pt());
    }
    tree_L3muons->Fill();


    //hltMuColl = "hltIterL2MuonCandidates"; hltIterL3OIL3Muons
    if (debug_) std::cout << " ---> L2 Tree " << std::endl;
    hltMuColl = "hltL2MuonCandidates";
    for (pat::TriggerObjectStandAlone obj : *triggerCollection) {
      if( ! obj.hasCollection(hltMuColl) ) continue;
      L2mu_eta.push_back(obj.eta());
      L2mu_phi.push_back(obj.phi());
      L2mu_pt.push_back(obj.pt());
    }
    tree_L2muons->Fill();


    for(View<pat::Muon>::const_iterator iMuon = thePATMuonHandle->begin(); iMuon != thePATMuonHandle->end(); ++iMuon)  {    
      mu_pt.push_back(iMuon->pt());
      mu_eta.push_back(iMuon->eta());
      mu_phi.push_back(iMuon->phi());
      mu_charge.push_back(iMuon->charge());
    }
    tree_muons->Fill();

    if ( (isMC_ || OnlyGen_) && pruned.isValid() ) {
      if (debug_) std::cout << " ---> GEN Tree " << std::endl;
      for (size_t i=0; i<pruned->size(); i++) {
        const reco::Candidate *gen_particle = &(*pruned)[i];    
        if (abs(gen_particle->pdgId())==13){
          GENmu_pt.push_back(gen_particle->pt());
          GENmu_eta.push_back(gen_particle->eta());
          GENmu_phi.push_back(gen_particle->phi());
          GENmu_charge.push_back(gen_particle->charge());
          GENmu_status.push_back(gen_particle->status());
          GENmu_mother.push_back(gen_particle->mother(0)->pdgId());
          int n_mothers = gen_particle->mother(0)->numberOfMothers();
          if (n_mothers>0){
            GENmu_grandmother.push_back(gen_particle->mother(0)->mother(0)->pdgId());
          } 
          else{
            GENmu_grandmother.push_back(0);
          }
          

          // std::cout<< "gen_particle->pt     :" << gen_particle->pt() << std::endl;
          // std::cout<< "gen_particle->eta    :" << gen_particle->eta() << std::endl;
          // std::cout<< "gen_particle->phi    :" << gen_particle->phi() << std::endl;
          // std::cout<< "gen_particle->charge :" << gen_particle->charge() << std::endl;
          // std::cout<< "gen_particle->status :" << gen_particle->status() << std::endl;
          // std::cout<< "gen_particle->mother :" << gen_particle->mother(0)->pdgId() << std::endl;
          // std::cout<< "gen_particle->grandmother :" << gen_particle->mother(0)->mother(0)->pdgId() << std::endl;
          //std::cout<< "gen_particle->finalD :" << gen_particle->findDecayedMother() << std::endl;
        }
      }
      tree_gen_muons->Fill();
    }
    

  if (debug_) std::cout << " ---> Trees Filled\n\n " << std::endl;


  }




  
   // Initialize the variables

   nB = 0; nMu = 0;

	 DiMu_mass = 0; DiMu_mass_err = 0;
   DiMu_pt = 0;  DiMu_eta = 0;  DiMu_phi = 0;
   DiMu_mu1_pt = 0;  DiMu_mu1_eta = 0;  DiMu_mu1_phi = 0;
   DiMu_mu2_pt = 0;  DiMu_mu2_eta = 0;  DiMu_mu2_phi = 0;

   DiMu_chi2 = 0; 
   DiMu_Prob = 0;

   DiMu_DecayVtxX = 0;     DiMu_DecayVtxY = 0;     DiMu_DecayVtxZ = 0;
   DiMu_DecayVtxXE = 0;    DiMu_DecayVtxYE = 0;    DiMu_DecayVtxZE = 0;
   DiMu_DecayVtxXYE = 0;   DiMu_DecayVtxXZE = 0;   DiMu_DecayVtxYZE = 0;
   lxy = 0; lxyerr = 0;
   lxy_pv = 0; lxy_pv_err = 0;

   nVtx = 0;
   priVtxX = 0;     priVtxY = 0;     priVtxZ = 0; 
   priVtxXE = 0;    priVtxYE = 0;    priVtxZE = 0; priVtxCL = 0;
   priVtxXYE = 0;   priVtxXZE = 0;   priVtxYZE = 0;    

   mu1C2 = 0;
   mu1NHits = 0; mu1NPHits = 0;
   mu2C2 = 0;
   mu2NHits = 0; mu2NPHits = 0;
   mu1dxy = 0; mu2dxy = 0; mu1dz = 0; mu2dz = 0; muon_dca = 0;

   HLT_Dim25 = 0; HLT_JpsiTrk_Bc = 0; HLT_JpsiTk = 0;
   HLT_DMu4_3_LM = 0;  HLT_DMu4_LM_Displaced = 0;
 
   mu1soft = 0; mu2soft = 0; mu1tight = 0; mu2tight = 0;
   mu1PF = 0; mu2PF = 0; mu1loose = 0; mu2loose = 0; 
   mu1Tracker = 0; mu2Tracker = 0; mu1Global = 0; mu2Global = 0;  


  mu1_L1_match = 0;
  mu2_L1_match = 0;
  mu1_L2_match = 0;
  mu2_L2_match = 0;
  mu1_L3_match = 0;
  mu2_L3_match = 0;


   L1mu_pt.clear(); L1mu_eta.clear(); L1mu_phi.clear();
   L2mu_pt.clear(); L2mu_eta.clear(); L2mu_phi.clear();
   L3mu_pt.clear(); L3mu_eta.clear(); L3mu_phi.clear();
   GENmu_pt.clear(); GENmu_eta.clear(); GENmu_phi.clear();

   L1mu_etaAtVtx.clear(); L1mu_phiAtVtx.clear();
   L1mu_quality.clear(); L1mu_charge.clear();
  
   GENmu_charge.clear(); GENmu_status.clear(); 
   GENmu_mother.clear(); GENmu_grandmother.clear();

   mu_pt.clear(); mu_eta.clear(); mu_phi.clear();
   mu_charge.clear();

   for(std::size_t i = 0; i < HLTPaths_.size(); ++i) {
    hltsVector[i] = 0;
    mu1_hltsVector[i] = 0;
    mu2_hltsVector[i] = 0;
   }

   for(std::size_t i = 0; i < L1Seeds_.size(); ++i) l1sVector[i] = 0;

}

bool MuMu::IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

bool MuMu::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
    if (ancestor == particle ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(ancestor,particle->mother(i))) return true;
    }
    return false;
}

double MuMu::GetLifetime(TLorentzVector b_p4, TVector3 production_vtx, TVector3 decay_vtx) {
   TVector3 pv_dv = decay_vtx - production_vtx;
   TVector3 b_p3  = b_p4.Vect();
   pv_dv.SetZ(0.);
   b_p3.SetZ(0.);
   Double_t lxy   = pv_dv.Dot(b_p3)/b_p3.Mag();
   return lxy*b_p4.M()/b_p3.Mag();
}


// bool MuMu::checkDeltaR(double phi1, double eta1, const std::vector<float>& phiVec, const std::vector<float>& etaVec, double threshold) {
//     for (size_t i = 0; i < phiVec.size(); ++i) {
//         float dPhi = phiVec[i] - phi1;
//         float dEta = etaVec[i] - eta1;
//         float dR = reco::deltaR2(eta1,phi1, etaVec[i], phiVec[i]);
        
//         if (dR < threshold) {
//             return true;
//         }
//     }
//     return false;
// }

// int MuMu::MatchedObjs(const float eta,const float phi,const std::vector<pat::TriggerObjectStandAlone>& trigObjs,const float maxDeltaR=0.1, const ){
//     std::vector<const pat::TriggerObjectStandAlone*> matchedObjs;
//     const float maxDR2 = maxDeltaR*maxDeltaR;
//     for(auto& trigObj : trigObjs){
//       const float dR2 = reco::deltaR2(eta,phi,trigObj.eta(),trigObj.phi());
//       if(dR2<maxDR2) return 1;
//     }
//     return 0;
// }
// ------------ method called once each job just before starting event loop  ------------

void 
MuMu::beginJob()
{

  std::cout << "Beginning analyzer job with value of isMC= " << isMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_        = fs->make<TTree>("ntuple",         "LowMass Dimuons ntuple");
  tree_muons   = fs->make<TTree>("ntuple_muons",  "muons ntuple");
  tree_L1muons = fs->make<TTree>("ntuple_L1muons","L1muons ntuple");
  tree_L2muons = fs->make<TTree>("ntuple_L2muons","L2muons ntuple");
  tree_L3muons = fs->make<TTree>("ntuple_L3muons","L3muons ntuple");


  tree_->Branch("nB",&nB,"nB/i");
  tree_->Branch("nMu",&nMu,"nMu/i");

  //tree_->Branch("B_charge", &B_charge);
  tree_->Branch("DiMu_mass", &DiMu_mass);
  tree_->Branch("DiMu_pt" , &DiMu_pt);
  tree_->Branch("DiMu_eta", &DiMu_eta);
  tree_->Branch("DiMu_phi", &DiMu_phi);

  /*tree_->Branch("B_k_charge1", &B_k_charge1);
  tree_->Branch("B_k_px", &B_k_px);
  tree_->Branch("B_k_py", &B_k_py);
  tree_->Branch("B_k_pz", &B_k_pz);
  tree_->Branch("B_k_px_track", &B_k_px_track);
  tree_->Branch("B_k_py_track", &B_k_py_track);
  tree_->Branch("B_k_pz_track", &B_k_pz_track);

  tree_->Branch("B_J_mass", &B_J_mass);
  tree_->Branch("B_J_massErr", &B_J_massErr);
  tree_->Branch("B_J_px", &B_J_px);
  tree_->Branch("B_J_py", &B_J_py);
  tree_->Branch("B_J_pz", &B_J_pz);

  tree_->Branch("B_chi2",    &B_chi2);
  tree_->Branch("B_J_chi2",  &B_J_chi2);*/

  tree_->Branch("DiMu_mu1_pt",  &DiMu_mu1_pt);
  tree_->Branch("DiMu_mu1_eta", &DiMu_mu1_eta);
  tree_->Branch("DiMu_mu1_phi", &DiMu_mu1_phi);
  tree_->Branch("mu1_charge",   &mu1_charge);

  tree_->Branch("DiMu_mu2_pt",  &DiMu_mu2_pt);
  tree_->Branch("DiMu_mu2_eta", &DiMu_mu2_eta);
  tree_->Branch("DiMu_mu2_phi", &DiMu_mu2_phi);
  tree_->Branch("mu2_charge",   &mu2_charge);


  tree_->Branch("DiMu_chi2",    &DiMu_chi2);
  tree_->Branch("DiMu_Prob",  &DiMu_Prob);
       
  tree_->Branch("DiMu_DecayVtxX",     &DiMu_DecayVtxX);
  tree_->Branch("DiMu_DecayVtxY",     &DiMu_DecayVtxY);
  tree_->Branch("DiMu_DecayVtxZ",     &DiMu_DecayVtxZ);
  tree_->Branch("DiMu_DecayVtxXE",    &DiMu_DecayVtxXE);
  tree_->Branch("DiMu_DecayVtxYE",    &DiMu_DecayVtxYE);
  tree_->Branch("DiMu_DecayVtxZE",    &DiMu_DecayVtxZE);
  tree_->Branch("DiMu_DecayVtxXYE",   &DiMu_DecayVtxXYE);
  tree_->Branch("DiMu_DecayVtxXZE",   &DiMu_DecayVtxXZE);
  tree_->Branch("DiMu_DecayVtxYZE",   &DiMu_DecayVtxYZE);

  tree_->Branch("lxy", &lxy);
  tree_->Branch("lxyerr", &lxyerr);

  tree_->Branch("lxy_pv", &lxy_pv);
  tree_->Branch("lxy_pv_err", &lxy_pv_err);

  tree_->Branch("priVtxX",  &priVtxX, "priVtxX/D");
  tree_->Branch("priVtxY",  &priVtxY, "priVtxY/D");
  tree_->Branch("priVtxZ",  &priVtxZ, "priVtxZ/D");
  tree_->Branch("priVtxXE", &priVtxXE, "priVtxXE/D");
  tree_->Branch("priVtxYE", &priVtxYE, "priVtxYE/D");
  tree_->Branch("priVtxZE", &priVtxZE, "priVtxZE/D");
  tree_->Branch("priVtxXYE",&priVtxXYE, "priVtxXYE/D");
  tree_->Branch("priVtxXZE",&priVtxXZE, "priVtxXZE/D");
  tree_->Branch("priVtxYZE",&priVtxYZE, "priVtxYZE/D");
  tree_->Branch("priVtxCL", &priVtxCL, "priVtxCL/D");

  tree_->Branch("nVtx",       &nVtx);
  tree_->Branch("run",        &run,       "run/I");
  tree_->Branch("event",        &event,     "event/L");
  tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");

  tree_muons->Branch("run",      &run,       "run/I");
  tree_muons->Branch("event",    &event,     "event/L");
  tree_muons->Branch("lumiblock",&lumiblock,"lumiblock/I");

  tree_L1muons->Branch("run",      &run,       "run/I");
  tree_L1muons->Branch("event",    &event,     "event/L");
  tree_L1muons->Branch("lumiblock",&lumiblock,"lumiblock/I");

  tree_L2muons->Branch("run",      &run,       "run/I");
  tree_L2muons->Branch("event",    &event,     "event/L");
  tree_L2muons->Branch("lumiblock",&lumiblock,"lumiblock/I");

  tree_L3muons->Branch("run",      &run,       "run/I");
  tree_L3muons->Branch("event",    &event,     "event/L");
  tree_L3muons->Branch("lumiblock",&lumiblock,"lumiblock/I");
    
  // *************************
 
  tree_->Branch("mu1C2",&mu1C2);  
  tree_->Branch("mu1NHits",&mu1NHits);
  tree_->Branch("mu1NPHits",&mu1NPHits);
  tree_->Branch("mu2C2",&mu2C2);  
  tree_->Branch("mu2NHits",&mu2NHits);
  tree_->Branch("mu2NPHits",&mu2NPHits);
  tree_->Branch("mu1dxy",&mu1dxy);
  tree_->Branch("mu2dxy",&mu2dxy);
  tree_->Branch("mu1dz",&mu1dz);
  tree_->Branch("mu2dz",&mu2dz);
  tree_->Branch("muon_dca",&muon_dca);

  tree_->Branch("mu1_L1_match", &mu1_L1_match);
  tree_->Branch("mu1_L2_match", &mu1_L2_match);
  tree_->Branch("mu1_L3_match", &mu1_L3_match);

  tree_->Branch("mu2_L1_match", &mu2_L1_match);
  tree_->Branch("mu2_L2_match", &mu2_L2_match);
  tree_->Branch("mu2_L3_match", &mu2_L3_match);

  tree_->Branch("HLT_Dim25",&HLT_Dim25);
  tree_->Branch("HLT_JpsiTrk_Bc",&HLT_JpsiTrk_Bc);
  tree_->Branch("HLT_JpsiTk",&HLT_JpsiTk);
  tree_->Branch("HLT_DMu4_3_LM",&HLT_DMu4_3_LM);
  tree_->Branch("HLT_DMu4_LM_Displaced",&HLT_DMu4_LM_Displaced);

  tree_L1muons->Branch("L1mu_pt", &L1mu_pt); 
  tree_L1muons->Branch("L1mu_eta", &L1mu_eta);
  tree_L1muons->Branch("L1mu_phi", &L1mu_phi);
  tree_L1muons->Branch("L1mu_etaAtVtx", &L1mu_etaAtVtx);
  tree_L1muons->Branch("L1mu_phiAtVtx", &L1mu_phiAtVtx);
  tree_L1muons->Branch("L1mu_charge", &L1mu_charge);
  tree_L1muons->Branch("L1mu_quality", &L1mu_quality);

  tree_muons->Branch("mu_pt", &mu_pt); 
  tree_muons->Branch("mu_eta", &mu_eta);
  tree_muons->Branch("mu_phi", &mu_phi);
  tree_muons->Branch("mu_charge", &mu_charge);
  
  tree_L2muons->Branch("L2mu_pt", &L2mu_pt); 
  tree_L2muons->Branch("L2mu_eta", &L2mu_eta);
  tree_L2muons->Branch("L2mu_phi", &L2mu_phi);
  
  tree_L3muons->Branch("L3mu_pt", &L3mu_pt); 
  tree_L3muons->Branch("L3mu_eta", &L3mu_eta);
  tree_L3muons->Branch("L3mu_phi", &L3mu_phi);

  tree_->Branch("mu1soft",&mu1soft);
  tree_->Branch("mu2soft",&mu2soft);
  tree_->Branch("mu1tight",&mu1tight);
  tree_->Branch("mu2tight",&mu2tight);
  tree_->Branch("mu1PF",&mu1PF);
  tree_->Branch("mu2PF",&mu2PF);
  tree_->Branch("mu1loose",&mu1loose);
  tree_->Branch("mu2loose",&mu2loose);
  tree_->Branch("mu1Tracker",&mu1Tracker);
  tree_->Branch("mu2Tracker",&mu2Tracker);
  tree_->Branch("mu1Global",&mu1Global);
  tree_->Branch("mu2Global",&mu2Global);

  // gen
  if (isMC_) {
     tree_gen_muons   = fs->make<TTree>("ntuple_gen_muons",  "gen muons ntuple");
     tree_gen_muons->Branch("run",      &run,       "run/I");
     tree_gen_muons->Branch("event",    &event,     "event/L");
     tree_gen_muons->Branch("lumiblock",&lumiblock,"lumiblock/I");
     tree_gen_muons->Branch("GENmu_pt", &GENmu_pt);
     tree_gen_muons->Branch("GENmu_eta", &GENmu_eta);
     tree_gen_muons->Branch("GENmu_phi", &GENmu_phi);
     tree_gen_muons->Branch("GENmu_charge", &GENmu_charge);
     tree_gen_muons->Branch("GENmu_status", &GENmu_status);
     tree_gen_muons->Branch("GENmu_mother", &GENmu_mother);
     tree_gen_muons->Branch("GENmu_grandmother", &GENmu_grandmother);

     tree_->Branch("gen_bc_p4",     "TLorentzVector",  &gen_bc_p4);
     tree_->Branch("gen_jpsi_p4",   "TLorentzVector",  &gen_jpsi_p4);
     tree_->Branch("gen_pion3_p4",  "TLorentzVector",  &gen_pion3_p4);
     tree_->Branch("gen_muon1_p4",  "TLorentzVector",  &gen_muon1_p4);
     tree_->Branch("gen_muon2_p4",  "TLorentzVector",  &gen_muon2_p4);
     tree_->Branch("gen_bc_vtx",    "TVector3",        &gen_bc_vtx);
     tree_->Branch("gen_jpsi_vtx",  "TVector3",        &gen_jpsi_vtx);
     tree_->Branch("gen_bc_ct",     &gen_bc_ct,        "gen_bc_ct/D");

     tree_->Branch("mu1_GEN_match", &mu1_GEN_match);
     tree_->Branch("mu2_GEN_match", &mu2_GEN_match);

  }
  
  
  for(std::size_t i = 0; i < HLTPaths_.size(); ++i){ 
    tree_->Branch(HLTPaths_[i].c_str(), &hltsVector[i]);
    tree_->Branch(("mu1_"+HLTPaths_[i]).c_str(), &mu1_hltsVector[i]);
    tree_->Branch(("mu2_"+HLTPaths_[i]).c_str(), &mu2_hltsVector[i]);
  }
  for(std::size_t i = 0; i < L1Seeds_.size(); ++i){ 
    tree_->Branch(L1Seeds_[i].c_str(), &l1sVector[i]);
  } 

}




// ------------ method called once each job just after ending the event loop  ------------
void MuMu::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuMu);
