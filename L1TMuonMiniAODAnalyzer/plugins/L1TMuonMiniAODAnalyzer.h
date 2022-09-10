// -*- C++ -*-
//
// Package:    L1TMuonMiniAODAnalyzer/L1TMuonMiniAODAnalyzer
// Class:      L1TMuonMiniAODAnalyzer
//
/**\class L1TMuonMiniAODAnalyzer L1TMuonMiniAODAnalyzer.cc L1TMuonMiniAODAnalyzer/L1TMuonMiniAODAnalyzer/plugins/L1TMuonMiniAODAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Efe Yigitbasi
//         Created:  Sat, 10 Sep 2022 11:08:53 GMT
//
//

// system include files
#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include "TLorentzVector.h"
#include "TFile.h"
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

// Data formats
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"

// muon track extrapolation
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuonSetup.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


class L1TMuonMiniAODAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit L1TMuonMiniAODAnalyzer(const edm::ParameterSet&);
  ~L1TMuonMiniAODAnalyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  // void beginRun(const edm::Run&, const edm::EventSetup&);
  // void endRun(const edm::Run&, const edm::EventSetup&);
  virtual void InitandClearStuff();
  bool PassTriggerLeg(std::string triggerlegstring, std::string triggerlegstringalt,const pat::Muon *muonit, const edm::Event&);

  // void fillTree();
  // void makeTree();


  // ----------member data ---------------------------
  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

  edm::EDGetTokenT<TriggerResults> metfilterspatToken_;
  edm::EDGetTokenT<TriggerResults> metfiltersrecoToken_;
  edm::EDGetTokenT<std::vector<Vertex> > verticesToken_;

  edm::EDGetTokenT<std::vector< pat::PackedCandidate>> pfcandsToken_;
  edm::EDGetTokenT<std::vector< pat::Muon> > muonToken_;
  edm::EDGetTokenT<GenParticleCollection> genpartToken_;
  edm::EDGetTokenT<std::vector< pat::PackedGenParticle>> packedgenpartToken_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigobjectToken_;
  edm::EDGetTokenT<BXVector<GlobalAlgBlk>> l1GtToken_;
  edm::EDGetTokenT<l1t::MuonBxCollection>l1MuonToken_;

  const PropagateToMuonSetup muPropagatorSetup1st_;
  const PropagateToMuonSetup muPropagatorSetup2nd_;

  PropagateToMuon muPropagator1st_;
  PropagateToMuon muPropagator2nd_;



  Float_t MuonPtCut_;
  Bool_t SaveTree_, IsMC_, Debug_;

  TTree* outputTree;

  unsigned long _eventNb;
  unsigned long _runNb;
  unsigned long _lumiBlock;
  unsigned long _bx;

  //Nb of primary vertices
  int _n_PV;
  // Float_t _LV_x,_LV_y,_LV_z;
  // Float_t _LV_errx,_LV_erry,_LV_errz;
  // Float_t _PUV1_x,_PUV1_y,_PUV1_z;
  int trueNVtx;

  //MINIAOD original MET filters decisions
  bool Flag_goodVertices;
  bool Flag_globalTightHalo2016Filter;
  bool Flag_globalSuperTightHalo2016Filter;
  bool Flag_BadPFMuonFilter;
  bool Flag_BadPFMuonDzFilter;

  //Muons
  vector<Float_t>  muon_eta;
  vector<Float_t>  muon_etaAtSt1;
  vector<Float_t>  muon_etaAtSt2;
  vector<Float_t>  muon_phi;
  vector<Float_t>  muon_phiAtSt1;
  vector<Float_t>  muon_phiAtSt2;
  vector<Float_t>  muon_pt;
  vector<Float_t>  muon_ptCorr;
  vector <int>     muon_charge;

  vector<Float_t>  muon_dz;
  vector<Float_t>  muon_dzError;
  vector<Float_t>  muon_dxy;
  vector<Float_t>  muon_dxyError;
  vector<Float_t>  muon_3dIP;
  vector<Float_t>  muon_3dIPError;

  vector<Bool_t>  muon_PassTightID;
  vector<Bool_t>  muon_PassLooseID;
  vector<Bool_t> muon_isSAMuon ;
  int muon_size;

  //Event variables (reco)
  //For dilepton studies
  Float_t _mll;
  Float_t _ptll;
  Float_t _pzll;
  Float_t _yll;
  Float_t _dphill;
  Float_t _phill;
  Float_t _costhCSll;

  //Triggers
  bool HLT_IsoMu27;
  bool HLT_IsoMu24;
  bool HLT_IsoTkMu24;
  bool HLT_Mu50;

  //Trigger matching variables: is a reco object matched to a trigger filter
  vector < bool >hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered;


  //L1 muon
  vector <int> l1mu_qual;
  vector <int> l1mu_charge;
  vector <Float_t> l1mu_pt;
  vector <Float_t> l1mu_pt_dxy;
  vector <int> l1mu_dxy;
  vector <Float_t> l11mu_eta;
  vector <Float_t> l11mu_etaAtVtx;
  vector <Float_t> l1mu_phi;
  vector <Float_t> l1mu_phiAtVtx;
  vector <int> l1mu_tfIdx;
  vector <int> l1mu_bx;
  int l1mu_size;

//
// constants, enums and typedefs
//

const int  N_METFilters=18;
enum METFilterIndex{
  idx_Flag_goodVertices,
  idx_Flag_globalTightHalo2016Filter,
  idx_Flag_globalSuperTightHalo2016Filter,
  idx_Flag_HBHENoiseFilter,
  idx_Flag_HBHENoiseIsoFilter,
  idx_Flag_EcalDeadCellTriggerPrimitiveFilter,
  idx_Flag_BadPFMuonFilter,
  idx_Flag_BadPFMuonDzFilter,
  idx_Flag_hfNoisyHitsFilter,
  idx_Flag_BadChargedCandidateFilter,
  idx_Flag_eeBadScFilter,
  idx_Flag_ecalBadCalibFilter,
  idx_Flag_ecalLaserCorrFilter,
  idx_Flag_EcalDeadCellBoundaryEnergyFilter,
  idx_PassecalBadCalibFilter_Update,
  idx_PassecalLaserCorrFilter_Update,
  idx_PassEcalDeadCellBoundaryEnergyFilter_Update,
  idx_PassBadChargedCandidateFilter_Update
};


//
// static data member definitions
//

