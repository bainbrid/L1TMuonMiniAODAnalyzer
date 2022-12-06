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
#include <fmt/printf.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"


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
using namespace edm;
// using namespace std;
using namespace reco;


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
    bool PassTriggerLeg(std::string triggerlegstring,const pat::Muon *muonit, const edm::Event&);
    virtual bool GetMETFilterDecision(const edm::Event& iEvent, edm::Handle<TriggerResults> METFilterResults, TString studiedfilter);

    // void fillTree();
    // void makeTree();


    // ----------member data ---------------------------
    edm::EDGetTokenT<TriggerResults> metfilterspatToken_;
    edm::EDGetTokenT<TriggerResults> metfiltersrecoToken_;
    edm::EDGetTokenT<std::vector<Vertex> > verticesToken_;

    edm::EDGetTokenT<std::vector< pat::PackedCandidate>> pfcandsToken_;
    edm::EDGetTokenT<std::vector< pat::Muon> > muonToken_;
    edm::EDGetTokenT<std::vector< pat::Muon> > dispMuonToken_;
    edm::EDGetTokenT<GenParticleCollection> genpartToken_;
    edm::EDGetTokenT<std::vector< pat::PackedGenParticle>> packedgenpartToken_;
    edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigobjectToken_;
    // edm::EDGetTokenT<BXVector<GlobalAlgBlk>> l1GtToken_;
    edm::EDGetTokenT<l1t::MuonBxCollection>l1MuonToken_;
    edm::EDGetTokenT<GlobalAlgBlkBxCollection> ugtToken_;
    edm::ESGetToken<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd> l1GtMenuToken_;

    Float_t MuonPtCut_;
    Bool_t SaveTree_, IsMC_, Debug_;

    const PropagateToMuonSetup muPropagatorSetup1st_;
    const PropagateToMuonSetup muPropagatorSetup2nd_;

    PropagateToMuon muPropagator1st_;
    PropagateToMuon muPropagator2nd_;

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

    GlobalAlgBlk const *results_;
    unsigned long long cache_id_;

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

    //Trigger matching variables: is a reco object matched to a trigger filter
    vector<Bool_t> muon_isIsoHLTMuon;
    vector<Bool_t> muon_isHLTMuon;
    vector<Bool_t> muon_isTauTo3MuMuon;
    vector<Bool_t> muon_isDoubleMuForBsMuon;
    vector<Bool_t> muon_isDoubleMuLowMassMuon;
    vector<Bool_t> muon_isSingleMuMuon;

    int muon_size;

    // Displaced Muons
    vector<Float_t>  dispMuon_eta;
    vector<Float_t>  dispMuon_etaAtSt1;
    vector<Float_t>  dispMuon_etaAtSt2;
    vector<Float_t>  dispMuon_phi;
    vector<Float_t>  dispMuon_phiAtSt1;
    vector<Float_t>  dispMuon_phiAtSt2;
    vector<Float_t>  dispMuon_pt;
    vector<Float_t>  dispMuon_ptCorr;
    vector <int>     dispMuon_charge;

    vector<Float_t>  dispMuon_dz;
    vector<Float_t>  dispMuon_dzError;
    vector<Float_t>  dispMuon_dxy;
    vector<Float_t>  dispMuon_dxyError;
    vector<Float_t>  dispMuon_3dIP;
    vector<Float_t>  dispMuon_3dIPError;

    vector<Bool_t>  dispMuon_PassTightID;
    vector<Bool_t>  dispMuon_PassLooseID;
    vector<Bool_t> dispMuon_isSAMuon ;

    //Trigger matching variables: is a reco object matched to a trigger filter
    vector<Bool_t> dispMuon_isDispL2HLTDiMuon;
    vector<Bool_t> dispMuon_isDispHLTDiMuon;
    vector<Bool_t> dispMuon_isL2DispHLTMuon;
    vector<Bool_t> dispMuon_isNoBPTXHLTMuon;
    vector<Bool_t> dispMuon_isL2NoVtxHLTMuon;
    vector<Bool_t> dispMuon_isL2CosmicSeedHLTMuon;
    vector<Bool_t> dispMuon_isL1NoBPTXMuon;
    vector<Bool_t> dispMuon_isL1DispMuon;

    int dispMuon_size;

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
    bool HLT_DoubleMu43NoFiltersNoVtx;
    bool HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm;
    bool HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm;
    bool HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed;
    bool HLT_DoubleL2Mu23NoVtx_2Cha;
    bool HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm;
    bool HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm;
    bool HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm;
    bool HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX;
    bool HLT_DoubleMu4_3_Bs;
    bool HLT_DoubleMu4_3_Jpsi;
    bool HLT_DoubleMu4_3_Photon4_BsToMMG;
    bool HLT_DoubleMu4_LowMass_Displaced;
    bool HLT_DoubleMu4_3_LowMass;
    bool HLT_DoubleMu3_TkMu_DsTau3Mu;
    bool HLT_DoubleMu3_Trk_Tau3mu;

    //L1 muon
    vector <int> l1mu_qual;
    vector <int> l1mu_charge;
    vector <Float_t> l1mu_pt;
    vector <Float_t> l1mu_pt_dxy;
    vector <int> l1mu_dxy;
    vector <Float_t> l1mu_eta;
    vector <Float_t> l1mu_etaAtVtx;
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

};