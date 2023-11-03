#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuonSetup.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TTree.h"

using namespace edm;
using namespace std;
using namespace reco;

class L1TMuonMiniAODAnalyzerFlat : public edm::one::EDAnalyzer<edm::one::SharedResources> {

public:

  explicit L1TMuonMiniAODAnalyzerFlat(const edm::ParameterSet&);
  ~L1TMuonMiniAODAnalyzerFlat() override;
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:

  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void efficiencies(const edm::Event&, const edm::EventSetup&);
  void purities(const edm::Event&, const edm::EventSetup&);
  void endJob() override;
  virtual void InitandClearStuff();

  // Tokens
  edm::EDGetTokenT<std::vector<Vertex> > verticesToken_;
  edm::EDGetTokenT<std::vector<pat::Muon> > muonToken_;
  edm::EDGetTokenT<GenParticleCollection> genpartToken_;
  edm::EDGetTokenT<l1t::MuonBxCollection> l1MuonToken_;
  
  // Misc
  const PropagateToMuonSetup muPropagatorSetup1st_;
  const PropagateToMuonSetup muPropagatorSetup2nd_;
  PropagateToMuon muPropagator1st_;
  PropagateToMuon muPropagator2nd_;
  
  // Configurables
  Bool_t IsMC_;
  Float_t genMuonPtThreshold_;
  Float_t genMuonEtaThreshold_;
  Float_t recoMuonPtThreshold_;
  Float_t recoMuonEtaThreshold_;
  Float_t trgMuonPtThreshold_;
  Float_t trgMuonEtaThreshold_;
  Float_t genRecoMatchingMaxDR_;
  Bool_t correctMuonCharge_;
  Float_t trgMuonMatchingMaxDR_;
  Bool_t correctTrgCharge_;
  int minTrgQuality_;
  Bool_t purities_;

  // Ntuple
  TTree* outputTree;
  
  // Scalars
  unsigned long eventNb;
  unsigned long runNb;
  unsigned long lumiBlock;
  unsigned long bx;
  int nPV;
  int trueNVtx;
  
  // GEN electrons
  int gen_idx;
  int gen_acc;
  Float_t gen_pt;
  Float_t gen_eta;
  Float_t gen_phi;
  int gen_charge;
  Float_t gen_dxy;
  int gen_pdgid;
  int gen_mom_pdgid;
  int gen_gran_pdgid;
  int gen_muon_idx;
  int gen_muon_match;
  Float_t gen_muon_dr;
  
  // Muons
  Float_t muon_pt;
  Float_t muon_eta;
  Float_t muon_phi;
  int muon_charge;
  Float_t muon_etaAtSt1;
  Float_t muon_etaAtSt2;
  Float_t muon_phiAtSt1;
  Float_t muon_phiAtSt2;
  Float_t muon_dz;
  Float_t muon_dzError;
  Float_t muon_dxy;
  Float_t muon_dxyError;
  Float_t muon_3dIP;
  Float_t muon_3dIPError;
  Bool_t muon_PassTightID;
  Bool_t muon_PassLooseID;
  Bool_t muon_isSAMuon ;
  int muon_trg_idx;
  int muon_trg_match;
  Float_t muon_trg_dr;
  
  // TRG
  Float_t trg_pt;
  Float_t trg_pt_dxy; //@@ KBMTF pT estimate
  Float_t trg_eta;
  Float_t trg_phi;
  int trg_charge;
  Float_t trg_etaAtVtx;
  Float_t trg_phiAtVtx;
  int trg_qual;
  int trg_dxy;
  int trg_tfIdx;
  int trg_bx;

};
