#include "L1TMuonMiniAODAnalyzerFlat.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

////////////////////////////////////////////////////////////////////////////////
//
L1TMuonMiniAODAnalyzerFlat::L1TMuonMiniAODAnalyzerFlat(const edm::ParameterSet& iConfig)
  :
  verticesToken_(consumes<std::vector<Vertex> >(iConfig.getParameter<edm::InputTag>("Vertices"))),
  muonToken_(consumes< std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("Muons"))),
  genpartToken_(consumes<GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticles"))),
  l1MuonToken_(consumes<l1t::MuonBxCollection>(iConfig.getParameter<edm::InputTag>("L1T"))),
  muPropagatorSetup1st_(iConfig.getParameter<edm::ParameterSet>("muProp1st"), consumesCollector()),
  muPropagatorSetup2nd_(iConfig.getParameter<edm::ParameterSet>("muProp2nd"), consumesCollector()),
  IsMC_(iConfig.getParameter<bool>("IsMC")),
  genMuonPtThreshold_(0.),
  genMuonEtaThreshold_(10.),
  recoMuonPtThreshold_(0.),
  recoMuonEtaThreshold_(10.),
  recoMatchingMaxDR_(10.),
  correctMuonCharge_(true),
  trgMatchingMaxDR_(10.),
  correctTrgCharge_(true),
  minTrgQuality_(0)
{

  edm::Service<TFileService> fs;
  outputTree = fs->make<TTree>("tree","tree");

  if ( iConfig.exists("genMuonPtThreshold") && iConfig.exists("genMuonEtaThreshold") ) {
    genMuonPtThreshold_ = iConfig.getParameter<double>("genMuonPtThreshold");
    genMuonEtaThreshold_ = iConfig.getParameter<double>("genMuonEtaThreshold");
  }

  if ( iConfig.exists("recoMuonPtThreshold") && iConfig.exists("recoMuonEtaThreshold") ) {
    recoMuonPtThreshold_ = iConfig.getParameter<double>("recoMuonPtThreshold");
    recoMuonEtaThreshold_ = iConfig.getParameter<double>("recoMuonEtaThreshold");
  }

  if ( iConfig.exists("recoMatchingMaxDR") ) {
    recoMatchingMaxDR_ = iConfig.getParameter<double>("recoMatchingMaxDR");
  }

  if ( iConfig.exists("correctMuonCharge") ) {
    correctMuonCharge_ = iConfig.getParameter<bool>("correctMuonCharge");
  }

  if ( iConfig.exists("trgMatchingMaxDR") ) {
    trgMatchingMaxDR_ = iConfig.getParameter<double>("trgMatchingMaxDR");
  }

  if ( iConfig.exists("correctTrgCharge") ) {
    correctTrgCharge_ = iConfig.getParameter<bool>("correctTrgCharge");
  }

  if ( iConfig.exists("minTrgQuality") ) {
    minTrgQuality_ = iConfig.getParameter<int>("minTrgQuality");
  }
  
}

////////////////////////////////////////////////////////////////////////////////
//
L1TMuonMiniAODAnalyzerFlat::~L1TMuonMiniAODAnalyzerFlat() {}
//

////////////////////////////////////////////////////////////////////////////////
//
void L1TMuonMiniAODAnalyzerFlat::analyze(const edm::Event& iEvent,
					 const edm::EventSetup& iSetup) {
  
  InitandClearStuff();

  muPropagator1st_ = muPropagatorSetup1st_.init(iSetup);
  muPropagator2nd_ = muPropagatorSetup2nd_.init(iSetup);

  runNb = iEvent.id().run();
  eventNb = iEvent.id().event();
  lumiBlock = iEvent.luminosityBlock();
  bx = iEvent.bunchCrossing();

  edm::Handle<std::vector<Vertex> > theVertices;
  iEvent.getByToken(verticesToken_,theVertices) ;
  nPV = theVertices->size();
  Vertex::Point PV(0,0,0);
  if(nPV){ PV = theVertices->begin()->position();}

  edm::Handle< std::vector<pat::Muon> > thePatMuons;
  iEvent.getByToken(muonToken_,thePatMuons);
  
  edm::Handle<l1t::MuonBxCollection> l1muoncoll;
  iEvent.getByToken(l1MuonToken_ , l1muoncoll);
  
  if (IsMC_) {

    std::vector<int> matched_muon;
    std::vector<int> matched_trg;
    
    // Iterate through GEN particles
    edm::Handle< std::vector<reco::GenParticle> > genParticles;
    iEvent.getByToken(genpartToken_,genParticles);
    for ( size_t igen = 0; igen < genParticles->size(); igen++ ) {
      reco::GenParticlePtr gen(genParticles, igen);
      if ( gen.isNull() || !(gen.isAvailable()) ) { 
	std::cerr << "[L1TMuonMiniAODAnalyzerFlat::analyze] ERROR! reco::GenParticlePtr:"
		  << " gen.isNull(): " << gen.isNull()
		  << " gen.isAvailable(): " << gen.isAvailable()
		  << std::endl;
	continue;
      }
      
      // Is GEN particle a muon?
      bool is_gen_muon = 
	( std::abs(gen->pdgId()) == 13 )
	&& ( gen->isLastCopy() );
      if ( !is_gen_muon ) { continue; }
      
      // Store 
      gen_idx = (int)igen;
      gen_pt = gen->pt();
      gen_eta = gen->eta();
      gen_phi = gen->phi();
      gen_charge = gen->charge();
      //gen_dxy = gen->dxy();
      gen_pdgid = gen->pdgId();
      //gen_mom_pdgid = gen->pt();
      //gen_gran_pdgid = gen->pt();
      
      // Is GEN particle in acceptance?
      bool acc_gen = 
	( gen->pt() > genMuonPtThreshold_ )
	&& ( std::fabs(gen->eta()) < genMuonEtaThreshold_ );
      gen_acc = acc_gen; //@@ continue below prevents Fill
      if (!acc_gen) { continue; }

      // Loop through RECO muons
      int temp_idx = -1;
      float temp_dr = recoMatchingMaxDR_;
      unsigned int muon_size = (unsigned int)thePatMuons->size();
      for( unsigned int imuon = 0; imuon < muon_size; ++imuon ) {
	edm::Ptr<pat::Muon> muon(thePatMuons,imuon);
	
	// Only consider RECO muons in acceptance
	bool acc_reco = 
	  ( muon->pt() > recoMuonPtThreshold_ )
	  && ( std::fabs(muon->eta()) < recoMuonEtaThreshold_ );
	if (acc_reco) {
	  
	  // Only consider RECO muons passing Loose ID
	  bool loose_id = muon->passed(reco::Muon::CutBasedIdLoose);
	  if (loose_id) {
	    
	    // If already matched to a GEN muon, ignore
	    if ( std::find(matched_muon.begin(), matched_muon.end(), imuon) != matched_muon.end() ) { continue; }
	    
	    float dr = reco::deltaR( gen->eta(), gen->phi(), muon->eta(), muon->phi() );
	    int dq = gen->charge() * muon->charge();
	    if ( dr < temp_dr && (dq==1 || !correctMuonCharge_) ) {
	      //@@ matched!
	      temp_dr = dr;
	      temp_idx = (int)imuon;
	    }
	  }
	}
      }
      
      // If matched, continue
      if ( temp_idx >= 0 ) {
	
	// Store matched result
	matched_muon.push_back(temp_idx);
	reco_idx = temp_idx;
	reco_match = (temp_idx>=0);
	reco_dr = temp_dr;
	
	// Ptr to matched muon
	edm::Ptr<pat::Muon> reco(thePatMuons,reco_idx);
	
	// Store muons
	reco_pt = reco->pt();
	reco_eta = reco->eta();
	reco_phi = reco->phi();
	reco_charge = reco->charge();
	
	reco_PassLooseID = reco->passed(reco::Muon::CutBasedIdLoose);
	reco_PassTightID = reco->passed(reco::Muon::CutBasedIdMediumPrompt);
	reco_isSAMuon = reco->isStandAloneMuon();
	
	if ( !(reco->innerTrack()).isNull() ) {
	  reco_dz = reco->innerTrack()->dz(PV);
	  reco_dzError = reco->innerTrack()->dzError();
	  reco_dxy = reco->innerTrack()->dxy(PV);
	  reco_dxyError = reco->innerTrack()->dxyError();
	  reco_3dIP = reco->dB(pat::Muon::PV3D);
	  reco_3dIPError = reco->edB(pat::Muon::PV3D);
	}
	
	//      // extrapolation of muon track coordinates
	//      TrajectoryStateOnSurface stateAtMuSt1 = muPropagator1st_.extrapolate(*reco);
	//      if (stateAtMuSt1.isValid()) {
	//        reco_etaAtSt1 = stateAtMuSt1.globalPosition().eta();
	//        reco_phiAtSt1 = stateAtMuSt1.globalPosition().phi();
	//      }
	//      
	//      TrajectoryStateOnSurface stateAtMuSt2 = muPropagator2nd_.extrapolate(*reco);
	//      if (stateAtMuSt2.isValid()) {
	//        reco_etaAtSt2 = stateAtMuSt2.globalPosition().eta();
	//        reco_phiAtSt2 = stateAtMuSt2.globalPosition().phi();
	//      }

	int tmp_idx = -1;
	float tmp_dr = trgMatchingMaxDR_;
	l1t::MuonBxCollection::const_iterator tmp;
	for ( int ibx = l1muoncoll->getFirstBX(); ibx <= l1muoncoll->getLastBX(); ++ibx ) {
	  l1t::MuonBxCollection::const_iterator iter = l1muoncoll->begin(ibx);
	  l1t::MuonBxCollection::const_iterator end = l1muoncoll->end(ibx);
	  for( ; iter != end; ++iter ) {
	    if (iter->pt() < 0) { continue; } // Check if correct BX?
	    int itrg = ibx*100 + std::distance(iter,end);

	    // If already matched to a muon, ignore
	    if ( std::find(matched_trg.begin(), matched_trg.end(), itrg) == matched_trg.end() ) { 
	      float dr = reco::deltaR( reco->eta(), reco->phi(), iter->etaAtVtx(), iter->phiAtVtx() );
	      int dq = iter->charge() * reco->charge();
	      if ( (dr<tmp_dr) && 
		   (dq==1 || !correctTrgCharge_) && 
		   (iter->hwQual()>=minTrgQuality_) ) {
		//@@ matched!
		tmp_dr = dr;
		tmp_idx = (int)itrg;
		tmp = iter;
	      }
	    }
	  }
	}

	// If matched, continue
	if ( tmp_idx >= 0 ) {
	    
	  // Store matched result
	  matched_trg.push_back(tmp_idx);
	  trg_idx = tmp_idx;
	  trg_match = (tmp_idx>=0);
	  trg_dr = tmp_dr;

	  // Store TRG
	  trg_pt = tmp->pt();
	  trg_pt_dxy = tmp->ptUnconstrained();
	  trg_eta = tmp->eta();
	  trg_phi = tmp->phi();
	  trg_charge = tmp->charge();
	  trg_etaAtVtx = tmp->etaAtVtx();
	  trg_phiAtVtx = tmp->phiAtVtx();
	  trg_qual = tmp->hwQual();
	  trg_dxy = tmp->hwDXY();
	  trg_tfIdx = tmp->tfMuonIndex();
	  //trg_bx = ibx;

	} // TRG matched to reco?
      } // GEN matched to reco?

      outputTree->Fill();
      
    } // Loop over GEN particles
    
  } // IsMC?
  
}

////////////////////////////////////////////////////////////////////////////////
//
void L1TMuonMiniAODAnalyzerFlat::beginJob() {

  // Scalars
  outputTree->Branch("event",&eventNb,"event/l");
  outputTree->Branch("run",&runNb,"run/l");
  outputTree->Branch("ls",&lumiBlock,"ls/l");
  outputTree->Branch("bx",&bx,"bx/l");
  outputTree->Branch("nPV",&nPV,"nPV/I");
  if (IsMC_) { outputTree->Branch("trueNVtx", &trueNVtx,"trueNVtx/I"); }

  // GEN
  outputTree->Branch("gen_idx",&gen_idx,"gen_idx/I");
  outputTree->Branch("gen_acc",&gen_acc,"gen_acc/I");
  outputTree->Branch("gen_pt",&gen_pt,"gen_pt/F");
  outputTree->Branch("gen_eta",&gen_eta,"gen_eta/F");
  outputTree->Branch("gen_phi",&gen_phi,"gen_phi/F");
  outputTree->Branch("gen_charge",&gen_charge,"gen_charge/I");
  outputTree->Branch("gen_dxy", &gen_dxy,"gen_dxy/F");
  outputTree->Branch("gen_pdgid",&gen_pdgid,"gen_pdgid/I");
  //outputTree->Branch("gen_mom_pdgid",&gen_mom_pdgid);
  //outputTree->Branch("gen_gran_pdgid",&gen_gran_pdgid);
  outputTree->Branch("reco_idx",&reco_idx,"reco_idx/I");
  outputTree->Branch("reco_match",&reco_match,"reco_match/I");
  outputTree->Branch("reco_dr",&reco_dr,"reco_dr/F");

  // RECO
  outputTree->Branch("reco_pt",&reco_pt);
  outputTree->Branch("reco_eta",&reco_eta);
  outputTree->Branch("reco_phi",&reco_phi);
  outputTree->Branch("reco_charge",&reco_charge);
  outputTree->Branch("reco_etaAtSt1",&reco_etaAtSt1);
  outputTree->Branch("reco_etaAtSt2",&reco_etaAtSt2);
  outputTree->Branch("reco_phiAtSt1",&reco_phiAtSt1);
  outputTree->Branch("reco_phiAtSt2",&reco_phiAtSt2);
  outputTree->Branch("reco_dz",&reco_dz);
  outputTree->Branch("reco_dzError",&reco_dzError);
  outputTree->Branch("reco_dxy",&reco_dxy);
  outputTree->Branch("reco_dxyError",&reco_dxyError);
  outputTree->Branch("reco_3dIP",&reco_3dIP);
  outputTree->Branch("reco_3dIPError",&reco_3dIPError);
  outputTree->Branch("reco_PassTightID",&reco_PassTightID);
  outputTree->Branch("reco_PassLooseID",&reco_PassLooseID);
  outputTree->Branch("reco_isSAMuon",&reco_isSAMuon);
  outputTree->Branch("trg_idx",&trg_idx,"trg_idx/I");
  outputTree->Branch("trg_match",&trg_match,"trg_match/I");
  outputTree->Branch("trg_dr",&trg_dr,"trg_dr/F");

  // TRG
  outputTree->Branch("trg_pt",&trg_pt);
  outputTree->Branch("trg_pt_dxy",&trg_pt_dxy); // KBMTF pT estimate
  outputTree->Branch("trg_eta",&trg_eta);
  outputTree->Branch("trg_phi",&trg_phi);
  outputTree->Branch("trg_charge",&trg_charge);
  outputTree->Branch("trg_etaAtVtx",&trg_etaAtVtx);
  outputTree->Branch("trg_phiAtVtx",&trg_phiAtVtx);
  outputTree->Branch("trg_qual",&trg_qual);
  outputTree->Branch("trg_dxy",&trg_dxy);
  outputTree->Branch("trg_tfIdx",&trg_tfIdx);
  outputTree->Branch("trg_bx",&trg_bx);

}

////////////////////////////////////////////////////////////////////////////////
//
void L1TMuonMiniAODAnalyzerFlat::endJob() {}

////////////////////////////////////////////////////////////////////////////////
//
void L1TMuonMiniAODAnalyzerFlat::InitandClearStuff() {

  // Scalars
  nPV = 0;
  trueNVtx = 0;

  // GEN
  gen_idx = 0;
  gen_acc = 0;
  gen_pt = 0.;
  gen_eta = 0.;
  gen_phi = 0.;
  gen_charge = 0;
  gen_dxy = 0.;
  gen_pdgid = 0;
  gen_mom_pdgid = 0;
  gen_gran_pdgid = 0;
  reco_idx = -1;
  reco_match = 0;
  reco_dr = 10.;
  
  // RECO
  reco_pt = 0.;
  reco_eta = 0.;
  reco_phi = 0.;
  reco_charge = 0;
  reco_etaAtSt1 = 0.;
  reco_etaAtSt2 = 0.;
  reco_phiAtSt1 = 0.;
  reco_phiAtSt2 = 0.;
  reco_dz = 0.;
  reco_dzError = 0.;
  reco_dxy = 0.;
  reco_dxyError = 0.;
  reco_3dIP = 0.;
  reco_3dIPError = 0.;
  reco_PassTightID = 0;
  reco_PassLooseID = 0;
  reco_isSAMuon = 0;
  trg_idx = -1;
  trg_match = 0;
  trg_dr = 10.;
  
  // L1 TRG
  trg_pt = 0.;
  trg_pt_dxy = 0.;
  trg_eta = 0.;
  trg_phi = 0.;
  trg_charge = 0;
  trg_etaAtVtx = 0.;
  trg_phiAtVtx = 0.;
  trg_qual = 0;
  trg_dxy = 0;
  trg_tfIdx = 0;
  trg_bx = 0;

}

void L1TMuonMiniAODAnalyzerFlat::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// Define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TMuonMiniAODAnalyzerFlat);
