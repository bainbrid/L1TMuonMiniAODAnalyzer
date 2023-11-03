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
  trgMuonPtThreshold_(0.),
  trgMuonEtaThreshold_(10.),
  genRecoMatchingMaxDR_(10.),
  correctMuonCharge_(true),
  trgMuonMatchingMaxDR_(10.),
  correctTrgCharge_(true),
  minTrgQuality_(0),
  purities_(0)
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

  if ( iConfig.exists("trgMuonPtThreshold") && iConfig.exists("trgMuonEtaThreshold") ) {
    trgMuonPtThreshold_ = iConfig.getParameter<double>("trgMuonPtThreshold");
    trgMuonEtaThreshold_ = iConfig.getParameter<double>("trgMuonEtaThreshold");
  }

  if ( iConfig.exists("genRecoMatchingMaxDR") ) {
    genRecoMatchingMaxDR_ = iConfig.getParameter<double>("genRecoMatchingMaxDR");
  }

  if ( iConfig.exists("correctMuonCharge") ) {
    correctMuonCharge_ = iConfig.getParameter<bool>("correctMuonCharge");
  }

  if ( iConfig.exists("trgMatchingMaxDR") ) {
    trgMuonMatchingMaxDR_ = iConfig.getParameter<double>("trgMatchingMaxDR");
  }

  if ( iConfig.exists("correctTrgCharge") ) {
    correctTrgCharge_ = iConfig.getParameter<bool>("correctTrgCharge");
  }

  if ( iConfig.exists("minTrgQuality") ) {
    minTrgQuality_ = iConfig.getParameter<int>("minTrgQuality");
  }
  
  if ( iConfig.exists("calcPurity") ) {
    purities_ = iConfig.getParameter<bool>("calcPurity");
  }
  std::cout << "[L1TMuonMiniAODAnalyzerFlat::L1TMuonMiniAODAnalyzerFlat]"
	    << " purities_ = " << purities_
	    << std::endl;
  
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
  if (purities_) { purities(iEvent,iSetup); }
  else           { efficiencies(iEvent,iSetup); }
}

////////////////////////////////////////////////////////////////////////////////
//
void L1TMuonMiniAODAnalyzerFlat::efficiencies(const edm::Event& iEvent,
					      const edm::EventSetup& iSetup) {
  
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
  
  if (!IsMC_) { return; }

  edm::Handle< std::vector<reco::GenParticle> > genParticles;
  iEvent.getByToken(genpartToken_,genParticles);

  std::vector<int> matched_muon;
  std::vector<int> matched_trg;
  
  // Iterate through GEN particles
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
    float temp_dr = genRecoMatchingMaxDR_;
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
      gen_muon_idx = temp_idx;
      gen_muon_match = (gen_muon_idx>=0);
      gen_muon_dr = temp_dr;

      // Record matched result
      matched_muon.push_back(gen_muon_idx);
      
      // Ptr to matched muon
      edm::Ptr<pat::Muon> reco(thePatMuons,gen_muon_idx);
      
      // Store muons
      muon_pt = reco->pt();
      muon_eta = reco->eta();
      muon_phi = reco->phi();
      muon_charge = reco->charge();
      
      muon_PassLooseID = reco->passed(reco::Muon::CutBasedIdLoose);
      muon_PassTightID = reco->passed(reco::Muon::CutBasedIdMediumPrompt);
      muon_isSAMuon = reco->isStandAloneMuon();
      
      if ( !(reco->innerTrack()).isNull() ) {
	muon_dz = reco->innerTrack()->dz(PV);
	muon_dzError = reco->innerTrack()->dzError();
	muon_dxy = reco->innerTrack()->dxy(PV);
	muon_dxyError = reco->innerTrack()->dxyError();
	muon_3dIP = reco->dB(pat::Muon::PV3D);
	muon_3dIPError = reco->edB(pat::Muon::PV3D);
      }
      
      //      // extrapolation of muon track coordinates
      //      TrajectoryStateOnSurface stateAtMuSt1 = muPropagator1st_.extrapolate(*reco);
      //      if (stateAtMuSt1.isValid()) {
      //        muon_etaAtSt1 = stateAtMuSt1.globalPosition().eta();
      //        muon_phiAtSt1 = stateAtMuSt1.globalPosition().phi();
      //      }
      //      
      //      TrajectoryStateOnSurface stateAtMuSt2 = muPropagator2nd_.extrapolate(*reco);
      //      if (stateAtMuSt2.isValid()) {
      //        muon_etaAtSt2 = stateAtMuSt2.globalPosition().eta();
      //        muon_phiAtSt2 = stateAtMuSt2.globalPosition().phi();
      //      }
      
      // Iterate through trigger muons
      int tmp_idx = -1;
      float tmp_dr = trgMuonMatchingMaxDR_;
      l1t::MuonBxCollection::const_iterator tmp;
      l1t::MuonBxCollection::const_iterator iter = l1muoncoll->begin(0);
      l1t::MuonBxCollection::const_iterator end = l1muoncoll->end(0);
      for( ; iter != end; ++iter ) {
	int itrg = std::distance(iter,end);
	
	// If already matched to a muon, ignore
	if ( std::find(matched_trg.begin(), matched_trg.end(), itrg) != matched_trg.end() ) { continue; }
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
      
      // If matched, continue
      if ( tmp_idx >= 0 ) {

	// Store matched result
	muon_trg_idx = tmp_idx;
	muon_trg_match = (muon_trg_idx>=0);
	muon_trg_dr = tmp_dr;

	// Record matched result
	matched_trg.push_back(muon_trg_idx);
	
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
	
      } // TRG matched to reco?
    } // GEN matched to reco?
    
    outputTree->Fill();
    
  } // Loop over GEN particles
  
}

////////////////////////////////////////////////////////////////////////////////
//
void L1TMuonMiniAODAnalyzerFlat::purities(const edm::Event& iEvent,
					  const edm::EventSetup& iSetup) {
  
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

  if (!IsMC_) { return; }

  edm::Handle< std::vector<reco::GenParticle> > genParticles;
  iEvent.getByToken(genpartToken_,genParticles);

  std::vector<int> matched_muon;
  std::vector<int> matched_trg;

  // Iterate through trigger objects
  l1t::MuonBxCollection::const_iterator iter = l1muoncoll->begin(0);
  l1t::MuonBxCollection::const_iterator end = l1muoncoll->end(0);
  for( ; iter != end; ++iter ) {
    //int itrg = std::distance(iter,end);

    // Only consider TRG muons in acceptance and with HW quality
    bool acc_trg = 
      ( iter->pt() > trgMuonPtThreshold_ )
      && ( std::fabs(iter->eta()) < trgMuonEtaThreshold_ )
      && ( iter->hwQual()>=minTrgQuality_ );
    if (acc_trg) {

      // Store TRG
      trg_pt = iter->pt();
      trg_pt_dxy = iter->ptUnconstrained();
      trg_eta = iter->eta();
      trg_phi = iter->phi();
      trg_charge = iter->charge();
      trg_etaAtVtx = iter->etaAtVtx();
      trg_phiAtVtx = iter->phiAtVtx();
      trg_qual = iter->hwQual();
      trg_dxy = iter->hwDXY();
      trg_tfIdx = iter->tfMuonIndex();
      
      // Loop through RECO muons
      int tmp_idx = -1;
      float tmp_dr = trgMuonMatchingMaxDR_;
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
	    
	    // If already matched to a muon, ignore
	    if ( std::find(matched_trg.begin(), matched_trg.end(), imuon) != matched_trg.end() ) { continue; }
	    float dr = reco::deltaR( muon->eta(), muon->phi(), iter->etaAtVtx(), iter->phiAtVtx() );
	    int dq = iter->charge() * muon->charge();
	    if ( (dr<tmp_dr) && (dq==1 || !correctTrgCharge_) ) {
	      //@@ matched!
	      tmp_dr = dr;
	      tmp_idx = (int)imuon;
	    }
	  }
	}
      }
	
      // If matched, continue
      if ( tmp_idx >= 0 ) {
	
	// Store matched result
	muon_trg_idx = tmp_idx;
	muon_trg_match = (muon_trg_idx>=0);
	muon_trg_dr = tmp_dr;

	// Record matched result
	matched_trg.push_back(muon_trg_idx);
	
	// Ptr to matched muon
	edm::Ptr<pat::Muon> muon(thePatMuons,muon_trg_idx);
	
	// Store muons
	muon_pt = muon->pt();
	muon_eta = muon->eta();
	muon_phi = muon->phi();
	muon_charge = muon->charge();
	
	muon_PassLooseID = muon->passed(reco::Muon::CutBasedIdLoose);
	muon_PassTightID = muon->passed(reco::Muon::CutBasedIdMediumPrompt);
	muon_isSAMuon = muon->isStandAloneMuon();
	
	if ( !(muon->innerTrack()).isNull() ) {
	  muon_dz = muon->innerTrack()->dz(PV);
	  muon_dzError = muon->innerTrack()->dzError();
	  muon_dxy = muon->innerTrack()->dxy(PV);
	  muon_dxyError = muon->innerTrack()->dxyError();
	  muon_3dIP = muon->dB(pat::Muon::PV3D);
	  muon_3dIPError = muon->edB(pat::Muon::PV3D);
	}
	
	//      // extrapolation of muon track coordinates
	//      TrajectoryStateOnSurface stateAtMuSt1 = muPropagator1st_.extrapolate(*reco);
	//      if (stateAtMuSt1.isValid()) {
	//        muon_etaAtSt1 = stateAtMuSt1.globalPosition().eta();
	//        muon_phiAtSt1 = stateAtMuSt1.globalPosition().phi();
	//      }
	//      
	//      TrajectoryStateOnSurface stateAtMuSt2 = muPropagator2nd_.extrapolate(*reco);
	//      if (stateAtMuSt2.isValid()) {
	//        muon_etaAtSt2 = stateAtMuSt2.globalPosition().eta();
	//        muon_phiAtSt2 = stateAtMuSt2.globalPosition().phi();
	//      }

	// Iterate through GEN particles
	int temp_idx = -1;
	float temp_dr = genRecoMatchingMaxDR_;
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
	  
	  // Is GEN particle in acceptance?
	  bool acc_gen = 
	    ( gen->pt() > genMuonPtThreshold_ )
	    && ( std::fabs(gen->eta()) < genMuonEtaThreshold_ );
	  gen_acc = acc_gen; //@@ continue below prevents Fill
	  if (!acc_gen) { continue; }
	  
	  // If already matched to a GEN muon, ignore
	  if ( std::find(matched_muon.begin(), matched_muon.end(), igen) != matched_muon.end() ) { continue; }
	  float dr = reco::deltaR( gen->eta(), gen->phi(), muon->eta(), muon->phi() );
	  int dq = gen->charge() * muon->charge();
	  if ( dr < temp_dr && (dq==1 || !correctMuonCharge_) ) {
	    //@@ matched!
	    temp_dr = dr;
	    temp_idx = (int)igen;
	  }
	}
	
	// If matched, continue
	if ( temp_idx >= 0 ) {
	  
	  // Store matched result
	  gen_muon_idx = temp_idx;
	  gen_muon_match = (gen_muon_idx>=0);
	  gen_muon_dr = temp_dr;

	  // Record matched result
	  matched_muon.push_back(gen_muon_idx);
	  
	  // Store 
	  reco::GenParticlePtr gen(genParticles,gen_muon_idx);
	  gen_idx = (int)gen_muon_idx;
	  gen_pt = gen->pt();
	  gen_eta = gen->eta();
	  gen_phi = gen->phi();
	  gen_charge = gen->charge();
	  //gen_dxy = gen->dxy();
	  gen_pdgid = gen->pdgId();
	  //gen_mom_pdgid = gen->pt();
	  //gen_gran_pdgid = gen->pt();
	  
	} // GEN matched to reco?
      } // TRG matched to reco?

      outputTree->Fill();

    } // HW quality
  } // Loop over TRG muons
  
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
  outputTree->Branch("gen_muon_idx",&gen_muon_idx,"gen_muon_idx/I");
  outputTree->Branch("gen_muon_match",&gen_muon_match,"gen_muon_match/I");
  outputTree->Branch("gen_muon_dr",&gen_muon_dr,"gen_muon_dr/F");

  // RECO
  outputTree->Branch("muon_pt",&muon_pt);
  outputTree->Branch("muon_eta",&muon_eta);
  outputTree->Branch("muon_phi",&muon_phi);
  outputTree->Branch("muon_charge",&muon_charge);
  outputTree->Branch("muon_etaAtSt1",&muon_etaAtSt1);
  outputTree->Branch("muon_etaAtSt2",&muon_etaAtSt2);
  outputTree->Branch("muon_phiAtSt1",&muon_phiAtSt1);
  outputTree->Branch("muon_phiAtSt2",&muon_phiAtSt2);
  outputTree->Branch("muon_dz",&muon_dz);
  outputTree->Branch("muon_dzError",&muon_dzError);
  outputTree->Branch("muon_dxy",&muon_dxy);
  outputTree->Branch("muon_dxyError",&muon_dxyError);
  outputTree->Branch("muon_3dIP",&muon_3dIP);
  outputTree->Branch("muon_3dIPError",&muon_3dIPError);
  outputTree->Branch("muon_PassTightID",&muon_PassTightID);
  outputTree->Branch("muon_PassLooseID",&muon_PassLooseID);
  outputTree->Branch("muon_isSAMuon",&muon_isSAMuon);
  outputTree->Branch("muon_trg_idx",&muon_trg_idx,"muon_trg_idx/I");
  outputTree->Branch("muon_trg_match",&muon_trg_match,"muon_trg_match/I");
  outputTree->Branch("muon_trg_dr",&muon_trg_dr,"muon_trg_dr/F");

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
  gen_muon_idx = -1;
  gen_muon_match = 0;
  gen_muon_dr = 10.;
  
  // RECO
  muon_pt = 0.;
  muon_eta = 0.;
  muon_phi = 0.;
  muon_charge = 0;
  muon_etaAtSt1 = 0.;
  muon_etaAtSt2 = 0.;
  muon_phiAtSt1 = 0.;
  muon_phiAtSt2 = 0.;
  muon_dz = 0.;
  muon_dzError = 0.;
  muon_dxy = 0.;
  muon_dxyError = 0.;
  muon_3dIP = 0.;
  muon_3dIPError = 0.;
  muon_PassTightID = 0;
  muon_PassLooseID = 0;
  muon_isSAMuon = 0;
  muon_trg_idx = -1;
  muon_trg_match = 0;
  muon_trg_dr = 10.;
  
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
