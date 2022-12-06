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

#include "L1TMuonMiniAODAnalyzer.h"


//
// constructors and destructor
//
L1TMuonMiniAODAnalyzer::L1TMuonMiniAODAnalyzer(const edm::ParameterSet& iConfig)
    :
    metfilterspatToken_(consumes<TriggerResults>(iConfig.getParameter<edm::InputTag>("METFiltersPAT"))),
    metfiltersrecoToken_(consumes<TriggerResults>(iConfig.getParameter<edm::InputTag>("METFiltersRECO"))),
    verticesToken_(consumes<std::vector<Vertex> > (iConfig.getParameter<edm::InputTag>("Vertices"))),
    pfcandsToken_(consumes<std::vector< pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCandidates"))),
    muonToken_(consumes< std::vector< pat::Muon> >(iConfig.getParameter<edm::InputTag>("Muons"))),
    dispMuonToken_(consumes< std::vector< pat::Muon> >(iConfig.getParameter<edm::InputTag>("DispMuons"))),
    genpartToken_(consumes<GenParticleCollection> (iConfig.getParameter<edm::InputTag>("GenParticles"))),
    packedgenpartToken_(consumes<std::vector< pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag>("PackedGenParticles"))),
    trgresultsToken_(consumes<TriggerResults>(iConfig.getParameter<edm::InputTag>("Triggers"))),
    trigobjectToken_(consumes<pat::TriggerObjectStandAloneCollection>(edm::InputTag("slimmedPatTrigger"))),
    l1MuonToken_(consumes<l1t::MuonBxCollection>(edm::InputTag("gmtStage2Digis","Muon"))),
    ugtToken_(consumes<GlobalAlgBlkBxCollection>(edm::InputTag("gtStage2Digis"))),
    l1GtMenuToken_(esConsumes<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd>()),


    MuonPtCut_(iConfig.getParameter<double>("MuonPtCut")),
    SaveTree_(iConfig.getParameter<bool>("SaveTree")),
    IsMC_(iConfig.getParameter<bool>("IsMC")),
    Debug_(iConfig.getParameter<bool>("Debug")),

    muPropagatorSetup1st_(iConfig.getParameter<edm::ParameterSet>("muProp1st"), consumesCollector()),
    muPropagatorSetup2nd_(iConfig.getParameter<edm::ParameterSet>("muProp2nd"), consumesCollector()),
    results_(nullptr)

{
  //now do what ever initialization is needed
  // usesResource("TFileService"); // shared resources

  edm::Service<TFileService> fs;
  outputTree = fs->make<TTree>("tree","tree");

}

L1TMuonMiniAODAnalyzer::~L1TMuonMiniAODAnalyzer() {}

//
// member functions
//

// ------------ method called for each event  ------------
void L1TMuonMiniAODAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  InitandClearStuff();

  muPropagator1st_ = muPropagatorSetup1st_.init(iSetup);
  muPropagator2nd_ = muPropagatorSetup2nd_.init(iSetup);


  _runNb = iEvent.id().run();
  _eventNb = iEvent.id().event();
  _lumiBlock = iEvent.luminosityBlock();
  _bx=iEvent.bunchCrossing();

  //Triggers
  unsigned long long id = iSetup.get<L1TUtmTriggerMenuRcd>().cacheIdentifier();
  if (id != cache_id_) {
    cache_id_ = id;
    edm::ESHandle<L1TUtmTriggerMenu> menu;
    menu = iSetup.getHandle(l1GtMenuToken_);

    for (auto const &keyval : menu->getAlgorithmMap()) {
      std::string const &name = keyval.second.getName();
      unsigned int index = keyval.second.getIndex();
      //std::cerr << fmt::sprintf("bit %4d: %s", index, name) << std::endl;
      outputTree->SetAlias(name.c_str(), fmt::sprintf("L1uGT.m_algoDecisionInitial[%d]", index).c_str());
    }
  }

  edm::Handle<GlobalAlgBlkBxCollection> ugt;
  iEvent.getByToken(ugtToken_, ugt);
  if (ugt.isValid()) {
    results_ = &ugt->at(0, 0);
  }

  edm::Handle<TriggerResults> trigResults;
  iEvent.getByToken(trgresultsToken_, trigResults);
  if( !trigResults.failedToGet() ) {
    int N_Triggers = trigResults->size();
    edm::TriggerNames const& trigName = iEvent.triggerNames(*trigResults);
    for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
      if (trigResults.product()->accept(i_Trig)) {
	      TString TrigPath = trigName.triggerName(i_Trig);
	      if(TrigPath.Contains("HLT_IsoMu27_v"))HLT_IsoMu27 =true;
        if(TrigPath.Contains("HLT_IsoMu24_v"))HLT_IsoMu24 =true;
        if(TrigPath.Contains("HLT_IsoTkMu24_v"))HLT_IsoTkMu24 =true;
	      if(TrigPath.Contains("HLT_Mu50_v"))HLT_Mu50 =true;
	      if(TrigPath.Contains("HLT_DoubleMu43NoFiltersNoVtx_v"))HLT_DoubleMu43NoFiltersNoVtx =true;
	      if(TrigPath.Contains("HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm_v"))HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm =true;
	      if(TrigPath.Contains("HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm_v"))HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm =true;
	      if(TrigPath.Contains("HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed_v"))HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed =true;
	      if(TrigPath.Contains("HLT_DoubleL2Mu23NoVtx_2Cha_v"))HLT_DoubleL2Mu23NoVtx_2Cha =true;
	      if(TrigPath.Contains("HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm_v"))HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm =true;
	      if(TrigPath.Contains("HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm_v"))HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm =true;
	      if(TrigPath.Contains("HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm_v"))HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm =true;
	      if(TrigPath.Contains("HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX_v"))HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX =true;
	      if(TrigPath.Contains("HLT_DoubleMu4_3_Bs_v"))HLT_DoubleMu4_3_Bs =true;
	      if(TrigPath.Contains("HLT_DoubleMu4_3_Jpsi_v"))HLT_DoubleMu4_3_Jpsi =true;
	      if(TrigPath.Contains("HLT_DoubleMu4_3_Photon4_BsToMMG_v"))HLT_DoubleMu4_3_Photon4_BsToMMG =true;
	      if(TrigPath.Contains("HLT_DoubleMu4_LowMass_Displaced_v"))HLT_DoubleMu4_LowMass_Displaced =true;
	      if(TrigPath.Contains("HLT_DoubleMu4_3_LowMass_v"))HLT_DoubleMu4_3_LowMass =true;
	      if(TrigPath.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu_v"))HLT_DoubleMu3_TkMu_DsTau3Mu =true;
	      if(TrigPath.Contains("HLT_DoubleMu3_Trk_Tau3mu_v"))HLT_DoubleMu3_Trk_Tau3mu =true;
      }
    }
  }

  edm::Handle<l1t::MuonBxCollection> l1muoncoll;
  iEvent.getByToken(l1MuonToken_ , l1muoncoll);
  //  const  int nbx = IsMC_ ? 0:2;
  for(int i = l1muoncoll->getFirstBX() ; i<= l1muoncoll->getLastBX() ;i++){
    for( l1t::MuonBxCollection::const_iterator l1muonit= l1muoncoll->begin(i); l1muonit != l1muoncoll->end(i) ; ++l1muonit){
      if(l1muonit->pt() < 0) continue;
      l1mu_qual.push_back( l1muonit->hwQual() );
      l1mu_charge.push_back( l1muonit->charge() );
      l1mu_pt.push_back( l1muonit->pt() );
      l1mu_pt_dxy.push_back( l1muonit->ptUnconstrained() );
      l1mu_dxy.push_back( l1muonit->hwDXY() );
      l1mu_eta.push_back( l1muonit->eta() );
      l1mu_etaAtVtx.push_back( l1muonit->etaAtVtx() );
      l1mu_phi.push_back( l1muonit->phi() );
      l1mu_phiAtVtx.push_back( l1muonit->phiAtVtx() );
      l1mu_tfIdx.push_back(l1muonit->tfMuonIndex());

      l1mu_bx.push_back( i);
      l1mu_size++;

    }
  }

  //Vertices
  // NB there are usually two available collection of vertices in MINIAOD
  edm::Handle<std::vector<Vertex> > theVertices;
  iEvent.getByToken(verticesToken_,theVertices) ;
  _n_PV = theVertices->size();
  Vertex::Point PV(0,0,0);
  if(_n_PV){ PV = theVertices->begin()->position();}


  //MET filters are stored in TriggerResults::RECO or TriggerResults::PAT . Should take the latter if it exists
  edm::Handle<TriggerResults> METFilterResults;
  iEvent.getByToken(metfilterspatToken_, METFilterResults);
  if(!(METFilterResults.isValid())) iEvent.getByToken(metfiltersrecoToken_, METFilterResults);
  
  Flag_goodVertices= GetMETFilterDecision(iEvent,METFilterResults,"Flag_goodVertices");
  Flag_globalTightHalo2016Filter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_globalTightHalo2016Filter");
  Flag_globalSuperTightHalo2016Filter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_globalSuperTightHalo2016Filter");
  Flag_BadPFMuonFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_BadPFMuonFilter");
  Flag_BadPFMuonDzFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_BadPFMuonDzFilter");

  //Muons
  edm::Handle< std::vector<pat::Muon> > thePatMuons;
  iEvent.getByToken(muonToken_,thePatMuons);
  for( std::vector<pat::Muon>::const_iterator muon = (*thePatMuons).begin(); muon != (*thePatMuons).end(); muon++ ) {
    if((&*muon)->pt() <0) continue; //Loose cut  on uncorrected pt 

    //Rochester corrections: https://twiki.cern.ch/twiki/bin/viewauth/CMS/RochcorMuon#Rochester_Correction
    //https://indico.cern.ch/event/926898/contributions/3897122/attachments/2052816/3441285/roccor.pdf
    double ptmuoncorr= (&*muon)->pt();

    // if( !IsMC_) ptmuoncorr *= rc.kScaleDT( (&*muon)->charge(),  (&*muon)->pt(), (&*muon)->eta(),(&*muon)->phi());
    // else{
    //   if( (&*muon)->genLepton() !=0)  ptmuoncorr *= rc.kSpreadMC( (&*muon)->charge(),  (&*muon)->pt(), (&*muon)->eta(),(&*muon)->phi(), (&*muon)->genLepton()->pt() );
    //   else if(! ((&*muon)->innerTrack()).isNull()) ptmuoncorr *= rc.kSmearMC( (&*muon)->charge(),  (&*muon)->pt(), (&*muon)->eta(),(&*muon)->phi(),  (&*muon)->innerTrack()->hitPattern().trackerLayersWithMeasurement(), gRandom->Rndm());
    // }

    // bool passvetoid =  false;
    // passvetoid =  (&*muon)->isStandAloneMuon() || (&*muon)->pt()>3;
    // if(!passvetoid) continue;
    //Counting the number of muons, not all of them will be stored
    // if((&*muon)->pt()<MuonPtCut_)continue;

    // store all reco muons for now
    muon_size++;
    muon_eta.push_back((&*muon)->eta());
    muon_phi.push_back((&*muon)->phi());
    muon_pt.push_back((&*muon)->pt());
    muon_ptCorr.push_back( ptmuoncorr );
    muon_charge.push_back((&*muon)->charge());
    muon_PassTightID.push_back(  (&*muon)->passed(reco::Muon::CutBasedIdMediumPrompt )&& (&*muon)->passed(reco::Muon::PFIsoTight ) );
    muon_PassLooseID.push_back(  (&*muon)->passed(reco::Muon::CutBasedIdLoose )&& (&*muon)->passed(reco::Muon::PFIsoLoose ) );
    muon_isSAMuon.push_back( (&*muon)->isStandAloneMuon());
    if( !((&*muon)->innerTrack()).isNull()){
      muon_dz.push_back( (&*muon)->innerTrack()->dz(PV));
      muon_dzError.push_back( (&*muon)->innerTrack()->dzError());
      muon_dxy.push_back( (&*muon)->innerTrack()->dxy(PV));
      muon_dxyError.push_back( (&*muon)->innerTrack()->dxyError());
      muon_3dIP.push_back( (&*muon)->dB(pat::Muon::PV3D));
      muon_3dIPError.push_back((&*muon)->edB(pat::Muon::PV3D));
    }
    else{
      muon_dz.push_back( -999.);
      muon_dzError.push_back(-999.);
      muon_dxy.push_back(-999.);
      muon_dxyError.push_back(-999.);
      muon_3dIP.push_back(-999.);
      muon_3dIPError.push_back(-999.);
    }


    // extrapolation of muon track coordinates
    TrajectoryStateOnSurface stateAtMuSt1 = muPropagator1st_.extrapolate(*muon);
    if (stateAtMuSt1.isValid()) {
        muon_etaAtSt1.push_back(stateAtMuSt1.globalPosition().eta());
        muon_phiAtSt1.push_back(stateAtMuSt1.globalPosition().phi());
    } else {
        muon_etaAtSt1.push_back(-999);
        muon_phiAtSt1.push_back(-999);
    }

    TrajectoryStateOnSurface stateAtMuSt2 = muPropagator2nd_.extrapolate(*muon);
    if (stateAtMuSt2.isValid()) {
        muon_etaAtSt2.push_back(stateAtMuSt2.globalPosition().eta());
        muon_phiAtSt2.push_back(stateAtMuSt2.globalPosition().phi());
    } else {
        muon_etaAtSt2.push_back(-999);
        muon_phiAtSt2.push_back(-999);
    }



    muon_isIsoHLTMuon.push_back(PassTriggerLeg("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered",&*muon,iEvent));
    muon_isHLTMuon.push_back(PassTriggerLeg("hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q",&*muon,iEvent));
    muon_isTauTo3MuMuon.push_back(PassTriggerLeg("hltL1sDoubleMu0er2p0SQOSdEtaMax1p6orTripleMu21p50",&*muon,iEvent));
    if (PassTriggerLeg("hltL1sDoubleMuForBs",&*muon,iEvent) or
        PassTriggerLeg("hltL1sDoubleMuForBsToMMG",&*muon,iEvent)) {
      muon_isDoubleMuForBsMuon.push_back(true);
    } else {
      muon_isDoubleMuForBsMuon.push_back(false);
    }
    if (PassTriggerLeg("hltL1sDoubleMuForLowMassInclusive",&*muon,iEvent) or
        PassTriggerLeg("hltL1sDoubleMuForLowMassDisplaced",&*muon,iEvent)) {
      muon_isDoubleMuLowMassMuon.push_back(true);
    } else {
      muon_isDoubleMuLowMassMuon.push_back(false);
    }

    muon_isSingleMuMuon.push_back(PassTriggerLeg("hltL1sSingleMu22or25",&*muon,iEvent));

  }

  // Displaced Muons
  edm::Handle< std::vector<pat::Muon> > thePatDispMuons;
  iEvent.getByToken(dispMuonToken_,thePatDispMuons);
  for( std::vector<pat::Muon>::const_iterator muon = (*thePatDispMuons).begin(); muon != (*thePatDispMuons).end(); muon++ ) {
    if((&*muon)->pt() <0) continue; //Loose cut  on uncorrected pt 

    //Rochester corrections: https://twiki.cern.ch/twiki/bin/viewauth/CMS/RochcorMuon#Rochester_Correction
    //https://indico.cern.ch/event/926898/contributions/3897122/attachments/2052816/3441285/roccor.pdf
    double ptmuoncorr= (&*muon)->pt();

    // if( !IsMC_) ptmuoncorr *= rc.kScaleDT( (&*muon)->charge(),  (&*muon)->pt(), (&*muon)->eta(),(&*muon)->phi());
    // else{
    //   if( (&*muon)->genLepton() !=0)  ptmuoncorr *= rc.kSpreadMC( (&*muon)->charge(),  (&*muon)->pt(), (&*muon)->eta(),(&*muon)->phi(), (&*muon)->genLepton()->pt() );
    //   else if(! ((&*muon)->innerTrack()).isNull()) ptmuoncorr *= rc.kSmearMC( (&*muon)->charge(),  (&*muon)->pt(), (&*muon)->eta(),(&*muon)->phi(),  (&*muon)->innerTrack()->hitPattern().trackerLayersWithMeasurement(), gRandom->Rndm());
    // }

    // bool passvetoid =  false;
    // passvetoid =  (&*muon)->isStandAloneMuon() || (&*muon)->pt()>3;
    // if(!passvetoid) continue;
    //Counting the number of muons, not all of them will be stored
    // if((&*muon)->pt()<MuonPtCut_)continue;

    // store all reco muons for now
    dispMuon_size++;
    dispMuon_eta.push_back((&*muon)->eta());
    dispMuon_phi.push_back((&*muon)->phi());
    dispMuon_pt.push_back((&*muon)->pt());
    dispMuon_ptCorr.push_back( ptmuoncorr );
    dispMuon_charge.push_back((&*muon)->charge());
    dispMuon_PassTightID.push_back(  (&*muon)->passed(reco::Muon::CutBasedIdMediumPrompt )&& (&*muon)->passed(reco::Muon::PFIsoTight ) );
    dispMuon_PassLooseID.push_back(  (&*muon)->passed(reco::Muon::CutBasedIdLoose )&& (&*muon)->passed(reco::Muon::PFIsoLoose ) );
    dispMuon_isSAMuon.push_back( (&*muon)->isStandAloneMuon());
    if( !((&*muon)->innerTrack()).isNull()){
      dispMuon_dz.push_back( (&*muon)->innerTrack()->dz(PV));
      dispMuon_dzError.push_back( (&*muon)->innerTrack()->dzError());
      dispMuon_dxy.push_back( (&*muon)->innerTrack()->dxy(PV));
      dispMuon_dxyError.push_back( (&*muon)->innerTrack()->dxyError());
      dispMuon_3dIP.push_back( (&*muon)->dB(pat::Muon::PV3D));
      dispMuon_3dIPError.push_back((&*muon)->edB(pat::Muon::PV3D));
    }
    else{
      dispMuon_dz.push_back( -999.);
      dispMuon_dzError.push_back(-999.);
      dispMuon_dxy.push_back(-999.);
      dispMuon_dxyError.push_back(-999.);
      dispMuon_3dIP.push_back(-999.);
      dispMuon_3dIPError.push_back(-999.);
    }


    // extrapolation of muon track coordinates
    TrajectoryStateOnSurface stateAtMuSt1 = muPropagator1st_.extrapolate(*muon);
    if (stateAtMuSt1.isValid()) {
        dispMuon_etaAtSt1.push_back(stateAtMuSt1.globalPosition().eta());
        dispMuon_phiAtSt1.push_back(stateAtMuSt1.globalPosition().phi());
    } else {
        dispMuon_etaAtSt1.push_back(-999);
        dispMuon_phiAtSt1.push_back(-999);
    }

    TrajectoryStateOnSurface stateAtMuSt2 = muPropagator2nd_.extrapolate(*muon);
    if (stateAtMuSt2.isValid()) {
        dispMuon_etaAtSt2.push_back(stateAtMuSt2.globalPosition().eta());
        dispMuon_phiAtSt2.push_back(stateAtMuSt2.globalPosition().phi());
    } else {
        dispMuon_etaAtSt2.push_back(-999);
        dispMuon_phiAtSt2.push_back(-999);
    }


    if (PassTriggerLeg("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered",&*muon,iEvent) or
        PassTriggerLeg("hltL3fDimuonL1f0CosmicL2NV2Chaf10L3NVf0Veto1PromptSimple",&*muon,iEvent) or
        PassTriggerLeg("hltL2fL1sMuORL1f0DoubleL2NoVtx23Q2ChaCosmicSeed",&*muon,iEvent) or
        PassTriggerLeg("hltL2DoubleMu23NoVertexL2Filtered2Cha",&*muon,iEvent) or
        PassTriggerLeg("hltL3fL1DoubleMuf0L2NVf15f7L3SingleMuNVf0VetoDxy0p1cm",&*muon,iEvent) or
        PassTriggerLeg("hltL3fL1DoubleMuf0L2NVf15f7L3SingleMuNVf16DisplacedHybDxy0p1cm",&*muon,iEvent)) {
      dispMuon_isDispL2HLTDiMuon.push_back(true);
    } else {
      dispMuon_isDispL2HLTDiMuon.push_back(false);
    }

    if (PassTriggerLeg("hltL3fL1DoubleMuf0L2NVf15f7L3DoubleMuNVf10Displaced",&*muon,iEvent) or
        PassTriggerLeg("hltL3fL1DoubleMuf0L2NVf15f7L3SingleMuNVf16Displaced",&*muon,iEvent) or
        PassTriggerLeg("hltL3dTkfL1DoubleMuf0L2NVf15f7L3SingleMuNVf16Displaced",&*muon,iEvent) or
        PassTriggerLeg("hltL3dTkfL1DoubleMuf0L2NVf15f7L3DoubleMuNVf10Displaced",&*muon,iEvent) or
        PassTriggerLeg("hltL3fDimuonL1f0L2NVf16L3NoFiltersNoVtxFiltered43",&*muon,iEvent)) {
      dispMuon_isDispHLTDiMuon.push_back(true);
    } else {
      dispMuon_isDispHLTDiMuon.push_back(false);
    }

    if (PassTriggerLeg("hltL2fL1fMuf0DoubleL2MuNoVtx10Q",&*muon,iEvent) or
        PassTriggerLeg("hltL2fL1sMuORL1f0DoubleL2NoVtx23Q",&*muon,iEvent) or
        PassTriggerLeg("hltL2fL1DoubleMuf0L2DoubleMuNoVtx2ChaFiltered7",&*muon,iEvent) or
        PassTriggerLeg("hltL2fL1DoubleMuf0L2DoubleMuNoVtxFiltered7",&*muon,iEvent) or
        PassTriggerLeg("hltL2fL1DoubleMuf0L2DoubleMuDisplacedFiltered7",&*muon,iEvent) or
        PassTriggerLeg("hltL2fDimuonL1f0L2NoVtxFiltered16",&*muon,iEvent)) {
      dispMuon_isL2DispHLTMuon.push_back(true);
    } else {
      dispMuon_isL2DispHLTMuon.push_back(false);
    }

    dispMuon_isNoBPTXHLTMuon.push_back(PassTriggerLeg("hltL2fL1sMuOpenNotBptxORNoHaloMu3BXL1f0NoVtxCosmicSeedMeanTimerL2Filtered40Sta3",&*muon,iEvent));
    dispMuon_isL2NoVtxHLTMuon.push_back(PassTriggerLeg("hltL2fL1fMuf0DoubleL2MuNoVtx10Q", &*muon, iEvent));
    dispMuon_isL2CosmicSeedHLTMuon.push_back(PassTriggerLeg("hltL2CosmicMuons", &*muon, iEvent));

    dispMuon_isL1NoBPTXMuon.push_back(PassTriggerLeg("hltL1sSingleMuOpenEr1p4NotBptxOR3BXORL1sSingleMuOpenEr1p1NotBptxOR3BX",&*muon,iEvent));
    dispMuon_isL1DispMuon.push_back(PassTriggerLeg("hltL1sDoubleMu155SQOR157ORTripleMu444ORDoubleMu0upt",&*muon,iEvent));


  }


  if(SaveTree_)outputTree->Fill();

}

// ------------ method called once each job just before starting event loop  ------------
void L1TMuonMiniAODAnalyzer::beginJob() {
  // please remove this method if not needed

  outputTree->Branch("_eventNb",   &_eventNb,   "_eventNb/l");
  outputTree->Branch("_runNb",     &_runNb,     "_runNb/l");
  outputTree->Branch("_lumiBlock", &_lumiBlock, "_lumiBlock/l");
  outputTree->Branch("_bx", &_bx, "_bx/l");
  outputTree->Branch("_n_PV", &_n_PV, "_n_PV/I");

  if(IsMC_){
    outputTree->Branch("trueNVtx", &trueNVtx,"trueNVtx/I");
  }

  outputTree->Branch("L1uGT", "GlobalAlgBlk", &results_, 32000, 3);

  outputTree->Branch("Flag_goodVertices",&Flag_goodVertices,"Flag_goodVertices/O");
  outputTree->Branch("Flag_globalTightHalo2016Filter",&Flag_globalTightHalo2016Filter,"Flag_globalTightHalo2016Filter/O");
  outputTree->Branch("Flag_globalSuperTightHalo2016Filter",&Flag_globalSuperTightHalo2016Filter,"Flag_globalSuperTightHalo2016Filter/O");
  outputTree->Branch("Flag_BadPFMuonFilter",&Flag_BadPFMuonFilter,"Flag_BadPFMuonFilter/O");
  outputTree->Branch("Flag_BadPFMuonDzFilter",&Flag_BadPFMuonDzFilter,"Flag_BadPFMuonDzFilter/O");

  outputTree->Branch("muon_eta",&muon_eta);
  outputTree->Branch("muon_etaAtSt1",&muon_etaAtSt1);
  outputTree->Branch("muon_etaAtSt2",&muon_etaAtSt2);
  outputTree->Branch("muon_phi",&muon_phi);
  outputTree->Branch("muon_phiAtSt1",&muon_phiAtSt1);
  outputTree->Branch("muon_phiAtSt2",&muon_phiAtSt2);
  outputTree->Branch("muon_pt",&muon_pt);
  outputTree->Branch("muon_ptCorr",&muon_ptCorr);
  outputTree->Branch("muon_charge",&muon_charge);

  outputTree->Branch("muon_dz",&muon_dz);
  outputTree->Branch("muon_dzError",&muon_dzError);
  outputTree->Branch("muon_dxy",&muon_dxy);
  outputTree->Branch("muon_dxyError",&muon_dxyError);
  outputTree->Branch("muon_3dIP",&muon_3dIP);
  outputTree->Branch("muon_3dIPError",&muon_3dIPError);

  outputTree->Branch("muon_PassTightID",&muon_PassTightID);
  outputTree->Branch("muon_PassLooseID",&muon_PassLooseID);
  outputTree->Branch("muon_isSAMuon",&muon_isSAMuon);

  outputTree->Branch("muon_isIsoHLTMuon",&muon_isIsoHLTMuon);
  outputTree->Branch("muon_isHLTMuon",&muon_isHLTMuon);
  outputTree->Branch("muon_isTauTo3MuMuon",&muon_isTauTo3MuMuon);
  outputTree->Branch("muon_isDoubleMuForBsMuon",&muon_isDoubleMuForBsMuon);
  outputTree->Branch("muon_isDoubleMuLowMassMuon",&muon_isDoubleMuLowMassMuon);
  outputTree->Branch("muon_isSingleMuMuon",&muon_isSingleMuMuon);

  outputTree->Branch("muon_size", &muon_size, "muon_size/I");

  outputTree->Branch("dispMuon_eta",&dispMuon_eta);
  outputTree->Branch("dispMuon_etaAtSt1",&dispMuon_etaAtSt1);
  outputTree->Branch("dispMuon_etaAtSt2",&dispMuon_etaAtSt2);
  outputTree->Branch("dispMuon_phi",&dispMuon_phi);
  outputTree->Branch("dispMuon_phiAtSt1",&dispMuon_phiAtSt1);
  outputTree->Branch("dispMuon_phiAtSt2",&dispMuon_phiAtSt2);
  outputTree->Branch("dispMuon_pt",&dispMuon_pt);
  outputTree->Branch("dispMuon_ptCorr",&dispMuon_ptCorr);
  outputTree->Branch("dispMuon_charge",&dispMuon_charge);

  outputTree->Branch("dispMuon_dz",&dispMuon_dz);
  outputTree->Branch("dispMuon_dzError",&dispMuon_dzError);
  outputTree->Branch("dispMuon_dxy",&dispMuon_dxy);
  outputTree->Branch("dispMuon_dxyError",&dispMuon_dxyError);
  outputTree->Branch("dispMuon_3dIP",&dispMuon_3dIP);
  outputTree->Branch("dispMuon_3dIPError",&dispMuon_3dIPError);

  outputTree->Branch("dispMuon_PassTightID",&dispMuon_PassTightID);
  outputTree->Branch("dispMuon_PassLooseID",&dispMuon_PassLooseID);
  outputTree->Branch("dispMuon_isSAdispMuon",&dispMuon_isSAMuon);

  outputTree->Branch("dispMuon_isDispL2HLTDiMuon",&dispMuon_isDispL2HLTDiMuon);
  outputTree->Branch("dispMuon_isDispHLTDiMuon",&dispMuon_isDispHLTDiMuon);
  outputTree->Branch("dispMuon_isL2DispHLTMuon",&dispMuon_isL2DispHLTMuon);
  outputTree->Branch("dispMuon_isNoBPTXHLTMuon",&dispMuon_isNoBPTXHLTMuon);
  outputTree->Branch("dispMuon_isL2NoVtxHLTMuon",&dispMuon_isL2NoVtxHLTMuon);
  outputTree->Branch("dispMuon_isL2CosmicSeedHLTMuon",&dispMuon_isL2CosmicSeedHLTMuon);
  outputTree->Branch("dispMuon_isL1NoBPTXMuon",&dispMuon_isL1NoBPTXMuon);
  outputTree->Branch("dispMuon_isL1DispMuon",&dispMuon_isL1DispMuon);

  outputTree->Branch("dispMuon_size", &dispMuon_size, "dispMuon_size/I");


  outputTree->Branch("_mll", &_mll, "_mll/f");
  outputTree->Branch("_ptll", &_ptll, "_ptll/f");
  outputTree->Branch("_pzll", &_pzll, "_pzll/f");
  outputTree->Branch("_yll", &_yll, "_yll/f");
  outputTree->Branch("_dphill", &_dphill, "_dphill/f");
  outputTree->Branch("_phill", &_phill, "_phill/f");
  outputTree->Branch("_costhCSll", &_costhCSll, "_costhCSll/f");

  outputTree->Branch("HLT_IsoMu27",&HLT_IsoMu27,"HLT_IsoMu27/O");
  outputTree->Branch("HLT_IsoMu24",&HLT_IsoMu24,"HLT_IsoMu24/O");
  outputTree->Branch("HLT_IsoTkMu24",&HLT_IsoTkMu24,"HLT_IsoTkMu24/O");
  outputTree->Branch("HLT_Mu50",&HLT_Mu50,"HLT_Mu50/O");
  outputTree->Branch("HLT_DoubleMu43NoFiltersNoVtx",&HLT_DoubleMu43NoFiltersNoVtx,"HLT_MHLT_DoubleMu43NoFiltersNoVtxu50/O");
  outputTree->Branch("HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm",&HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm,"HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm/O");
  outputTree->Branch("HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm",&HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm,"HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm/O");
  outputTree->Branch("HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed",&HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed,"HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed/O");
  outputTree->Branch("HLT_DoubleL2Mu23NoVtx_2Cha",&HLT_DoubleL2Mu23NoVtx_2Cha,"HLT_DoubleL2Mu23NoVtx_2Cha/O");
  outputTree->Branch("HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm",&HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm,"HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm/O");
  outputTree->Branch("HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm",&HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm,"HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm/O");
  outputTree->Branch("HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm",&HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm,"HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm/O");
  outputTree->Branch("HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX",&HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX,"HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX/O");
  outputTree->Branch("HLT_DoubleMu4_3_Bs",&HLT_DoubleMu4_3_Bs,"HLT_DoubleMu4_3_Bs/O");
  outputTree->Branch("HLT_DoubleMu4_3_Jpsi",&HLT_DoubleMu4_3_Jpsi,"HLT_DoubleMu4_3_Jpsi/O");
  outputTree->Branch("HLT_DoubleMu4_3_Photon4_BsToMMG",&HLT_DoubleMu4_3_Photon4_BsToMMG,"HLT_DoubleMu4_3_Photon4_BsToMMG/O");
  outputTree->Branch("HLT_DoubleMu4_LowMass_Displaced",&HLT_DoubleMu4_LowMass_Displaced,"HLT_DoubleMu4_LowMass_Displaced/O");
  outputTree->Branch("HLT_DoubleMu4_3_LowMass",&HLT_DoubleMu4_3_LowMass,"HLT_DoubleMu4_3_LowMass/O");
  outputTree->Branch("HLT_DoubleMu3_TkMu_DsTau3Mu",&HLT_DoubleMu3_TkMu_DsTau3Mu,"HLT_DoubleMu3_TkMu_DsTau3Mu/O");
  outputTree->Branch("HLT_DoubleMu3_Trk_Tau3mu",&HLT_DoubleMu3_Trk_Tau3mu,"HLT_DoubleMu3_Trk_Tau3mu/O");

  outputTree->Branch("l1mu_qual",&l1mu_qual);
  outputTree->Branch("l1mu_charge",&l1mu_charge);
  outputTree->Branch("l1mu_pt",&l1mu_pt);
  outputTree->Branch("l1mu_pt_dxy",&l1mu_pt_dxy);
  outputTree->Branch("l1mu_dxy",&l1mu_dxy);
  outputTree->Branch("l1mu_eta",&l1mu_eta);
  outputTree->Branch("l1mu_etaAtVtx",&l1mu_etaAtVtx);
  outputTree->Branch("l1mu_phi",&l1mu_phi);
  outputTree->Branch("l1mu_phiAtVtx",&l1mu_phiAtVtx);
  outputTree->Branch("l1mu_tfIdx",&l1mu_tfIdx);
  outputTree->Branch("l1mu_bx",&l1mu_bx);
  outputTree->Branch("l1mu_size", &l1mu_size, "l1mu_size/I");


}

void L1TMuonMiniAODAnalyzer::endJob() {}

void L1TMuonMiniAODAnalyzer::InitandClearStuff() {

  _mll = 0;
  _ptll = 0;
  _pzll = 0;
  _yll = 0;
  _dphill = 0;
  _phill = 0;
  _costhCSll = 0;

  _n_PV = 0;
  trueNVtx = 0;

  results_ = nullptr;
  cache_id_ = 0;

  HLT_IsoMu27 = false;
  HLT_IsoMu24 = false;
  HLT_IsoTkMu24 = false;
  HLT_Mu50 = false;
  HLT_DoubleMu43NoFiltersNoVtx = false;
  HLT_DoubleL3Mu16_10NoVtx_DxyMin0p01cm = false;
  HLT_DoubleL3dTksMu16_10NoVtx_DxyMin0p01cm = false;
  HLT_DoubleL2Mu23NoVtx_2Cha_CosmicSeed = false;
  HLT_DoubleL2Mu23NoVtx_2Cha = false;
  HLT_DoubleL2Mu10NoVtx_2Cha_VetoL3Mu0DxyMax1cm = false;
  HLT_DoubleL2Mu_L3Mu16NoVtx_VetoL3Mu0DxyMax0p1cm = false;
  HLT_DoubleL2Mu10NoVtx_2Cha_CosmicSeed_VetoL3Mu0DxyMax1cm = false;
  HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX = false;
  HLT_DoubleMu4_3_Bs = false;
  HLT_DoubleMu4_3_Jpsi = false;
  HLT_DoubleMu4_3_Photon4_BsToMMG = false;
  HLT_DoubleMu4_LowMass_Displaced = false;
  HLT_DoubleMu4_3_LowMass = false;
  HLT_DoubleMu3_TkMu_DsTau3Mu = false;
  HLT_DoubleMu3_Trk_Tau3mu = false;


  Flag_goodVertices = false;
  Flag_globalTightHalo2016Filter = false;
  Flag_globalSuperTightHalo2016Filter = false;
  Flag_BadPFMuonFilter = false;
  Flag_BadPFMuonDzFilter = false;

  muon_eta.clear();
  muon_etaAtSt1.clear();
  muon_etaAtSt2.clear();
  muon_phi.clear();
  muon_phiAtSt1.clear();
  muon_phiAtSt2.clear();
  muon_pt.clear();
  muon_ptCorr.clear();
  muon_charge.clear();

  muon_dz.clear();
  muon_dzError.clear();
  muon_dxy.clear();
  muon_dxyError.clear();
  muon_3dIP.clear();
  muon_3dIPError.clear();

  muon_PassTightID.clear();
  muon_PassLooseID.clear();
  muon_isSAMuon.clear() ;

  muon_isIsoHLTMuon.clear();
  muon_isHLTMuon.clear();
  muon_isTauTo3MuMuon.clear();
  muon_isDoubleMuForBsMuon.clear();
  muon_isDoubleMuLowMassMuon.clear();
  muon_isSingleMuMuon.clear();

  muon_size = 0;

  dispMuon_eta.clear();
  dispMuon_etaAtSt1.clear();
  dispMuon_etaAtSt2.clear();
  dispMuon_phi.clear();
  dispMuon_phiAtSt1.clear();
  dispMuon_phiAtSt2.clear();
  dispMuon_pt.clear();
  dispMuon_ptCorr.clear();
  dispMuon_charge.clear();

  dispMuon_dz.clear();
  dispMuon_dzError.clear();
  dispMuon_dxy.clear();
  dispMuon_dxyError.clear();
  dispMuon_3dIP.clear();
  dispMuon_3dIPError.clear();

  dispMuon_PassTightID.clear();
  dispMuon_PassLooseID.clear();
  dispMuon_isSAMuon.clear() ;

  dispMuon_isDispL2HLTDiMuon.clear();
  dispMuon_isDispHLTDiMuon.clear();
  dispMuon_isL2DispHLTMuon.clear();
  dispMuon_isNoBPTXHLTMuon.clear();
  dispMuon_isL2NoVtxHLTMuon.clear();
  dispMuon_isL2CosmicSeedHLTMuon.clear();
  dispMuon_isL1NoBPTXMuon.clear();
  dispMuon_isL1DispMuon.clear();

  dispMuon_size = 0;

  l1mu_qual.clear();
  l1mu_charge.clear();
  l1mu_pt.clear();
  l1mu_pt_dxy.clear();
  l1mu_dxy.clear();
  l1mu_eta.clear();
  l1mu_etaAtVtx.clear();
  l1mu_phi.clear();
  l1mu_phiAtVtx.clear();
  l1mu_tfIdx.clear();
  l1mu_bx.clear();
  l1mu_size = 0;

}

bool L1TMuonMiniAODAnalyzer::PassTriggerLeg(std::string triggerlegstring, const pat::Muon *muonit, const edm::Event& iEvent ){

  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::InputTag triggerBits_("TriggerResults","","HLT");
  edm::InputTag  triggerObjects_("slimmedPatTrigger");
  iEvent.getByToken(trigobjectToken_, triggerObjects);
  iEvent.getByToken(trgresultsToken_, triggerBits);
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    obj.unpackNamesAndLabels(iEvent,*triggerBits);
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
      TString myfillabl=obj.filterLabels()[h];
      if( (myfillabl.Contains(triggerlegstring)) && deltaR(muonit->eta(),muonit->phi(), obj.eta(),obj.phi())<0.4 ) return true;
      if( (myfillabl.Contains(triggerlegstring)) && myfillabl.Contains("hltL1s") && deltaR(muonit->eta(),muonit->phi(), obj.eta(),obj.phi())<0.5 ) return true;
    }
  }
  return false;
}

bool L1TMuonMiniAODAnalyzer::GetMETFilterDecision(const edm::Event& iEvent,edm::Handle<TriggerResults> METFilterResults, TString studiedfilter){
  
  if( !METFilterResults.failedToGet() ) {
    int N_MetFilters = METFilterResults->size();
    const edm::TriggerNames & metfilterName = iEvent.triggerNames(*METFilterResults);
    for( int i_Metfilter = 0; i_Metfilter < N_MetFilters; ++i_Metfilter ) {
      TString MetfilterPath =metfilterName.triggerName(i_Metfilter);
      //      cout << MetfilterPath<<endl;
      if(MetfilterPath.Index(studiedfilter) >=0)  return METFilterResults.product()->accept(i_Metfilter);

    }
  }
   return true; 
}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1TMuonMiniAODAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TMuonMiniAODAnalyzer);
