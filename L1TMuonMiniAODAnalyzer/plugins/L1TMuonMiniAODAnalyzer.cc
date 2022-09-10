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

#include "L1TMuonMiniAODAnalyzer/L1TMuonMiniAODAnalyzer/plugin/L1TMuonMiniAODAnalyzer.h"


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
    genpartToken_(consumes<GenParticleCollection> (iConfig.getParameter<edm::InputTag>("GenParticles"))),
    packedgenpartToken_(consumes<std::vector< pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag>("PackedGenParticles"))),
    trgresultsToken_(consumes<TriggerResults>(iConfig.getParameter<edm::InputTag>("Triggers"))),
    trigobjectToken_(consumes<pat::TriggerObjectStandAloneCollection>(edm::InputTag("slimmedPatTrigger"))),
    l1MuonToken_(consumes<l1t::MuonBxCollection>(edm::InputTag("gmtStage2Digis","Muon"))),

    MuonPtCut_(iConfig.getParameter<double>("MuonPtCut")),
    SaveTree_(iConfig.getParameter<bool>("SaveTree")),
    IsMC_(iConfig.getParameter<bool>("IsMC")),
    Debug_(iConfig.getParameter<bool>("Debug")),

    muPropagatorSetup1st_(iConfig.getParameter<edm::ParameterSet>("muProp1st"), consumesCollector()),
    muPropagatorSetup2nd_(iConfig.getParameter<edm::ParameterSet>("muProp2nd"), consumesCollector())

{
  //now do what ever initialization is needed
  usesResource("TFileService"); // shared resources

  edm::Service<TFileService> fs;
  outputTree = fs->make<TTree>("tree","tree");

}

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
  edm::Handle<TriggerResults> trigResults;
  iEvent.getByToken(trgresultsToken_, trigResults);
  if( !trigResults.failedToGet() ) {
    int N_Triggers = trigResults->size();
    const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);
    for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
      if (trigResults.product()->accept(i_Trig)) {
	      TString TrigPath =trigName.triggerName(i_Trig);
	      if(TrigPath.Contains("HLT_IsoMu27_v"))HLT_IsoMu27 =true;
        if(TrigPath.Contains("HLT_IsoMu24_v"))HLT_IsoMu24 =true;
        if(TrigPath.Contains("HLT_IsoTkMu24_v"))HLT_IsoTkMu24 =true;
	      if(TrigPath.Contains("HLT_Mu50_v"))HLT_Mu50 =true;
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
      l1mu_tfIdx.push_back(it->tfMuonIndex());

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
    muon_ptcorr.push_back( ptmuoncorr );
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
    TrajectoryStateOnSurface stateAtMuSt1 = muPropagator1st_.extrapolate(muon);
    if (stateAtMuSt1.isValid()) {
        muon_etaAtSt1.push_back(stateAtMuSt1.globalPosition().eta());
        muon_phiAtSt1.push_back(stateAtMuSt1.globalPosition().phi());
    } else {
        muon_etaAtSt1.push_back(-999);
        muon_phiAtSt1.push_back(-999);
    }

    TrajectoryStateOnSurface stateAtMuSt2 = muPropagator2nd_.extrapolate(muon);
    if (stateAtMuSt2.isValid()) {
        muon_etaAtSt2.push_back(stateAtMuSt2.globalPosition().eta());
        muon_phiAtSt2.push_back(stateAtMuSt2.globalPosition().phi());
    } else {
        muon_etaAtSt2.push_back(-999);
        muon_phiAtSt2.push_back(-999);
    }



    hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered.push_back(PassTriggerLeg("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered","hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered",&*muon,iEvent));


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
  outputTree->Branch("muon_size", &muon_size, "muon_size/I");


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
  outputTree->Branch("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered",&hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered);

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


void L1TMuonMiniAODAnalyzer::InitandClearStuff() {

  _mll = 0;
  _ptll = 0;
  _pzll = 0;
  _yll = 0;
  _dphill = 0;
  _phill = 0;
  _costhCSll = 0;
  _nElesll = 0 ;

  _n_PV = 0;
  trueNVtx = 0;


  HLT_IsoMu27 = false;
  HLT_IsoMu24 = false;
  HLT_IsoTkMu24 = false;
  HLT_Mu50 = false;

  Flag_goodVertices = false;
  Flag_globalTightHalo2016Filter = false;
  Flag_globalSuperTightHalo2016Filter = false;
  Flag_BadPFMuonFilter = false;
  Flag_BadPFMuonDzFilter = false;

  hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered.clear();

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
  muon_size = 0;

  l1mu_qual.clear();
  l1mu_charge.clear();
  l1mu_pt.clear();
  l1mu_pt_dxy.clear();
  l1mu_dxy.clear();
  l11mu_eta.clear();
  l11mu_etaAtVtx.clear();
  l1mu_phi.clear();
  l1mu_phiAtVtx.clear();
  l1mu_tfIdx.clear();
  l1mu_bx.clear();
  l1mu_size = 0;

}

bool L1TMuonMiniAODAnalyzer::PassTriggerLeg(std::string triggerlegstring, std::string triggerlegstringalt,  const pat::Muon *muonit, const edm::Event& iEvent ){

  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::InputTag triggerBits_("TriggerResults","","HLT");
  edm::InputTag  triggerObjects_("slimmedPatTrigger");
  iEvent.getByToken(trigobjectToken_, triggerObjects);
  iEvent.getByToken(trgresultsToken_, triggerBits);
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    obj.unpackNamesAndLabels(iEvent,*triggerBits);
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h){
      string myfillabl=obj.filterLabels()[h];
      if( (myfillabl.find(triggerlegstring)!=string::npos  || myfillabl.find(triggerlegstringalt)!=string::npos)&& deltaR(muonit->eta(),muonit->phi(), obj.eta(),obj.phi())<0.4 ) return true;
      if( (myfillabl.find(triggerlegstring)!=string::npos  || myfillabl.find(triggerlegstringalt)!=string::npos)&& myfillabl.find("hltL1s")!=string::npos &&deltaR(muonit->eta(),muonit->phi(), obj.eta(),obj.phi())<0.5 ) return true;
    }
  }
  return false;
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
