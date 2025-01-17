import FWCore.ParameterSet.Config as cms

L1TMuonMiniAODAnalyzerFlat = cms.EDAnalyzer(
    'L1TMuonMiniAODAnalyzerFlat',
    
    Vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    Muons = cms.InputTag("slimmedMuons"),
    GenParticles = cms.InputTag("prunedGenParticles"),
    L1T = cms.InputTag("gmtStage2Digis","Muon"),

    IsMC = cms.bool(True),
    #
    genMuonPtThreshold = cms.double(3.),
    recoMuonPtThreshold = cms.double(0.),
    trgMuonPtThreshold = cms.double(0.),
    #
    genMuonEtaThreshold = cms.double(0.8),
    recoMuonEtaThreshold = cms.double(0.8),
    trgMuonEtaThreshold = cms.double(0.8),
    #
    genRecoMatchingMaxDR = cms.double(0.01),
    trgMuonMatchingMaxDR = cms.double(0.4),
    #
    correctMuonCharge = cms.bool(True),
    correctTrgCharge = cms.bool(True),
    #
    minTrgQuality = cms.int32(12),
    calcPurity = cms.bool(False),
    
    # muon track extrapolation to 1st station
    muProp1st = cms.PSet(
          useTrack = cms.string("tracker"),  # 'none' to use Candidate P4; or 'tracker', 'muon', 'global'
          useState = cms.string("atVertex"), # 'innermost' and 'outermost' require the TrackExtra
          useSimpleGeometry = cms.bool(True),
    useStation2 = cms.bool(False),
          fallbackToME1 = cms.bool(False),
          cosmicPropagationHypothesis = cms.bool(False),
          useMB2InOverlap = cms.bool(False),
          propagatorAlong = cms.ESInputTag("", "SteppingHelixPropagatorAlong"),
          propagatorAny = cms.ESInputTag("", "SteppingHelixPropagatorAny"),
          propagatorOpposite = cms.ESInputTag("", "SteppingHelixPropagatorOpposite")
    ),

    # muon track extrapolation to 2nd station
    muProp2nd = cms.PSet(
          useTrack = cms.string("tracker"),  # 'none' to use Candidate P4; or 'tracker', 'muon', 'global'
          useState = cms.string("atVertex"), # 'innermost' and 'outermost' require the TrackExtra
          useSimpleGeometry = cms.bool(True),
          useStation2 = cms.bool(True),
          fallbackToME1 = cms.bool(False),
          cosmicPropagationHypothesis = cms.bool(False),
          useMB2InOverlap = cms.bool(False),
          propagatorAlong = cms.ESInputTag("", "SteppingHelixPropagatorAlong"),
          propagatorAny = cms.ESInputTag("", "SteppingHelixPropagatorAny"),
          propagatorOpposite = cms.ESInputTag("", "SteppingHelixPropagatorOpposite")
    ),

)
