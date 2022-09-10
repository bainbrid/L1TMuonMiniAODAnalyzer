import FWCore.ParameterSet.Config as cms

L1TMuonMiniAODAnalyzer = cms.EDAnalyzer('L1TMuonMiniAODAnalyzer',

    METFiltersPAT = cms.InputTag("TriggerResults::PAT"),
    METFiltersRECO = cms.InputTag("TriggerResults::RECO"),
    Vertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
    PFCandidates=cms.InputTag("packedPFCandidates"),
    Muons=cms.InputTag("slimmedMuons"),
    GenParticles=cms.InputTag("prunedGenParticles"),
    PackedGenParticles=cms.InputTag("packedGenParticles"),
    Triggers = cms.InputTag("TriggerResults::HLT"),

    MuonPtCut=cms.double(0),
    IsMC=cms.bool(False),
    SaveTree=cms.bool(True),
    Debug=cms.bool(False)

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