import FWCore.ParameterSet.Config as cms

process = cms.Process("L1TMuonAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
patAlgosToolsTask = getPatAlgosToolsTask(process)


process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(

        # Test file
        #'/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2810000/0009bf67-3c7a-4e81-b50a-e3914b3d2ffa.root'

        # 10 files
        '/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/f365faac-9a86-45cc-b1b4-87f036b1a2fa.root',
        '/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/fad35f5f-e549-40d0-a48f-d2b5e2086041.root',
        '/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/c95208fb-fa54-40a4-a000-5a65877fb3c2.root',
        '/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/7019e6c1-2adb-4fb0-b04e-f97449a0dfa4.root',
        '/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/9d7dd0b0-a52b-41a7-8197-fe0f7d051af2.root',
        '/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/678a8d19-b4d8-48be-89d6-c98bca514b69.root',
        '/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/839d6b3e-e85b-4c44-9c64-8632971e300f.root',
        '/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/1acebde6-6bee-464f-8775-2b161f1380ec.root',
        '/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/5351effc-3ca0-44f0-9acc-e1e355860c43.root',
        '/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/23aa49a6-a411-4e54-9212-f5ec90c0cb39.root',

        # 10 files
        #'/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/163eb519-49a0-4dc0-9aad-55a53add0a5e.root',
        #'/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/209cd86b-b3a3-4d8f-88ff-1916399e27a7.root',
        #'/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/f2432d3e-b8c2-48cd-a206-d0bd38507d4b.root',
        #'/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/4c114f30-acb4-44d6-9390-34991c4d1f5a.root',
        #'/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/a3285029-4cde-42f2-9f87-71be739778a8.root',
        #'/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/48c83db8-b0f1-4343-b2e4-e85edd5c6fba.root',
        #'/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/0502f452-0df8-4d85-890d-443e015e31eb.root',
        #'/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/c7590f38-fbbf-473f-ac12-6a39ae090b79.root',
        #'/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/e8c46d5e-ba1f-411a-a0c4-ad8792f34a10.root',
        #'/store/mc/Run3Summer22MiniAODv3/InclusiveDileptonMinBias_TuneCP5Plus_13p6TeV_pythia8/MINIAODSIM/Pilot_124X_mcRun3_2022_realistic_v12-v5/2820000/033d3edc-1ccd-4817-889b-42b051cbcd93.root',
                                )
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(8),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("L1TMuonNtuple_run3_mc.root") )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag="124X_dataRun3_v9"

process.endjob_step = cms.EndPath(process.endOfProcess)

process.load('L1TMuonMiniAODAnalyzer.L1TMuonMiniAODAnalyzer.L1TMuonMiniAODAnalyzerFlat_cfi')

process.analysis_step = cms.Path(process.L1TMuonMiniAODAnalyzerFlat)

process.schedule = cms.Schedule(process.analysis_step, process.endjob_step)
