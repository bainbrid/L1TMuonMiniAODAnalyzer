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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000))
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(23374) )  

from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
patAlgosToolsTask = getPatAlgosToolsTask(process)


process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(

'/store/data/Run2022C/Muon/MINIAOD/PromptReco-v1/000/357/479/00000/00e393c9-cdda-4d21-b34f-b56715eeac7d.root'


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



process.TFileService = cms.Service("TFileService", fileName = cms.string("l1tMuonNtuple.root") )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag="124X_dataRun3_v9"

process.endjob_step = cms.EndPath(process.endOfProcess)


# import FWCore.PythonUtilities.LumiList as LumiList
# import FWCore.ParameterSet.Types as CfgTypes
# process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())

# myLumis = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')
# process.source.lumisToProcess.extend(myLumis)
# JSONfile = 'Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt'

process.load('L1TMuonMiniAODAnalyzer.L1TMuonMiniAODAnalyzer.L1TMuonMiniAODAnalyzer_cfi')

process.analysis_step = cms.Path(process.L1TMuonMiniAODAnalyzer)

process.schedule = cms.Schedule(process.analysis_step, process.endjob_step)