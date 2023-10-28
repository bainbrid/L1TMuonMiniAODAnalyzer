# L1TMuonMiniAODAnalyzer

This is a simple tool for analyzing MINIAOD files for making L1T Muon performance plots very quickly.

Structure is pretty straightforward. The EDAnalyzer just creates a flat ntuple using `slimmedMuons`, `slimmedDisplacedMuons`, and `gmtStage2Digis` while using some trigger information as flags.

Currently we associate offline prompt muons with Single/Double/Triple muons. For SingleMu we use the default SingleMu22/25, for Double/Triple muon objects for now we only use the newly added seeds for Run 3 (LowMass and BPH seeds for DoubleMu, TauTo3Mu seeds for TripleMu). The rest can/will be added later.

# Usage

In principle it should work in any 12_4/5/6_X release. Example using 12_4_8:

```
cmsrel CMSSW_12_4_8
cd CMSSW_12_4_8/src
cmsenv
voms-proxy-init --voms cms --valid 24:00:00
git cms-init

git clone git@github.com:bainbrid/L1TMuonMiniAODAnalyzer.git
scram b -j8

cmsRun L1TMuonMiniAODAnalyzer/pset_files/run3_data.py
```

Of course, for running over big datasets, it's better to use CRAB.  
