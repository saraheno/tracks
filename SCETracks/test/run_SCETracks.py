import FWCore.ParameterSet.Config as cms

process = cms.Process("RUN")

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring("file:/data/users/eno/CMSSW_7_6_1/src/step2.root"))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(5))

process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring("cout"),
    cout = cms.untracked.PSet(threshold = cms.untracked.string("ERROR")))

process.SCETracks = cms.EDAnalyzer("SCETracks")

process.path = cms.Path(process.SCETracks)
