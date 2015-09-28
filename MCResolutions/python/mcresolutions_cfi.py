import FWCore.ParameterSet.Config as cms

MCResolutions = cms.EDAnalyzer('MCResolutions',
   leptonTag         = cms.InputTag('GoodLeptons'),
	jetTag				= cms.InputTag('GoodJets'),
   btagTag           = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
   btagCut           = cms.double(0.890),
	genJetTag			= cms.InputTag('slimmedGenJets'),
	weightName			= cms.InputTag('WeightProducer','weight','ResponseTemplates'),
	deltaRMatch			= cms.double(0.25),
	deltaRMatchVeto	= cms.double(0.7),
	absPtVeto			= cms.double(20.),
	relPtVeto			= cms.double(0.01),
	GenJetPtCut			= cms.double(0.),
   fileName				= cms.string('MCJetResolution.root')
)
