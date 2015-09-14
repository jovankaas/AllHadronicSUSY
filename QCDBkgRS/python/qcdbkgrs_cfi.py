import FWCore.ParameterSet.Config as cms

## Default QCD Smearing method
QCDfromSmearing = cms.EDAnalyzer('QCDBkgRS',
	RebalanceJetPt = cms.double(10.),
	RebalanceMode = cms.string('MHT'),
	NSmearedJets = cms.int32(-1),
	SmearedJetPt = cms.double(30.),
	SmearCollection = cms.string('Reco'),
	PtBinEdges_scaling = cms.vdouble(0., 7000.),
	EtaBinEdges_scaling = cms.vdouble(0.0, 5.0),
	AdditionalSmearing = cms.vdouble(0.1),
	LowerTailScaling = cms.vdouble(1.0),
	UpperTailScaling = cms.vdouble(1.0),
	AdditionalSmearing_variation = cms.double(1.0),
	LowerTailScaling_variation = cms.double(1.0),
	UpperTailScaling_variation = cms.double(1.0),
	PtBinEdges = cms.vdouble(0, 20, 30, 50, 80, 120, 170, 230, 300, 380, 470, 570, 680, 800, 1000, 1300, 1700, 2200, 2800, 3500, 4300, 5200, 6500),
   #PtBinEdges = cms.vdouble(0, 5, 10, 15, 20, 25, 30, 50, 80, 120, 170, 230, 300, 380, 470, 570, 680, 800, 1000, 1300, 1700, 2200, 2800, 3500),
   EtaBinEdges = cms.vdouble(0, 0.5, 1.1, 1.7, 2.3, 3.2, 5.0),
   #EtaBinEdges = cms.vdouble(0, 0.3, 0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 2.3, 2.8, 3.2, 4.1, 5.0),
	genjetCollection = cms.InputTag('slimmedGenJets'),
   btagTag = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
   btagCut = cms.double(0.890),
	jetCollection = cms.InputTag('slimmedJets'),
	jetCollection_reb = cms.string('rebalancedJets'),
	jetCollection_smeared = cms.string('smearedJets'),
	genjetCollection_smeared = cms.string('smearedGenJets'),
   leptonTag = cms.InputTag('GoodLeptons'),
   VertexCollection  = cms.InputTag('offlineSlimmedPrimaryVertices'),
	uncertaintyName = cms.string(''),
	InputHisto1_HF = cms.string('hResponse1_HF'),
	InputHisto2_HF = cms.string('hResponse2_HF'),
	InputHisto3p_HF = cms.string('hResponse3p_HF'),
   InputHisto1_NoHF = cms.string('hResponse1_NoHF'),
	InputHisto2_NoHF = cms.string('hResponse2_NoHF'),
	InputHisto3p_NoHF = cms.string('hResponse3p_NoHF'),
   RebalanceCorrectionFile = cms.string('/nfs/dust/cms/user/csander/RA2/AdditionalInputFiles_8TeV/RebalanceCorrection.root'),                              
	ControlPlots = cms.bool(False),
   IsData       = cms.bool(False),
   IsMadgraph   = cms.bool(False),                              
	SmearingFile = cms.string('file.root'),
	OutputFile = cms.string('SmearedAndScaledResponse'),
	NRebin = cms.int32(1),
   weightName = cms.InputTag('WeightProducer','weight','RA2QCDSmearingClosure'),
	absoluteTailScaling = cms.bool(True),
	cleverPrescaleTreating = cms.bool(True),
   useRebalanceCorrectionFactors = cms.bool(False),                              
	A0RMS = cms.double(2.5),
	A1RMS = cms.double(10.0),
   #probExtreme = cms.double(3.e-5),
	probExtreme = cms.double(0.),
   HTSeedTag = cms.InputTag('HT'),
   HTSeedMin = cms.double(300.),
   NJetsSeedTag = cms.InputTag('NJets'),
   NJetsSeedMin = cms.int32(2),
	MHTmin = cms.double(0.),
	MHTmax = cms.double(1500.),
	HTmin = cms.double(0.),
	HTmax = cms.double(5000.),
	NbinsMHT = cms.int32(150),
	NbinsHT = cms.int32(100),
	Ntries = cms.int32(100),
	NJetsSave = cms.int32(3),
   HTSave = cms.double(500.),
   MHTSave = cms.double(0.),
	JetsHTPt = cms.double(30.),
	JetsHTEta = cms.double(2.4),
	JetsMHTPt = cms.double(30.),
	JetsMHTEta = cms.double(5.0),
	JetDeltaMin = cms.vdouble(0.5, 0.5, 0.3),
	MHTcut_low = cms.double(200.),
	MHTcut_medium = cms.double(350.),
	MHTcut_high = cms.double(500.),
	HTcut_low = cms.double(350.),
	HTcut_medium = cms.double(500.),
	HTcut_high = cms.double(800.),
	HTcut_veryhigh = cms.double(1000.),
	HTcut_extremehigh = cms.double(1200.)
)



