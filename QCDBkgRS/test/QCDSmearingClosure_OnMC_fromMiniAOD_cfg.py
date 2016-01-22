# For testing do:
# cmsRun QCDSmearingClosure_OnMC_fromMiniAOD_cfg.py

## --- Read parameters --------------------------------------------------
from TreeMaker.Utils.CommandLineParams import CommandLineParams
parameters = CommandLineParams()

runOnMC = parameters.value("is_mc",True)
dataSetName = parameters.value("dataset_name","/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM")
testFileName = parameters.value("test_filename","/store/mc/RunIISpring15DR74/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/6CE6635E-7016-E511-A22F-0025905938A8.root")
globalTag = parameters.value("global_tag","MCRUN2_74_V9")+"::All"
tagName = parameters.value("tag_name","PAT")
jecFile = parameters.value("jec_file","Summer15_25nsV2_MC")
lumi = parameters.value("lumi",10000.0)
doPUReweighting = parameters.value("pu_reweighting",False)
InputJetTag = parameters.value("jet_tag","GoodJets")
InputLeptonTag = parameters.value("lepton_tag","GoodLeptons")
OutputFileName = parameters.value("out_name","QCDSmearingClosure_OnMC.root")
mtcut = parameters.value("mTcut",100.0)

print "*** JOB SETUP ****************************************************"
print "  is_mc          : "+str(runOnMC)
print "  data_set       : "+dataSetName
print "  global_tag     : "+globalTag
print "  jecFile        : "+jecFile
print "  lumi           : "+str(lumi)
print "  pu_reweighting : "+str(doPUReweighting)
print "  jet_tag        : "+InputJetTag
print "  lepton_tag     : "+InputLeptonTag
print "  out_name       : "+OutputFileName
print "******************************************************************"

## --- The process needs to be defined AFTER reading sys.argv, ---------
## --- otherwise edmConfigHash fails -----------------------------------
import FWCore.ParameterSet.Config as cms
import sys,os
process = cms.Process("RA2QCDSmearingClosure")

## configure geometry & conditions
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = globalTag

## --- Log output ------------------------------------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr = cms.untracked.PSet(placeholder = cms.untracked.bool(True))
process.MessageLogger.cout = cms.untracked.PSet(INFO = cms.untracked.PSet(reportEvery = cms.untracked.int32(100)))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

## --- For genlevel info ------------------------------------------------
##process.options.SkipEvent = cms.untracked.vstring('ProductNotFound') 

## --- Input Source ----------------------------------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(testFileName))
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

## --- Setup WeightProducer --------------------------------------------
print "*** WeightProducer setup **************************************************"
from TreeMaker.WeightProducer.getWeightProducer_cff import getWeightProducer
process.WeightProducer = getWeightProducer(dataSetName)
process.WeightProducer.Lumi = cms.double(lumi)
process.WeightProducer.PU = cms.int32(0) # PU: 3 for S10, 2 for S7
process.WeightProducer.FileNamePUDataDistribution = cms.string("NONE")

###############################################################################
process.TFileService = cms.Service("TFileService",fileName = cms.string("QCDSmearingClosure_OnMC.root") )
###############################################################################

###############################################################################
process.load("AllHadronicSUSY.TruthNoiseFilter.truthnoisefilter_cfi")
process.TruthNoiseFilter.jetCollection = InputJetTag
###############################################################################

###############################################################################
process.load("AllHadronicSUSY.QCDBkgRS.qcdbkgrs_cfi")
###############################################################################

###############################################################################
# Rebalancing and Smearing configuration
###############################################################################
print "*** R+S Configuration **************************************************"
#process.QCDfromSmearing.SmearingFile = '/afs/desy.de/user/s/sonnevej/xxl/CMSSW_7_4_6_patch6/src/AllHadronicSUSY/MCResolutions/data/QCD_13TeV_MGMLM_Spring15_fineBins_bestMatching_DeadECALTP.root'
process.QCDfromSmearing.SmearingFile = 'QCD_13TeV_madgraph-MLM_finebins_bestMatching_angles_NeutrinosEverywhere_NeutrinosInBins_DeadECALTP.root'
#process.QCDfromSmearing.SmearingFile = '/afs/desy.de/user/s/sonnevej/xxl/CMSSW_7_4_6_patch6/src/AllHadronicSUSY/MCResolutions/data/QCD_13TeV_madgraph-MLM_finebins_bestMatching_angles_NeutrinosEverywhere_NeutrinosInBins_DeadECALTP_recoPTbins.root'
#process.QCDfromSmearing.SmearingFile = '/afs/desy.de/user/s/sonnevej/xxl/CMSSW_7_4_6_patch6/src/AllHadronicSUSY/MCResolutions/data/QCD_13TeV_madgraph-MLM_finebins_bestMatching_angles_NeutrinosEverywhere_NeutrinosInBins_DeadECALTP_recoPTbins.root'
process.QCDfromSmearing.jetCollection = InputJetTag
process.QCDfromSmearing.leptonTag = InputLeptonTag
process.QCDfromSmearing.uncertaintyName = ''
process.QCDfromSmearing.InputHistoPt_HF = 'h_b_JetAll_ResponsePt'
process.QCDfromSmearing.InputHistoEta_HF = 'h_b_JetAll_ResponseEta'
process.QCDfromSmearing.InputHistoPhi_HF = 'h_b_JetAll_ResponsePhi'
process.QCDfromSmearing.InputHistoPt_NoHF = 'h_nob_JetAll_ResponsePt'
process.QCDfromSmearing.InputHistoEta_NoHF = 'h_nob_JetAll_ResponseEta'
process.QCDfromSmearing.InputHistoPhi_NoHF = 'h_nob_JetAll_ResponsePhi'
#process.QCDfromSmearing.RebalanceCorrectionFile = '/afs/desy.de/user/s/sonnevej/dust/RA2b_input/RebalanceCorrectionFactors_madgraph_spring15_withoutPUReweighting_withNeutrinosEverywhere_DeadECALTP_recoPTbinspt10.root'
process.QCDfromSmearing.RebalanceCorrectionFile = 'RebalanceCorrectionFactors_madgraph_spring15_withoutPUReweighting_withNeutrinosEverywhere.root'
#'/nfs/dust/cms/user/csander/RA2/AdditionalInputFiles_13TeV/RebalanceCorrectionFactors_madgraph_spring15_withoutPUReweighting_pt10.root'
#process.QCDfromSmearing.BTagEfficiencyFile = '/afs/desy.de/user/s/sonnevej/dust/RA2b_input/B_Mis_TagEfficiencies_Spring15MadGraph_DeadECALTP_recoPTbins.root'
process.QCDfromSmearing.BTagEfficiencyFile = 'B_Mis_TagEfficiencies_Neutrinos_everywhere_Spring15MadGraph.root'
process.QCDfromSmearing.NRebin = 1
#process.QCDfromSmearing.SmearCollection = 'Reco'
process.QCDfromSmearing.SmearCollection = 'Gen'
process.QCDfromSmearing.PtBinEdges_scaling = cms.vdouble(0., 7000.)
process.QCDfromSmearing.EtaBinEdges_scaling = cms.vdouble(0.0, 5.0)
process.QCDfromSmearing.AdditionalSmearing = cms.vdouble(1.0)
process.QCDfromSmearing.LowerTailScaling = cms.vdouble(1.0)
process.QCDfromSmearing.UpperTailScaling = cms.vdouble(1.0)
process.QCDfromSmearing.SmearedJetPt = 0.
process.QCDfromSmearing.RebalanceJetPt = 10.
process.QCDfromSmearing.RebalanceMode = 'MHThigh'
#process.QCDfromSmearing.RebalanceMode = 'MHTall'
process.QCDfromSmearing.weightName = 'WeightProducer:weight:RA2QCDSmearingClosure'
process.QCDfromSmearing.ControlPlots = True
process.QCDfromSmearing.HTSeedMin = 0.
process.QCDfromSmearing.NJetsSeedMin = 2
process.QCDfromSmearing.Ntries = 100
###default for GenSmearing or full closure tests
process.QCDfromSmearing.NJetsSave = 3
process.QCDfromSmearing.HTSave = cms.double(500.)
process.QCDfromSmearing.MHTSave = cms.double(0.)
###default for full closure tests with clever prescale treating
#process.QCDfromSmearing.NJetsSave = 4
#process.QCDfromSmearing.HTSave = cms.double(500.)
#process.QCDfromSmearing.MHTSave = cms.double(200.)
process.QCDfromSmearing.cleverPrescaleTreating = False
process.QCDfromSmearing.useRebalanceCorrectionFactors = False
process.QCDfromSmearing.useBTagEfficiencyFactors = True
process.QCDfromSmearing.useCleverRebalanceCorrectionFactors = True
process.QCDfromSmearing.MHTcut_low = cms.double(200.)
process.QCDfromSmearing.MHTcut_medium = cms.double(350.)
process.QCDfromSmearing.MHTcut_high = cms.double(500.)
process.QCDfromSmearing.HTcut_low = cms.double(500.)
process.QCDfromSmearing.HTcut_medium = cms.double(800.)
process.QCDfromSmearing.HTcut_high = cms.double(1000.)
process.QCDfromSmearing.HTcut_veryhigh = cms.double(1200.)
process.QCDfromSmearing.HTcut_extremehigh = cms.double(1400.)
###############################################################################

VarsInt = cms.vstring()
VectorInt = cms.vstring()
VarsDouble = cms.vstring()
VarsBool = cms.vstring()
VectorDouble = cms.vstring()
RecoCandVector = cms.vstring()

# baseline producers
process.Baseline = cms.Sequence(
)

## --- NVertex producer -----------------------------------------------
print "*** NVertex producer **************************************************"
from TreeMaker.Utils.primaryvertices_cfi import primaryvertices
process.NVtx = primaryvertices.clone(
                                      VertexCollection  = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                      )
process.Baseline += process.NVtx
VarsInt.extend(['NVtx'])

## --- isotrack producer -----------------------------------------------
from TreeMaker.Utils.trackIsolationMaker_cfi import trackIsolationFilter

process.IsolatedElectronTracksVeto = trackIsolationFilter.clone(
                                                                doTrkIsoVeto= True,
                                                                vertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                                                pfCandidatesTag = cms.InputTag("packedPFCandidates"),
                                                                dR_ConeSize         = cms.double(0.3),
                                                                dz_CutValue         = cms.double(0.1),
                                                                minPt_PFCandidate   = cms.double(5.0),
                                                                isoCut              = cms.double(0.2),
                                                                pdgId               = cms.int32(11),
                                                                mTCut=mtcut,
                                                                )

process.IsolatedMuonTracksVeto = trackIsolationFilter.clone(
                                                            doTrkIsoVeto= True,
                                                            vertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                                            pfCandidatesTag = cms.InputTag("packedPFCandidates"),
                                                            dR_ConeSize         = cms.double(0.3),
                                                            dz_CutValue         = cms.double(0.1),
                                                            minPt_PFCandidate   = cms.double(5.0),
                                                            isoCut              = cms.double(0.2),
                                                            pdgId               = cms.int32(13),
                                                            mTCut=mtcut,
                                                            )

process.IsolatedPionTracksVeto = trackIsolationFilter.clone(
                                                            doTrkIsoVeto= True,
                                                            vertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                                            pfCandidatesTag = cms.InputTag("packedPFCandidates"),
                                                            dR_ConeSize         = cms.double(0.3),
                                                            dz_CutValue         = cms.double(0.1),
                                                            minPt_PFCandidate   = cms.double(10.0),
                                                            isoCut              = cms.double(0.1),
                                                            pdgId               = cms.int32(211),
                                                            mTCut=mtcut,
                                                            )

process.Baseline += process.IsolatedElectronTracksVeto
process.Baseline += process.IsolatedMuonTracksVeto
process.Baseline += process.IsolatedPionTracksVeto

## --- good leptons producer -------------------------------------------
from TreeMaker.Utils.leptonproducer_cfi import leptonproducer
process.GoodLeptons = leptonproducer.clone(
                                           MuonTag          = cms.InputTag('slimmedMuons'),
                                           ElectronTag      = cms.InputTag('slimmedElectrons'),
                                           PrimaryVertex    = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                           minElecPt        = cms.double(10),
                                           maxElecEta       = cms.double(2.5),
                                           minMuPt          = cms.double(10),
                                           maxMuEta         = cms.double(2.4),
                                           UseMiniIsolation = cms.bool(True),
                                           muIsoValue       = cms.double(0.2),
                                           elecIsoValue     = cms.double(0.1), # only has an effect when used with miniIsolation
                                           METTag           = cms.InputTag('slimmedMETs'),
                                           )
process.Baseline += process.GoodLeptons
VarsInt.extend(['GoodLeptons'])

## --- good photon producer -----------------------------------------------
process.GoodPhotons = cms.EDProducer("PhotonIDisoProducer",
                                     photonCollection = cms.untracked.InputTag("slimmedPhotons"),
                                     electronCollection = cms.untracked.InputTag("slimmedElectrons"),
                                     conversionCollection = cms.untracked.InputTag("reducedEgamma","reducedConversions",tagName),
                                     beamspotCollection = cms.untracked.InputTag("offlineBeamSpot","","RECO"),
                                     ecalRecHitsInputTag_EE = cms.InputTag("reducedEgamma","reducedEERecHits"),
                                     ecalRecHitsInputTag_EB = cms.InputTag("reducedEgamma","reducedEBRecHits"),
                                     rhoCollection = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                                     genParCollection = cms.untracked.InputTag("prunedGenParticles"), 
                                     debug = cms.untracked.bool(False)
                                     )
process.Baseline += process.GoodPhotons

## ----------------------------------------------------------------------------------------------
## JECs
## ----------------------------------------------------------------------------------------------

# get the JECs (disabled by default)
# this requires the user to download the .db file from this twiki
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
if len(jecFile)>0:
   JECPatch = cms.string('sqlite_file:'+jecFile+'.db')
   if os.getenv('GC_CONF'):
      JECPatch = cms.string('sqlite_file:/nfs/dust/cms/user/csander/RA2/AdditionalInputFiles_13TeV/'+jecFile+'.db')

   process.load("CondCore.DBCommon.CondDBCommon_cfi")
   from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup
   process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
                              connect = JECPatch,
                              toGet   = cms.VPSet(
                                                  cms.PSet(
                                                           record = cms.string("JetCorrectionsRecord"),
                                                           tag    = cms.string("JetCorrectorParametersCollection_"+jecFile+"_AK4PFchs"),
                                                           label  = cms.untracked.string("AK4PFchs")
                                                           ),
                                                  cms.PSet(
                                                           record = cms.string("JetCorrectionsRecord"),
                                                           tag    = cms.string("JetCorrectorParametersCollection_"+jecFile+"_AK4PF"),
                                                           label  = cms.untracked.string("AK4PF")
                                                           )
                                                  )
                              )
   process.es_prefer_jec = cms.ESPrefer("PoolDBESSource","jec")

   from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
   process.patJetCorrFactorsReapplyJEC = patJetCorrFactorsUpdated.clone(
                                                                        src     = cms.InputTag("slimmedJets"),
                                                                        levels  = ['L1FastJet',
                                                                                   'L2Relative',
                                                                                   'L3Absolute'],
                                                                        payload = 'AK4PFchs' # Make sure to choose the appropriate levels and payload here!
                                                                        )
                                                                        #   if residual: process.patJetCorrFactorsReapplyJEC.levels.append('L2L3Residual')
                                                                        
   from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
   process.patJetsReapplyJEC = patJetsUpdated.clone(
                                                    jetSource = cms.InputTag("slimmedJets"),
                                                    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
                                                    )

## --- good jets producer -----------------------------------------------
from TreeMaker.Utils.goodjetsproducer_cfi import GoodJetsProducer
process.GoodJets = GoodJetsProducer.clone(
                                          TagMode = cms.bool(False),
                                          #JetTag= cms.InputTag('slimmedJets'),
                                          JetTag= cms.InputTag('patJetsReapplyJEC'),
                                          maxJetEta = cms.double(5.0),
                                          maxMuFraction = cms.double(2),
                                          minNConstituents = cms.double(2),
                                          maxNeutralFraction = cms.double(0.90),
                                          maxPhotonFraction = cms.double(0.95),
                                          minChargedMultiplicity = cms.double(0),
                                          minChargedFraction = cms.double(0),
                                          maxChargedEMFraction = cms.double(0.99),
                                          jetPtFilter = cms.double(30),
                                          ExcludeLepIsoTrackPhotons = cms.bool(True),
                                          JetConeSize = cms.double(0.4),
                                          MuonTag = cms.InputTag('GoodLeptons:IdIsoMuon'),
                                          ElecTag = cms.InputTag('GoodLeptons:IdIsoElectron'),
                                          IsoElectronTrackTag = cms.InputTag('IsolatedElectronTracksVeto'),
                                          IsoMuonTrackTag = cms.InputTag('IsolatedMuonTracksVeto'),
                                          IsoPionTrackTag = cms.InputTag('IsolatedPionTracksVeto'),
                                          PhotonTag = cms.InputTag('GoodPhotons','bestPhoton'),
                                          ### TEMPORARY ###
                                          VetoHF = cms.bool(False),
                                          VetoEta = cms.double(3.0)
                                          )
process.Baseline += process.GoodJets
VarsBool.extend(['GoodJets:JetID(JetID)'])

## --- HT jets producer ------------------------------------------------
print "*** HT jets producer **************************************************"
from TreeMaker.Utils.subJetSelection_cfi import SubJetSelection
process.HTJets = SubJetSelection.clone(
                                       JetTag   = cms.InputTag('GoodJets'),
                                       MinPt    = cms.double(30),
                                       MaxEta   = cms.double(2.4),
)
process.Baseline += process.HTJets

## --- HT producer -----------------------------------------------------
print "*** HT producer **************************************************"
from TreeMaker.Utils.htdouble_cfi import htdouble
process.HT = htdouble.clone(
                            JetTag  = cms.InputTag('HTJets'),
)
process.Baseline += process.HT
VarsDouble.extend(['HT'])
   
## --- NJets producer --------------------------------------------------
print "*** NJets producer **************************************************"
from TreeMaker.Utils.njetint_cfi import njetint
process.NJets = njetint.clone(
                              JetTag  = cms.InputTag('HTJets'),
)
process.Baseline += process.NJets
VarsInt.extend(['NJets'])

## --- NBJets producer -------------------------------------------------
print "*** NBJets producer **************************************************"
from TreeMaker.Utils.btagint_cfi import btagint
process.BTags = btagint.clone(
                              JetTag         = cms.InputTag('HTJets'),
                              BTagInputTag = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
                              BTagCutValue = cms.double(0.890)
)
process.Baseline += process.BTags
VarsInt.extend(['BTags'])

## --- MHT jets producer -----------------------------------------------
print "*** MHT jets producer **************************************************"
from TreeMaker.Utils.subJetSelection_cfi import SubJetSelection
process.MHTJets = SubJetSelection.clone(
                                        JetTag  = cms.InputTag('GoodJets'),
                                        MinPt   = cms.double(30),
                                        MaxEta  = cms.double(5.0),
                                        )
process.Baseline += process.MHTJets

## --- MHT producer ----------------------------------------------------
print "*** MHT producer **************************************************"
from TreeMaker.Utils.metdouble_cfi import metdouble
process.MET = metdouble.clone(
                              JetTag  = cms.InputTag('MHTJets'),
                              )
process.Baseline += process.MET
VarsDouble.extend(['MET:minDeltaPhiN'])

print "*** MHT producer **************************************************"
from TreeMaker.Utils.mhtdouble_cfi import mhtdouble
process.MHT = mhtdouble.clone(
                              JetTag  = cms.InputTag('MHTJets'),
)
process.Baseline += process.MHT
VarsDouble.extend(['MHT:Pt(MHT)'])

## --- DeltaPhi producer -----------------------------------------------
print "*** DeltaPhi producer **************************************************"
from TreeMaker.Utils.deltaphidouble_cfi import deltaphidouble
process.DeltaPhi = deltaphidouble.clone(
                                        DeltaPhiJets  = cms.InputTag('HTJets'),
                                        MHTJets       = cms.InputTag('MHTJets'),
)
process.Baseline += process.DeltaPhi
VarsDouble.extend(['DeltaPhi:Jet1Pt','DeltaPhi:Jet2Pt','DeltaPhi:Jet3Pt','DeltaPhi:Jet4Pt'])
VarsDouble.extend(['DeltaPhi:Jet1Eta','DeltaPhi:Jet2Eta','DeltaPhi:Jet3Eta','DeltaPhi:Jet4Eta'])
VarsDouble.extend(['DeltaPhi:DeltaPhi1','DeltaPhi:DeltaPhi2','DeltaPhi:DeltaPhi3'])

## --- Setup of TreeMaker ----------------------------------------------

VarsDouble.extend(['WeightProducer:weight(Weight)'])

print "*** Treemaker setup **************************************************"
from TreeMaker.TreeMaker.treeMaker import TreeMaker
process.RA2TreeMaker = TreeMaker.clone(
                                       TreeName       = cms.string("PreSelection"),
                                       VarsRecoCand   = RecoCandVector,
                                       VarsDouble     = VarsDouble,
                                       VarsBool     = VarsBool,
                                       VectorDouble   = VectorDouble,
                                       VarsInt        = VarsInt,
                                       VectorInt      = VectorInt,
)
###############################################################################

### Very simple tail and core scalings for uncertainty checks
process.QCDfromSmearingCoreUP = process.QCDfromSmearing.clone()
process.QCDfromSmearingCoreUP.AdditionalSmearing = cms.vdouble(1.1)
process.QCDfromSmearingCoreUP.uncertaintyName = 'CoreUP'
process.QCDfromSmearingCoreDN = process.QCDfromSmearing.clone()
process.QCDfromSmearingCoreDN.AdditionalSmearing = cms.vdouble(0.9)
process.QCDfromSmearingCoreDN.uncertaintyName = 'CoreDN'
process.QCDfromSmearingTailUP = process.QCDfromSmearing.clone()
process.QCDfromSmearingTailUP.LowerTailScaling = cms.vdouble(1.2)
process.QCDfromSmearingTailUP.UpperTailScaling = cms.vdouble(1.2)
process.QCDfromSmearingTailUP.uncertaintyName = 'TailUP'
process.QCDfromSmearingTailDN = process.QCDfromSmearing.clone()
process.QCDfromSmearingTailDN.LowerTailScaling = cms.vdouble(0.8)
process.QCDfromSmearingTailDN.UpperTailScaling = cms.vdouble(0.8)
process.QCDfromSmearingTailDN.uncertaintyName = 'TailDN'

###############################################################################
process.dump   = cms.EDAnalyzer("EventContentAnalyzer")
###############################################################################

## --- MET Filters -----------------------------------------------
from RecoMET.METFilters.metFilters_cff import EcalDeadCellTriggerPrimitiveFilter
process.ECALDeadCellFilter = EcalDeadCellTriggerPrimitiveFilter.clone()
process.prediction = cms.Path(
                              process.ECALDeadCellFilter *
                              process.patJetCorrFactorsReapplyJEC *
                              process.patJetsReapplyJEC *
                              process.Baseline *
                              #process.TruthNoiseFilter *
                              process.WeightProducer *
                              process.QCDfromSmearing
                              #* process.dump
)

#process.predictionCoreUP = cms.Path(
#                                    process.patJetCorrFactorsReapplyJEC *
#                                    process.patJetsReapplyJEC *
#                                    process.Baseline *
#                                    process.WeightProducer *
#                                    process.QCDfromSmearingCoreUP
#                                    )

#process.predictionCoreDN = cms.Path(
#                                    process.patJetCorrFactorsReapplyJEC *
#                                    process.patJetsReapplyJEC *
#                                    process.Baseline *
#                                    process.WeightProducer *
#                                    process.QCDfromSmearingCoreDN
#                                    )

#process.predictionTailUP = cms.Path(
#                                    process.patJetCorrFactorsReapplyJEC *
#                                    process.patJetsReapplyJEC *
#                                    process.Baseline *
#                                    process.WeightProducer *
#                                    process.QCDfromSmearingTailUP
#                                    )

#process.predictionTailDN = cms.Path(
#                                    process.patJetCorrFactorsReapplyJEC *
#                                    process.patJetsReapplyJEC *
#                                    process.Baseline *
#                                    process.WeightProducer *
#                                    process.QCDfromSmearingTailDN
#                                    )

process.mc = cms.Path(
                      process.patJetCorrFactorsReapplyJEC *
                      process.patJetsReapplyJEC *
                      process.Baseline *
                      #process.TruthNoiseFilter *
                      process.WeightProducer *
                      process.RA2TreeMaker
                      #* process.dump
)
