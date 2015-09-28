# Expects a file name as argument e.g.
# cmsRun mcresolutions_cfg.py data_set=/QCD_HT-500To1000_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM, globaltag=PHYS14_25_V1

## --- Read parameters --------------------------------------------------
from AllHadronicSUSY.Utils.CommandLineParams import CommandLineParams
parameters = CommandLineParams()

runOnMC = parameters.value("is_mc",True)
dataSetName = parameters.value("dataset_name","/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM")
testFileName = parameters.value("test_filename","/store/mc/RunIISpring15DR74/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/6CE6635E-7016-E511-A22F-0025905938A8.root")
globalTag = parameters.value("global_tag","MCRUN2_74_V9")+"::All"
tagName = parameters.value("tag_name","PAT")
jecFile = parameters.value("jec_file","Summer15_25nsV2_MC")
lumi = parameters.value("lumi",10000.0)
doPUReweighting = parameters.value("pu_reweighting",False)
mcResoJetTag = parameters.value("jet_tag","GoodJets")
mcResoLeptonTag = parameters.value("lepton_tag","GoodLeptons")
mcResoFileName = parameters.value("out_name","MCResoDefault.root")
mtcut = parameters.value("mTcut",100.0)

print "*** JOB SETUP ****************************************************"
print "  is_mc          : "+str(runOnMC)
print "  data_set       : "+dataSetName
print "  global_tag     : "+globalTag
print "  jecFile        : "+jecFile
print "  lumi           : "+str(lumi)
print "  pu_reweighting : "+str(doPUReweighting)
print "  jet_tag        : "+mcResoJetTag
print "  lepton_tag     : "+mcResoLeptonTag
print "  out_name       : "+mcResoFileName
print "******************************************************************"

## --- The process needs to be defined AFTER reading sys.argv, ---------
## --- otherwise edmConfigHash fails -----------------------------------
import FWCore.ParameterSet.Config as cms
import sys,os
process = cms.Process("ResponseTemplates")

## configure geometry & conditions
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = globalTag

## --- Log output ------------------------------------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr = cms.untracked.PSet(placeholder = cms.untracked.bool(True))
process.MessageLogger.cout = cms.untracked.PSet(INFO = cms.untracked.PSet(reportEvery = cms.untracked.int32(1000)))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

## --- Input Source ----------------------------------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(testFileName))
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

## --- Setup WeightProducer --------------------------------------------
from TreeMaker.WeightProducer.getWeightProducer_cff import getWeightProducer
process.WeightProducer = getWeightProducer(dataSetName)
process.WeightProducer.Lumi = cms.double(lumi)
process.WeightProducer.PU = cms.int32(0) # PU: 3 for S10, 2 for S7
process.WeightProducer.FileNamePUDataDistribution = cms.string("NONE")

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

## --- good photon producer -----------------------------------------------
process.GoodPhotons = cms.EDProducer("PhotonIDisoProducer",
                                     photonCollection = cms.untracked.InputTag("slimmedPhotons"),
                                     electronCollection = cms.untracked.InputTag("slimmedElectrons"),
                                     conversionCollection = cms.untracked.InputTag("reducedEgamma","reducedConversions",tagName),
                                     beamspotCollection = cms.untracked.InputTag("offlineBeamSpot","","RECO"),
                                     ecalRecHitsInputTag_EE = cms.InputTag("reducedEgamma","reducedEERecHits"),
                                     ecalRecHitsInputTag_EB = cms.InputTag("reducedEgamma","reducedEBRecHits"),
                                     rhoCollection = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
                                     debug = cms.untracked.bool(False)
                                     )

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


## --- Setup MC Truth templates writer ---------------------------------
from AllHadronicSUSY.MCResolutions.mcresolutions_cfi import MCResolutions
process.MCReso = MCResolutions.clone()
process.MCReso.jetTag = mcResoJetTag
process.MCReso.fileName = mcResoFileName
process.MCReso.leptonTag = mcResoLeptonTag

## --- Setup dump event content ----------------------------------------
process.dump   = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(
                     #process.dump *
                     process.IsolatedElectronTracksVeto *
                     process.IsolatedMuonTracksVeto *
                     process.IsolatedPionTracksVeto *
                     process.patJetCorrFactorsReapplyJEC *
                     process.patJetsReapplyJEC *
                     process.GoodLeptons *
                     process.GoodPhotons *
                     process.GoodJets *
                     process.WeightProducer *
                     process.MCReso
                     )
