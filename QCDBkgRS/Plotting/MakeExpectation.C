#include <TChain.h>
#include "TProofServ.h"
#include "TProof.h"
#include "TString.h"
//#include "EffMaker.h"

/*
 * Instructions:
 * first time use: make sure that you linked the ~/.proof folder to sonas (or some large storage) sinze the temp files might be very large!
 * check also  https://twiki.cern.ch/twiki/bin/viewauth/CMS/HamburgWikiComputingNAFPOD for more information
 * log in to naf
 * load module use -a /afs/desy.de/group/cms/modulefiles/
 * export SCRAM_ARCH=slc6_amd64_gcc491
 * load moduels:
 * module load root
 * module load pod
 * start the POD server: 
 * pod-server start
 * pod-submit -n 80 -r ge
 * wait, check how many servers are running with:
 * pod-info -n
 * when enough have started
 * change to cmsenv root (needed for the proof to funciton with this trees)
 * cmsenv (in enviroment)
 * run root
 * root -l MakePrediction.C
*/
void MakeExpectation()
{
//  	TString connect = gSystem->GetFromPipe("pod-info -c");
// 	TProof *proof = TProof::Open("adraeger@nafhh-cms03.desy.de:21001");
//      	TProof *proof = TProof::Open("workers=10");
// // 	TProof *proof = TProof::Open("");
   // 	TProof *proof = TProof::Open("");
    	TString connect = gSystem->GetFromPipe("pod-info -c");
    	TProof *proof = TProof::Open(connect);
	// analyse the trees
// 	TChain *Effchain = new TChain("TreeMaker2/PreSelection");
	// 	Effchain->Add("/nfs/dust/cms/user/adraeger/phys14/mc/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/*root");
	// 	Effchain->Add("/nfs/dust/cms/user/adraeger/phys14/mc/WJetsToLNu_HT-200to400/*root");
	// 	Effchain->Add("/nfs/dust/cms/user/adraeger/phys14/mc/WJetsToLNu_HT-400to600/*root");
	// 	Effchain->Add("/nfs/dust/cms/user/adraeger/phys14/mc/WJetsToLNu_HT-600toInf/*root");
// 	TString connect = gSystem.GetFromPipe("pod-info -c");
// 	TProof *proof = TProof::Open(connect);
	TChain *Effchain = new TChain("RA2TreeMaker/PreSelection");
// 	Effchain->Add("/nfs/dust/cms/user/sonnevej/QCD_13TeV_RandS_withreb_finebins_bestMatching_angles_withNeutrinos_everywhere_DeadECALTP_recoPTbins_newsoftware/*.root");
	Effchain->Add("/nfs/dust/cms/user/sonnevej/QCD_13TeV_RandS_withreb_finebins_bestMatching_angles_withNeutrinos_everywhere_DeadECALTP_recoPTbins_newsoftware/*.root");
// 	Effchain->Add("/afs/desy.de/user/a/adraeger/xxl/CMSSW_7_4_15/src/AllHadronicSUSY/QCDBkgRS/test/QCDSmearingClosure_OnMC.root");
// 	Effchain->Add("/nfs/dust/cms/user/sonnevej/QCD_13TeV_GenSmear_finebins_bestMatching_angles_withNeutrinos_everywhere_DeadECALTP_HBHEnoise_newTreeMaker/*.root");
// 	Effchain->Add("/nfs/dust/cms/user/sonnevej/QCD_13TeV_RandS_withreb_finebins_bestMatching_angles_withNeutrinos_everywhere_DeadECALTP_recoPTbins_newsoftware/*300to500*.root");
	

	//    	Effchain->Add("/nfs/dust/cms/user/adraeger/phys14/mc/SMS-T1bbbb_2J_mGl-1000_mLSP-900/*root");
	//   	Effchain->Add("/nfs/dust/cms/user/adraeger/phys14/mc/SMS-T1bbbb_2J_mGl-1500_mLSP-100/*root");
	//   	Effchain->Add("/nfs/dust/cms/user/adraeger/phys14/mc/SMS-T1qqqq_2J_mGl-1000_mLSP-800/*root");
	//   	Effchain->Add("/nfs/dust/cms/user/adraeger/phys14/mc/SMS-T1qqqq_2J_mGl-1400_mLSP-100/*root");
	//    	Effchain->Add("/nfs/dust/cms/user/adraeger/phys14/mc/SMS-T1tttt_2J_mGl-1200_mLSP-800/*root");
	//   	Effchain->Add("/nfs/dust/cms/user/adraeger/phys14/mc/SMS-T1tttt_2J_mGl-1500_mLSP-100/*root");
//  	Effchain->Add("/nfs/dust/cms/user/adraeger/phys14/mc/QCD_HT_1000ToInf_13TeV-madgraph_v1/*root");
//  	Effchain->Add("/nfs/dust/cms/user/adraeger/phys14/mc/QCD_HT_1000ToInf_13TeV-madgraph_v1_ext1/*root");
//  	Effchain->Add("/nfs/dust/cms/user/adraeger/phys14/mc/QCD_HT_250To500_13TeV-madgraph_v1/*root");
//  	Effchain->Add("/nfs/dust/cms/user/adraeger/phys14/mc/QCD_HT_250To500_13TeV-madgraph_v1_ext1/*root");
//  	Effchain->Add("/nfs/dust/cms/user/adraeger/phys14/mc/QCD_HT-500To1000_13TeV-madgraph_v1/*root");
//  	Effchain->Add("/nfs/dust/cms/user/adraeger/phys14/mc/QCD_HT-500To1000_13TeV-madgraph_v1_ext1/*root");
       	Effchain->SetProof();
	//	Effchain->Process("EffMaker.C+g",0,800000);
	Effchain->Process("Expectation.C++g");
//Effchain->Process("Prediction.C++g",0,1000);
//       	Effchain->SetProof(0);
       	delete proof;
				
				
				// plotting
				std::cout<<"Expectation done.\n";
}
