#include <TChain.h>
#include "TProofServ.h"
#include "TProof.h"
#include "TString.h"
//#include "EffMaker.h"

/*
 * Instructions:
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

using namespace std;

TCanvas* DrawComparison(TH1F* prediction, TH1F* selection, TString Title, TString LumiTitle, TString xTitle, TString yTitle, bool isData)
{
	double MinX = selection->GetXaxis()->GetXmin();
	double MaxX = selection->GetXaxis()->GetXmax();
	double MaxY = selection->GetBinContent(selection->GetMaximumBin());
	double YRangeMax = 2*pow(10,(int(log10(MaxY))/1)+2);
	double MinY = selection->GetBinContent(selection->GetMinimumBin());
	if (MinY < 0.01) MinY = 0.01;
	double YRangeMin = 0.5*pow(10,(int(log10(MinY))/1)-1);
	TString titlePrediction;
	TString titleSelection;
	TString RatioTitle;
	
	if( isData ){
		titlePrediction = "Pred. from Data";
		titleSelection = "Data";
		RatioTitle = "(Pred-Data)/Data";
	}
	else {
		//titlePrediction = "Data-driven Pred. from MC";
		titlePrediction = "Smeared Generator Jets";
		titleSelection = "MC Expectation";
		//RatioTitle = "(Pred-MC)/MC";
		RatioTitle = "(Gen-MC)/MC";
	}
	
	static Int_t c_LightBrown   = TColor::GetColor( "#D9D9CC" );
	static Int_t c_LightGray    = TColor::GetColor( "#DDDDDD" );
	
	selection->SetAxisRange(MinX, MaxX, "X");
	selection->GetYaxis()->SetRangeUser(0.005, YRangeMax);
	selection->SetMarkerStyle(20);
	selection->SetMarkerSize(0.9);
	selection->SetMarkerColor(kBlack);
	selection->SetXTitle(xTitle);
	selection->SetYTitle(yTitle);
	prediction->SetAxisRange(MinX, MaxX, "X");
	prediction->GetYaxis()->SetRangeUser(YRangeMin, YRangeMax);
	// prediction->SetFillColor(c_LightBrown);
	prediction->SetFillColor(c_LightGray);
	prediction->SetTitle("");
	prediction->SetXTitle(xTitle);
	prediction->SetYTitle(yTitle);
	
	TCanvas *c = new TCanvas(prediction->GetName(), "Comparison and ratio of two histos", 700, 700);
	
	TPad *pad1 = new TPad("pad1a", "pad1a", 0, 0.35, 1, 1);
	pad1->SetLogy();
	pad1->SetBottomMargin(0);
	pad1->Draw();
	pad1->cd();
	
	prediction->DrawCopy("hist");
	selection->Draw("esame");
	prediction->SetFillColor(kAzure-3);
	prediction->SetFillStyle(3354);
	prediction->DrawCopy("e2same");
	
	prediction->SetFillStyle(1001);
	//  prediction->SetFillColor(c_LightBrown);
	prediction->SetFillColor(c_LightGray);
	
	//  TLegend* leg1 = new TLegend(0.48, 0.63, 0.95, 0.83);
	TLegend* leg1 = new TLegend(0.44, 0.63, 0.91, 0.83);
	leg1->SetFillStyle(0);
	leg1->SetLineStyle(1);
	leg1->SetTextFont(42);
	// leg1->SetTextSize(0.04);
	leg1->SetTextSize(0.045);
	leg1->AddEntry(prediction, titlePrediction, "lf");
	leg1->AddEntry(selection, titleSelection, "lep");
	leg1->Draw("same");
	
	TPaveText* pt = new TPaveText(0.11, 0.98, 0.95, 0.86, "NDC");
	pt->SetBorderSize(0);
	pt->SetFillStyle(0);
	pt->SetTextAlign(12);
	pt->SetTextSize(0.045);
	pt->AddText(Title);
	pt->AddText(LumiTitle);
	pt->Draw();
	
	c->cd();
	TPad *pad2 = new TPad("pad2a", "pad2a", 0, 0, 1, 0.35);
	pad2->SetTopMargin(0);
	pad2->Draw();
	pad2->cd();
	TH1F* r = new TH1F(*prediction);
	r->SetTitle("");
	r->SetLabelSize(0.08, "XYZ");
	r->SetLabelOffset(0.01, "XYZ");
	// r->SetTitleSize(0.09, "XYZ");
	r->SetTitleSize(0.125, "XYZ");
	r->SetTitleOffset(0.95, "X");
	r->SetTitleOffset(0.53, "Y");
	// r->SetTitleOffset(0.65, "Y");
	r->SetTickLength(0.05);
	r->SetYTitle(RatioTitle);
	r->SetStats(0);
	r->SetMarkerStyle(20);
	r->SetMarkerSize(0.9);
	r->SetMarkerColor(kBlack);
	r->Reset();
	r->Add(prediction, 1);
	r->Add(selection, -1);
	r->Divide(selection);
	r->SetMaximum(1.2);
	r->SetMinimum(-1.2);
	r->Draw("ep");
	TLine l;
	l.DrawLine(MinX, 0., MaxX, 0.);
	c->cd();
	return c;
}

void MakePrediction()
{
//  	TString connect = gSystem->GetFromPipe("pod-info -c");
// 	TProof *proof = TProof::Open("adraeger@nafhh-cms03.desy.de:21001");
//      	TProof *proof = TProof::Open("workers=10");
// // 	TProof *proof = TProof::Open("");
  //  	TProof *proof = TProof::Open("");
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
	TChain *Effchain = new TChain("QCDfromSmearing/QCDPrediction");
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
	Effchain->Process("Prediction.C++g");
//Effchain->Process("Prediction.C++g",0,1000);
//       	Effchain->SetProof(0);
       	delete proof;
				
				
				// plotting
				std::cout<<"Prediction done plotting starting...\n";
				
				gROOT->SetStyle("Plain");
				//gStyle->SetPalette(51, 0);
				//gStyle->SetHatchesLineWidth(1.2);
				
				// For the canvas:
				gStyle->SetCanvasColor(0);
				//gStyle->SetCanvasBorderMode(0);
				
				// For the Pad:
				gStyle->SetPadColor(0);
				gStyle->SetPadTickX(1);
				gStyle->SetPadTickY(1);
				gStyle->SetPadBorderSize(2);
				//gStyle->SetPadBorderMode(0);
				
				// For the frame:
				gStyle->SetFrameBorderMode(0);
				
				// For the histo:
				// gStyle->SetMarkerSize(0.7);
				// gStyle->SetMarkerStyle(20);
				// gStyle->SetMarkerColor(1);
				
				// For the statistics box:
				gStyle->SetOptStat(0);
				//gStyle->SetOptFit(1011);
				
				// Margins:
				gStyle->SetPadBottomMargin(0.25);
				gStyle->SetPadTopMargin(0.15);
				gStyle->SetPadLeftMargin(0.15);
				gStyle->SetPadRightMargin(0.1);
				
				// For the Global title:
				gStyle->SetOptTitle(0);
				gStyle->SetTitleColor(1);
				gStyle->SetTitleFillColor(10);
				gStyle->SetTitleTextColor(1);
				gStyle->SetTitleFont(42);
				gStyle->SetTitleFontSize(0.05);
				gStyle->SetTitleBorderSize(0);
				
				// For the axis
				gStyle->SetNdivisions(510, "X");
				gStyle->SetNdivisions(510, "Y");
				gStyle->SetTickLength(0.03);
				
				// For the axis titles:
				gStyle->SetTitleOffset(1.4, "X");
				//gStyle->SetTitleOffset(1.25, "Y");
				gStyle->SetTitleOffset(1.2, "Y");
				gStyle->SetTitleOffset(0.5, "Z");
				// gStyle->SetTitleSize(0.05, "XYZ");
				gStyle->SetTitleSize(0.061, "XYZ");
				gStyle->SetTitleFont(42, "XYZ");
				//gStyle->SetTitleX(0.15);
				//gStyle->SetTitleY(0.99);
				
				// For the axis labels:
				gStyle->SetLabelSize(0.04, "XYZ");
				gStyle->SetLabelOffset(0.01, "XYZ");
				gStyle->SetLabelFont(42, "XYZ");
				
				// For the legend
				gStyle->SetLegendBorderSize(0);
				
				gROOT->ForceStyle();
				bool isData = false;
				TString LumiTitle;
				TString yTitle = "Events";
				if( isData ) LumiTitle = "CMS preliminary, L = x.yz fb^{  -1}, #sqrt{s} = 13 TeV";
				else LumiTitle = "CMS Simulation, L = 10 fb^{  -1}, #sqrt{s} = 13 TeV";
				
				TFile *pre = new TFile("Prediction.root","OPEN");
				TFile *exp = new TFile("Expectation.root","OPEN");
				
				TFile *outPutFile = new TFile("Compare.root","RECREATE");
				outPutFile->cd();
				std::cout<<"OutputFileCreated\n";
								
				DrawComparison((TH1F*)pre->Get("HT"), (TH1F*)exp->Get("HT"), (TString)"Title", LumiTitle, "H_{T} [GeV]", yTitle,isData)->Write();
				DrawComparison((TH1F*)pre->Get("MHT"), (TH1F*)exp->Get("MHT"), (TString)"Title", LumiTitle, "#slash{H}_{T} [GeV]", yTitle,isData)->Write();
				DrawComparison((TH1F*)pre->Get("NJets"), (TH1F*)exp->Get("NJets"), (TString)"Title", LumiTitle, "N_{Jets}", yTitle,isData)->Write();
				DrawComparison((TH1F*)pre->Get("BTags"), (TH1F*)exp->Get("BTags"), (TString)"Title", LumiTitle, "B_{Tags}", yTitle,isData)->Write();
				DrawComparison((TH1F*)pre->Get("Bin"), (TH1F*)exp->Get("Bin"), (TString)"Title", LumiTitle, "Search Region", yTitle,isData)->Write();
				
				
				outPutFile->Close();
}
