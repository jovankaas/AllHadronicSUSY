#include <TROOT.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TChain.h>
#include <TPad.h>
#include <TStyle.h>
#include <TFile.h>
#include <TPostScript.h>
#include <TLegend.h>
#include <TMath.h>
#include <TString.h>
#include <TArrayF.h>
#include <TLine.h>
#include <TPaveText.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

void DeltaPhiPlotter()
{
   
   gROOT->SetStyle("Plain");
   //gStyle->SetPalette(51, 0);
   
   gStyle->SetHatchesLineWidth(2.0);
   
   // For the canvas:
   gStyle->SetCanvasColor(0);
   gStyle->SetCanvasBorderMode(0);
   
   // For the Pad:
   gStyle->SetPadColor(0);
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   gStyle->SetPadBorderSize(0);
   gStyle->SetPadBorderMode(0);
   
   // For the frame:
   gStyle->SetFrameBorderMode(0);
   
   // For the histo:
   gStyle->SetMarkerSize(0.9);
   gStyle->SetMarkerStyle(20);
   gStyle->SetMarkerColor(1);
   
   // For the statistics box:
   gStyle->SetOptStat(0);
   //gStyle->SetOptFit(1011);
   
   // Margins:
   gStyle->SetPadBottomMargin(0.15);
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
   // gStyle->SetNdivisions(505, "X");
   //gStyle->SetNdivisions(510, "Y");
   gStyle->SetNdivisions(505, "Y");
   gStyle->SetTickLength(0.03, "XY");
   
   // For the axis titles:
   gStyle->SetTitleOffset(1.35, "X");
   gStyle->SetTitleOffset(1.52, "Y");
   gStyle->SetTitleOffset(0.5, "Z");
   //gStyle->SetTitleSize(0.045, "XY");
   gStyle->SetTitleSize(0.05, "XY");
   gStyle->SetTitleFont(42, "XYZ");
   
   // For the axis labels:
   //gStyle->SetLabelSize(0.035, "XYZ");
   gStyle->SetLabelSize(0.04, "XYZ");
   gStyle->SetLabelOffset(0.01, "XYZ");
   gStyle->SetLabelFont(42, "XYZ");
   
   // For the legend
   gStyle->SetLegendBorderSize(0);
   
   gROOT->ForceStyle();

   
   TH2F* hh_DeltaPhiRecoGenJet1_GenHT_LowMHT = new TH2F("h_DeltaPhiRecoGenJet1_GenHT_LowMHT", "h_DeltaPhiRecoGenJet1_GenHT_LowMHT", 100, 0, 0.5, 80, 0, 4000 );
   TH2F* hh_DeltaPhiRecoGenJet2_GenHT_LowMHT = new TH2F("h_DeltaPhiRecoGenJet2_GenHT_LowMHT", "h_DeltaPhiRecoGenJet2_GenHT_LowMHT", 100, 0, 0.5, 80, 0, 4000 );
   TH2F* hh_DeltaPhiRecoGenJet3_GenHT_LowMHT = new TH2F("h_DeltaPhiRecoGenJet3_GenHT_LowMHT", "h_DeltaPhiRecoGenJet3_GenHT_LowMHT", 100, 0, 0.5, 80, 0, 4000 );
   
   TH2F* hh_DeltaPhiRecoGenJet1_GenHT_HighMHT = new TH2F("h_DeltaPhiRecoGenJet1_GenHT_HighMHT", "h_DeltaPhiRecoGenJet1_GenHT_HighMHT", 100, 0, 0.5, 80, 0, 4000 );
   TH2F* hh_DeltaPhiRecoGenJet2_GenHT_HighMHT = new TH2F("h_DeltaPhiRecoGenJet2_GenHT_HighMHT", "h_DeltaPhiRecoGenJet2_GenHT_HighMHT", 100, 0, 0.5, 80, 0, 4000 );
   TH2F* hh_DeltaPhiRecoGenJet3_GenHT_HighMHT = new TH2F("h_DeltaPhiRecoGenJet3_GenHT_HighMHT", "h_DeltaPhiRecoGenJet3_GenHT_HighMHT", 100, 0, 0.5, 80, 0, 4000 );
   
   TH2F* hh_AddRelJetActivity07GenJet1_GenHT_LowMHT = new TH2F("h_AddRelJetActivity07GenJet1_GenHT_LowMHT", "h_AddRelJetActivity07GenJet1_GenHT_LowMHT", 100, 0, 1., 80, 0, 4000 );
   TH2F* hh_AddRelJetActivity07GenJet2_GenHT_LowMHT = new TH2F("h_AddRelJetActivity07GenJet2_GenHT_LowMHT", "h_AddRelJetActivity07GenJet2_GenHT_LowMHT", 100, 0, 1., 80, 0, 4000 );
   TH2F* hh_AddRelJetActivity07GenJet3_GenHT_LowMHT = new TH2F("h_AddRelJetActivity07GenJet3_GenHT_LowMHT", "h_AddRelJetActivity07GenJet3_GenHT_LowMHT", 100, 0, 1., 80, 0, 4000 );
   
   TH2F* hh_AddRelJetActivity07GenJet1_GenHT_HighMHT = new TH2F("h_AddRelJetActivity07GenJet1_GenHT_HighMHT", "h_AddRelJetActivity07GenJet1_GenHT_HighMHT", 100, 0, 1., 80, 0, 4000 );
   TH2F* hh_AddRelJetActivity07GenJet2_GenHT_HighMHT = new TH2F("h_AddRelJetActivity07GenJet2_GenHT_HighMHT", "h_AddRelJetActivity07GenJet2_GenHT_HighMHT", 100, 0, 1., 80, 0, 4000 );
   TH2F* hh_AddRelJetActivity07GenJet3_GenHT_HighMHT = new TH2F("h_AddRelJetActivity07GenJet3_GenHT_HighMHT", "h_AddRelJetActivity07GenJet3_GenHT_HighMHT", 100, 0, 1., 80, 0, 4000 );

   string root_file;

   ifstream myfile ("filelists_phys14/filelist.txt");
   if (myfile.is_open()) {
      while( myfile.good() ) {
         getline (myfile,root_file);
         cout << root_file << endl;
                     
         TH2F* h_DeltaPhiRecoGenJet1_GenHT_LowMHT_temp;
         TH2F* h_DeltaPhiRecoGenJet2_GenHT_LowMHT_temp;
         TH2F* h_DeltaPhiRecoGenJet3_GenHT_LowMHT_temp;
         TH2F* h_DeltaPhiRecoGenJet1_GenHT_HighMHT_temp;
         TH2F* h_DeltaPhiRecoGenJet2_GenHT_HighMHT_temp;
         TH2F* h_DeltaPhiRecoGenJet3_GenHT_HighMHT_temp;

         TH2F* h_AddRelJetActivity07GenJet1_GenHT_LowMHT_temp;
         TH2F* h_AddRelJetActivity07GenJet2_GenHT_LowMHT_temp;
         TH2F* h_AddRelJetActivity07GenJet3_GenHT_LowMHT_temp;
         TH2F* h_AddRelJetActivity07GenJet1_GenHT_HighMHT_temp;
         TH2F* h_AddRelJetActivity07GenJet2_GenHT_HighMHT_temp;
         TH2F* h_AddRelJetActivity07GenJet3_GenHT_HighMHT_temp;
         
         TString path = root_file;

         TFile* input_file = TFile::Open(path, "READ");
         input_file->cd("QCDfromSmearing");

         gDirectory->GetObject("h_DeltaPhiRecoGenJet1_GenHT_LowMHT;1", h_DeltaPhiRecoGenJet1_GenHT_LowMHT_temp);
         hh_DeltaPhiRecoGenJet1_GenHT_LowMHT->Add(h_DeltaPhiRecoGenJet1_GenHT_LowMHT_temp);
         gDirectory->GetObject("h_DeltaPhiRecoGenJet2_GenHT_LowMHT;1", h_DeltaPhiRecoGenJet2_GenHT_LowMHT_temp);
         hh_DeltaPhiRecoGenJet2_GenHT_LowMHT->Add(h_DeltaPhiRecoGenJet2_GenHT_LowMHT_temp);
         gDirectory->GetObject("h_DeltaPhiRecoGenJet3_GenHT_LowMHT;1", h_DeltaPhiRecoGenJet3_GenHT_LowMHT_temp);
         hh_DeltaPhiRecoGenJet3_GenHT_LowMHT->Add(h_DeltaPhiRecoGenJet3_GenHT_LowMHT_temp);

         gDirectory->GetObject("h_DeltaPhiRecoGenJet1_GenHT_HighMHT;1", h_DeltaPhiRecoGenJet1_GenHT_HighMHT_temp);
         hh_DeltaPhiRecoGenJet1_GenHT_HighMHT->Add(h_DeltaPhiRecoGenJet1_GenHT_HighMHT_temp);
         gDirectory->GetObject("h_DeltaPhiRecoGenJet2_GenHT_HighMHT;1", h_DeltaPhiRecoGenJet2_GenHT_HighMHT_temp);
         hh_DeltaPhiRecoGenJet2_GenHT_HighMHT->Add(h_DeltaPhiRecoGenJet2_GenHT_HighMHT_temp);
         gDirectory->GetObject("h_DeltaPhiRecoGenJet3_GenHT_HighMHT;1", h_DeltaPhiRecoGenJet3_GenHT_HighMHT_temp);
         hh_DeltaPhiRecoGenJet3_GenHT_HighMHT->Add(h_DeltaPhiRecoGenJet3_GenHT_HighMHT_temp);

         gDirectory->GetObject("h_AddRelJetActivity07GenJet1_GenHT_LowMHT;1", h_AddRelJetActivity07GenJet1_GenHT_LowMHT_temp);
         hh_AddRelJetActivity07GenJet1_GenHT_LowMHT->Add(h_AddRelJetActivity07GenJet1_GenHT_LowMHT_temp);
         gDirectory->GetObject("h_AddRelJetActivity07GenJet2_GenHT_LowMHT;1", h_AddRelJetActivity07GenJet2_GenHT_LowMHT_temp);
         hh_AddRelJetActivity07GenJet2_GenHT_LowMHT->Add(h_AddRelJetActivity07GenJet2_GenHT_LowMHT_temp);
         gDirectory->GetObject("h_AddRelJetActivity07GenJet3_GenHT_LowMHT;1", h_AddRelJetActivity07GenJet3_GenHT_LowMHT_temp);
         hh_AddRelJetActivity07GenJet3_GenHT_LowMHT->Add(h_AddRelJetActivity07GenJet3_GenHT_LowMHT_temp);
         
         gDirectory->GetObject("h_AddRelJetActivity07GenJet1_GenHT_HighMHT;1", h_AddRelJetActivity07GenJet1_GenHT_HighMHT_temp);
         hh_AddRelJetActivity07GenJet1_GenHT_HighMHT->Add(h_AddRelJetActivity07GenJet1_GenHT_HighMHT_temp);
         gDirectory->GetObject("h_AddRelJetActivity07GenJet2_GenHT_HighMHT;1", h_AddRelJetActivity07GenJet2_GenHT_HighMHT_temp);
         hh_AddRelJetActivity07GenJet2_GenHT_HighMHT->Add(h_AddRelJetActivity07GenJet2_GenHT_HighMHT_temp);
         gDirectory->GetObject("h_AddRelJetActivity07GenJet3_GenHT_HighMHT;1", h_AddRelJetActivity07GenJet3_GenHT_HighMHT_temp);
         hh_AddRelJetActivity07GenJet3_GenHT_HighMHT->Add(h_AddRelJetActivity07GenJet3_GenHT_HighMHT_temp);

         input_file->Close();
         
      }
      myfile.close();
   }

   // ---------------------------------------------------- //
   TH1F* h_DeltaPhiRecoGenJet1_GenHT_LowMHT_proj = new TH1F();
   h_DeltaPhiRecoGenJet1_GenHT_LowMHT_proj = (TH1F*) hh_DeltaPhiRecoGenJet1_GenHT_LowMHT->ProjectionY();
   TH1F* h_DeltaPhiRecoGenJet1_GenHT_HighMHT_proj = new TH1F();
   h_DeltaPhiRecoGenJet1_GenHT_HighMHT_proj = (TH1F*) hh_DeltaPhiRecoGenJet1_GenHT_HighMHT->ProjectionY();

   TH1F* h_DeltaPhiRecoGenJet2_GenHT_LowMHT_proj = new TH1F();
   h_DeltaPhiRecoGenJet2_GenHT_LowMHT_proj = (TH1F*) hh_DeltaPhiRecoGenJet2_GenHT_LowMHT->ProjectionY();
   TH1F* h_DeltaPhiRecoGenJet2_GenHT_HighMHT_proj = new TH1F();
   h_DeltaPhiRecoGenJet2_GenHT_HighMHT_proj = (TH1F*) hh_DeltaPhiRecoGenJet2_GenHT_HighMHT->ProjectionY();

   TH1F* h_DeltaPhiRecoGenJet3_GenHT_LowMHT_proj = new TH1F();
   h_DeltaPhiRecoGenJet3_GenHT_LowMHT_proj = (TH1F*) hh_DeltaPhiRecoGenJet3_GenHT_LowMHT->ProjectionY();
   TH1F* h_DeltaPhiRecoGenJet3_GenHT_HighMHT_proj = new TH1F();
   h_DeltaPhiRecoGenJet3_GenHT_HighMHT_proj = (TH1F*) hh_DeltaPhiRecoGenJet3_GenHT_HighMHT->ProjectionY();

   TH1F* h_AddRelJetActivity07GenJet1_GenHT_LowMHT_proj = new TH1F();
   h_AddRelJetActivity07GenJet1_GenHT_LowMHT_proj = (TH1F*) hh_AddRelJetActivity07GenJet1_GenHT_LowMHT->ProjectionY();
   TH1F* h_AddRelJetActivity07GenJet1_GenHT_HighMHT_proj = new TH1F();
   h_AddRelJetActivity07GenJet1_GenHT_HighMHT_proj = (TH1F*) hh_AddRelJetActivity07GenJet1_GenHT_HighMHT->ProjectionY();
   
   TH1F* h_AddRelJetActivity07GenJet2_GenHT_LowMHT_proj = new TH1F();
   h_AddRelJetActivity07GenJet2_GenHT_LowMHT_proj = (TH1F*) hh_AddRelJetActivity07GenJet2_GenHT_LowMHT->ProjectionY();
   TH1F* h_AddRelJetActivity07GenJet2_GenHT_HighMHT_proj = new TH1F();
   h_AddRelJetActivity07GenJet2_GenHT_HighMHT_proj = (TH1F*) hh_AddRelJetActivity07GenJet2_GenHT_HighMHT->ProjectionY();
   
   TH1F* h_AddRelJetActivity07GenJet3_GenHT_LowMHT_proj = new TH1F();
   h_AddRelJetActivity07GenJet3_GenHT_LowMHT_proj = (TH1F*) hh_AddRelJetActivity07GenJet3_GenHT_LowMHT->ProjectionY();
   TH1F* h_AddRelJetActivity07GenJet3_GenHT_HighMHT_proj = new TH1F();
   h_AddRelJetActivity07GenJet3_GenHT_HighMHT_proj = (TH1F*) hh_AddRelJetActivity07GenJet3_GenHT_HighMHT->ProjectionY();

   for (int i = 0; i <= hh_DeltaPhiRecoGenJet1_GenHT_LowMHT->GetYaxis()->GetNbins(); ++i) {
      
      TH1F h1_low = *((TH1F*) hh_DeltaPhiRecoGenJet1_GenHT_LowMHT->ProjectionX("px", i, i));
      double mean_low1 = h1_low.GetMean();
      double error_low1 = h1_low.GetMeanError();
      h_DeltaPhiRecoGenJet1_GenHT_LowMHT_proj->SetBinContent(i, mean_low1);
      h_DeltaPhiRecoGenJet1_GenHT_LowMHT_proj->SetBinError(i, error_low1);
      
      TH1F h1_high = *((TH1F*) hh_DeltaPhiRecoGenJet1_GenHT_HighMHT->ProjectionX("px", i, i));
      double mean_high1 = h1_high.GetMean();
      double error_high1 = h1_high.GetMeanError();
      h_DeltaPhiRecoGenJet1_GenHT_HighMHT_proj->SetBinContent(i, mean_high1);
      h_DeltaPhiRecoGenJet1_GenHT_HighMHT_proj->SetBinError(i, error_high1);

      TH1F h2_low = *((TH1F*) hh_DeltaPhiRecoGenJet2_GenHT_LowMHT->ProjectionX("px", i, i));
      double mean_low2 = h2_low.GetMean();
      double error_low2 = h2_low.GetMeanError();
      h_DeltaPhiRecoGenJet2_GenHT_LowMHT_proj->SetBinContent(i, mean_low2);
      h_DeltaPhiRecoGenJet2_GenHT_LowMHT_proj->SetBinError(i, error_low2);
      
      TH1F h2_high = *((TH1F*) hh_DeltaPhiRecoGenJet2_GenHT_HighMHT->ProjectionX("px", i, i));
      double mean_high2 = h2_high.GetMean();
      double error_high2 = h2_high.GetMeanError();
      h_DeltaPhiRecoGenJet2_GenHT_HighMHT_proj->SetBinContent(i, mean_high2);
      h_DeltaPhiRecoGenJet2_GenHT_HighMHT_proj->SetBinError(i, error_high2);
      
      TH1F h3_low = *((TH1F*) hh_DeltaPhiRecoGenJet3_GenHT_LowMHT->ProjectionX("px", i, i));
      double mean_low3 = h3_low.GetMean();
      double error_low3 = h3_low.GetMeanError();
      h_DeltaPhiRecoGenJet3_GenHT_LowMHT_proj->SetBinContent(i, mean_low3);
      h_DeltaPhiRecoGenJet3_GenHT_LowMHT_proj->SetBinError(i, error_low3);
      
      TH1F h3_high = *((TH1F*) hh_DeltaPhiRecoGenJet3_GenHT_HighMHT->ProjectionX("px", i, i));
      double mean_high3 = h3_high.GetMean();
      double error_high3 = h3_high.GetMeanError();
      h_DeltaPhiRecoGenJet3_GenHT_HighMHT_proj->SetBinContent(i, mean_high3);
      h_DeltaPhiRecoGenJet3_GenHT_HighMHT_proj->SetBinError(i, error_high3);

   }

   for (int i = 0; i <= hh_AddRelJetActivity07GenJet1_GenHT_LowMHT->GetYaxis()->GetNbins(); ++i) {
      
      TH1F h1_low = *((TH1F*) hh_AddRelJetActivity07GenJet1_GenHT_LowMHT->ProjectionX("px", i, i));
      double mean_low1 = h1_low.GetMean();
      double error_low1 = h1_low.GetMeanError();
      h_AddRelJetActivity07GenJet1_GenHT_LowMHT_proj->SetBinContent(i, mean_low1);
      h_AddRelJetActivity07GenJet1_GenHT_LowMHT_proj->SetBinError(i, error_low1);
      
      TH1F h1_high = *((TH1F*) hh_AddRelJetActivity07GenJet1_GenHT_HighMHT->ProjectionX("px", i, i));
      double mean_high1 = h1_high.GetMean();
      double error_high1 = h1_high.GetMeanError();
      h_AddRelJetActivity07GenJet1_GenHT_HighMHT_proj->SetBinContent(i, mean_high1);
      h_AddRelJetActivity07GenJet1_GenHT_HighMHT_proj->SetBinError(i, error_high1);
      
      TH1F h2_low = *((TH1F*) hh_AddRelJetActivity07GenJet2_GenHT_LowMHT->ProjectionX("px", i, i));
      double mean_low2 = h2_low.GetMean();
      double error_low2 = h2_low.GetMeanError();
      h_AddRelJetActivity07GenJet2_GenHT_LowMHT_proj->SetBinContent(i, mean_low2);
      h_AddRelJetActivity07GenJet2_GenHT_LowMHT_proj->SetBinError(i, error_low2);
      
      TH1F h2_high = *((TH1F*) hh_AddRelJetActivity07GenJet2_GenHT_HighMHT->ProjectionX("px", i, i));
      double mean_high2 = h2_high.GetMean();
      double error_high2 = h2_high.GetMeanError();
      h_AddRelJetActivity07GenJet2_GenHT_HighMHT_proj->SetBinContent(i, mean_high2);
      h_AddRelJetActivity07GenJet2_GenHT_HighMHT_proj->SetBinError(i, error_high2);
      
      TH1F h3_low = *((TH1F*) hh_AddRelJetActivity07GenJet3_GenHT_LowMHT->ProjectionX("px", i, i));
      double mean_low3 = h3_low.GetMean();
      double error_low3 = h3_low.GetMeanError();
      h_AddRelJetActivity07GenJet3_GenHT_LowMHT_proj->SetBinContent(i, mean_low3);
      h_AddRelJetActivity07GenJet3_GenHT_LowMHT_proj->SetBinError(i, error_low3);
      
      TH1F h3_high = *((TH1F*) hh_AddRelJetActivity07GenJet3_GenHT_HighMHT->ProjectionX("px", i, i));
      double mean_high3 = h3_high.GetMean();
      double error_high3 = h3_high.GetMeanError();
      h_AddRelJetActivity07GenJet3_GenHT_HighMHT_proj->SetBinContent(i, mean_high3);
      h_AddRelJetActivity07GenJet3_GenHT_HighMHT_proj->SetBinError(i, error_high3);
      
   }

   // ---------------------------------------------------- //

   TCanvas *c1 = new TCanvas("c1", "", 800., 800.);

   h_DeltaPhiRecoGenJet1_GenHT_LowMHT_proj->SetXTitle("gen H_{T} [GeV]");
   h_DeltaPhiRecoGenJet1_GenHT_LowMHT_proj->SetYTitle("mean #Delta R(genJet1, recoJet)");
   h_DeltaPhiRecoGenJet1_GenHT_LowMHT_proj->SetMarkerColor(3);
   h_DeltaPhiRecoGenJet1_GenHT_LowMHT_proj->Draw();
   h_DeltaPhiRecoGenJet1_GenHT_HighMHT_proj->SetMarkerColor(2);
   h_DeltaPhiRecoGenJet1_GenHT_HighMHT_proj->Draw("same");
   
   TLegend* leg1 = new TLegend(0.44, 0.63, 0.91, 0.83);
   leg1->SetFillStyle(0);
   leg1->SetLineStyle(1);
   leg1->SetTextFont(42);
   // leg1->SetTextSize(0.04);
   leg1->SetTextSize(0.045);
   leg1->AddEntry(h_DeltaPhiRecoGenJet1_GenHT_LowMHT_proj, "MHT < 100 GeV", "lep");
   leg1->AddEntry(h_DeltaPhiRecoGenJet1_GenHT_HighMHT_proj, "MHT > 100 GeV", "lep");
   leg1->Draw("same");

   c1->Print("DeltaRRecoGenJet1_GenHT.pdf");

   TCanvas *c2 = new TCanvas("c2", "", 800., 800.);
   h_DeltaPhiRecoGenJet2_GenHT_LowMHT_proj->SetXTitle("gen H_{T} [GeV]");
   h_DeltaPhiRecoGenJet2_GenHT_LowMHT_proj->SetYTitle("mean #Delta R(genJet2, recoJet)");
   h_DeltaPhiRecoGenJet2_GenHT_LowMHT_proj->SetMarkerColor(3);
   h_DeltaPhiRecoGenJet2_GenHT_LowMHT_proj->Draw();
   h_DeltaPhiRecoGenJet2_GenHT_HighMHT_proj->SetMarkerColor(2);
   h_DeltaPhiRecoGenJet2_GenHT_HighMHT_proj->Draw("same");
   leg1->Draw("same");
   c2->Print("DeltaRRecoGenJet2_GenHT.pdf");

   TCanvas *c3 = new TCanvas("c3", "", 800., 800.);
   h_DeltaPhiRecoGenJet3_GenHT_LowMHT_proj->SetXTitle("gen H_{T} [GeV]");
   h_DeltaPhiRecoGenJet3_GenHT_LowMHT_proj->SetYTitle("mean #Delta R(genJet3, recoJet)");
   h_DeltaPhiRecoGenJet3_GenHT_LowMHT_proj->SetMarkerColor(3);
   h_DeltaPhiRecoGenJet3_GenHT_LowMHT_proj->Draw();
   h_DeltaPhiRecoGenJet3_GenHT_HighMHT_proj->SetMarkerColor(2);
   h_DeltaPhiRecoGenJet3_GenHT_HighMHT_proj->Draw("same");
   leg1->Draw("same");
   c3->Print("DeltaRRecoGenJet3_GenHT.pdf");

   TCanvas *d1 = new TCanvas("d1", "", 800., 800.);
   h_AddRelJetActivity07GenJet1_GenHT_LowMHT_proj->SetXTitle("gen H_{T} [GeV]");
   h_AddRelJetActivity07GenJet1_GenHT_LowMHT_proj->SetYTitle("Add. recoJet activity dR<0.7 (genJet1)");
   h_AddRelJetActivity07GenJet1_GenHT_LowMHT_proj->SetMarkerColor(3);
   h_AddRelJetActivity07GenJet1_GenHT_LowMHT_proj->Draw();
   h_AddRelJetActivity07GenJet1_GenHT_HighMHT_proj->SetMarkerColor(2);
   h_AddRelJetActivity07GenJet1_GenHT_HighMHT_proj->Draw("same");
   leg1->Draw("same");
   d1->Print("AddRelJetActivity07GenJet1_GenHT.pdf");
   
   TCanvas *d2 = new TCanvas("d2", "", 800., 800.);
   h_AddRelJetActivity07GenJet2_GenHT_LowMHT_proj->SetXTitle("gen H_{T} [GeV]");
   h_AddRelJetActivity07GenJet2_GenHT_LowMHT_proj->SetYTitle("Add. recoJet activity dR<0.7 (genJet2)");
   h_AddRelJetActivity07GenJet2_GenHT_LowMHT_proj->SetMarkerColor(3);
   h_AddRelJetActivity07GenJet2_GenHT_LowMHT_proj->Draw();
   h_AddRelJetActivity07GenJet2_GenHT_HighMHT_proj->SetMarkerColor(2);
   h_AddRelJetActivity07GenJet2_GenHT_HighMHT_proj->Draw("same");
   leg1->Draw("same");
   d2->Print("AddRelJetActivity07GenJet2_GenHT.pdf");
   
   TCanvas *d3 = new TCanvas("d3", "", 800., 800.);
   h_AddRelJetActivity07GenJet3_GenHT_LowMHT_proj->SetXTitle("gen H_{T} [GeV]");
   h_AddRelJetActivity07GenJet3_GenHT_LowMHT_proj->SetYTitle("Add. recoJet activity dR<0.7 (genJet3)");
   h_AddRelJetActivity07GenJet3_GenHT_LowMHT_proj->SetMarkerColor(3);
   h_AddRelJetActivity07GenJet3_GenHT_LowMHT_proj->Draw();
   h_AddRelJetActivity07GenJet3_GenHT_HighMHT_proj->SetMarkerColor(2);
   h_AddRelJetActivity07GenJet3_GenHT_HighMHT_proj->Draw("same");
   leg1->Draw("same");
   d3->Print("AddRelJetActivity07GenJet3_GenHT.pdf");
   
   // ---------------------------------------------------- //

}


