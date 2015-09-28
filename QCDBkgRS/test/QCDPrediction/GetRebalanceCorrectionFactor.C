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

//#include "/afs/desy.de/user/c/csander/xxl-af-cms/CMSSW_7_2_3_patch1/src/macros/PlottingUtils.C"

using namespace std;

void GetRebalanceCorrectionFactor()
{
   //gROOT->ProcessLine(".L /afs/desy.de/user/c/csander/xxl-af-cms/CMSSW_7_4_6_patch1/src/macros/PlottingUtils.C+");

   //SetPlotStyle();

   // --- define output file for ps ---//
   TString pt = "pt10";
   double x_min = 10.;

   TString outfile = "RebalanceCorrectionFactor/RebalanceCorrectionFactors_madgraph_spring15_withoutPUReweighting_" + pt;

   TH2F* RebCorrection_vsReco = new TH2F("RebCorrection_vsReco", "Jet pt", 1000, 0., 1000., 100, 0., 3.);
   TH2F* RebCorrection_vsReco_b = new TH2F("RebCorrection_vsReco_b", "Jet pt", 1000, 0., 1000., 100, 0., 3.);
   
   string root_file;

   // madgraph
   ifstream myfile ("filelists_phys14/filelist.txt");
   if (myfile.is_open()) {
      while( myfile.good() ) {
         getline (myfile,root_file);
         cout << root_file << endl;
                     
         TH2F* RebCorrection_vsReco_temp; 
         TH2F* RebCorrection_vsReco_b_temp;
         
         TString path = root_file;

         TFile* input_file = TFile::Open(path, "READ");
         input_file->cd("QCDfromSmearing");

         gDirectory->GetObject("RebCorrection_vsReco;1", RebCorrection_vsReco_temp);
         RebCorrection_vsReco->Add(RebCorrection_vsReco_temp);
         gDirectory->GetObject("RebCorrection_vsReco_b;1", RebCorrection_vsReco_b_temp);
         RebCorrection_vsReco_b->Add(RebCorrection_vsReco_b_temp);

         input_file->Close();
         
      }
      myfile.close();
   }

   // ---------------------------------------------------- //
   // madgraph
   TH1F* correction_vsReco = new TH1F();
   correction_vsReco = (TH1F*) RebCorrection_vsReco->ProjectionX();
   correction_vsReco->Reset();
   for (int i = 0; i <= RebCorrection_vsReco->GetXaxis()->GetNbins(); ++i) {
      TH1F h = *((TH1F*) RebCorrection_vsReco->ProjectionY("py", i, i));
            
      double mean = h.GetMean();
      double error = h.GetMeanError();

      cout << "i: " << i << " " << "mean: " << mean << " " << "error: " << error << endl;
            
      correction_vsReco->SetBinContent(i, mean);
      correction_vsReco->SetBinError(i, error);
   }

   TH1F* correction_vsReco_b = new TH1F();
   correction_vsReco_b = (TH1F*) RebCorrection_vsReco_b->ProjectionX();
   correction_vsReco_b->Reset();
   for (int i = 0; i <= RebCorrection_vsReco_b->GetXaxis()->GetNbins(); ++i) {
      TH1F h = *((TH1F*) RebCorrection_vsReco_b->ProjectionY("py", i, i));
      
      double mean = h.GetMean();
      double error = h.GetMeanError();
      
      cout << "i: " << i << " " << "mean: " << mean << " " << "error: " << error << endl;
      
      correction_vsReco_b->SetBinContent(i, mean);
      correction_vsReco_b->SetBinError(i, error);
   }

   // ---------------------------------------------------- //

   TCanvas *c = new TCanvas("c", "", 800, 800);
   correction_vsReco->SetMinimum(0.8);
   correction_vsReco->SetMaximum(1.2);
   correction_vsReco->SetAxisRange(x_min, correction_vsReco->GetXaxis()->GetXmax());
   correction_vsReco->SetXTitle("reco jet p_{T} [GeV]");
   correction_vsReco->SetYTitle("reb jet p_{T} / gen jet p_{T} ");
   correction_vsReco->Draw();

   c->Print(outfile + "_vsReco.ps");  
   c->Print(outfile + "_vsReco.png"); 

   TCanvas *d = new TCanvas("d", "", 800, 800);
   correction_vsReco_b->SetMinimum(0.8);
   correction_vsReco_b->SetMaximum(1.2);
   correction_vsReco_b->SetAxisRange(x_min, correction_vsReco->GetXaxis()->GetXmax());
   correction_vsReco_b->SetXTitle("reco b-jet p_{T} [GeV]");
   correction_vsReco_b->SetYTitle("reb b-jet p_{T} / gen jet p_{T} ");
   correction_vsReco_b->Draw();
   
   d->Print(outfile + "_vsReco_b.ps");
   d->Print(outfile + "_vsReco_b.png");
   // ---------------------------------------------------- //

   TFile* RebalanceCorrection = new TFile(outfile + ".root", "RECREATE");
   
   correction_vsReco->Write();
   correction_vsReco_b->Write();

   RebCorrection_vsReco->Write();
   RebCorrection_vsReco_b->Write();

   RebalanceCorrection->Write();
}


