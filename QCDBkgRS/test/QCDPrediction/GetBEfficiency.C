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

void GetBEfficiency()
{
   //gROOT->ProcessLine(".L /afs/desy.de/user/c/csander/xxl-af-cms/CMSSW_7_4_6_patch1/src/macros/PlottingUtils.C+");

   //SetPlotStyle();

   // --- define output file for ps ---//
   double x_min = 10.;

   TString outfile = "BTagEfficiency/BTagEfficiencies_Spring15MadGraph";

   // The number of true b's as a function of pt and eta
   std::vector<TH1F*> BTrue_vs_RecoPt_Eta;
   // The number of true b-tagged b's as a function of reco pt and eta
   std::vector<TH1F*> BTrue_BTag_vs_RecoPt_Eta;
   // The number of true b's that are b-tagged at reco level as a
   // function of reco pt and eta
   std::vector<TH1F*> BTagEfficiency_vs_RecoPt_Eta;

   // 12 = the number Eta bins -- Change if this changes!
   BTrue_vs_RecoPt_Eta.resize(12);
   BTrue_BTag_vs_RecoPt_Eta.resize(12);
   BTagEfficiency_vs_RecoPt_Eta.resize(12);

   for(int e_eta=0; e_eta<12; ++e_eta){
      char hname[100];
      // Book histograms b-tag efficiencies
      sprintf(hname, "BTrue_vs_RecoPt_Eta%i", e_eta);
      BTrue_vs_RecoPt_Eta.at(e_eta) = new TH1F(hname, hname, 200, 0., 1000.);
      sprintf(hname, "BTrue_BTag_vs_RecoPt_Eta%i", e_eta);
      BTrue_BTag_vs_RecoPt_Eta.at(e_eta) = new TH1F(hname, hname, 200, 0., 1000.);
      sprintf(hname, "BTagEfficiency_vs_RecoPt_Eta%i", e_eta);
      BTagEfficiency_vs_RecoPt_Eta.at(e_eta) = new TH1F(hname, hname, 200, 0., 1000.);
   }


   
   string root_file;

   // madgraph
   ifstream myfile ("filelists_phys14/filelist.txt");
   if (myfile.is_open()) {
      while( myfile.good() ) {
         getline (myfile,root_file);
         cout << root_file << endl;
         if(root_file == ""){
             continue;
         }
                     
         TH1F* BTrue_vs_RecoPt_temp; 
         TH1F* BTrue_BTag_vs_RecoPt_temp; 
         
         TString path = root_file;

         TFile* input_file = TFile::Open(path, "READ");


         for(int e_eta=0; e_eta<12; ++e_eta){
            char hname[100];
            // Book histograms with numbers of b-tagged events
            sprintf(hname, "h_trueb_RecoPt_Eta%i;1", e_eta);
            gDirectory->GetObject(hname, BTrue_vs_RecoPt_temp);
            BTrue_vs_RecoPt_Eta.at(e_eta)->Add(BTrue_vs_RecoPt_temp);
            cout << "Added trueb Histograms" << e_eta << endl;

            sprintf(hname, "h_trueb_btag_RecoPt_Eta%i;1", e_eta);
            gDirectory->GetObject(hname, BTrue_BTag_vs_RecoPt_temp);
            BTrue_BTag_vs_RecoPt_Eta.at(e_eta)->Add(BTrue_BTag_vs_RecoPt_temp);
            cout << "Added trueb and btag Histograms" << e_eta << endl;

         }

         cout << "Added Histograms" << endl;


         input_file->Close();
         
      }
      myfile.close();
   }


   cout << "Will divide Histograms" << endl;
   for(int e_eta=0; e_eta<12; ++e_eta){
        BTagEfficiency_vs_RecoPt_Eta.at(e_eta) = new TH1F(*BTrue_BTag_vs_RecoPt_Eta.at(e_eta));
        cout << "Created BTagEfficiency Vector" << e_eta << endl;
        // Divide # btags at reco level by # true b's at generator level:
        BTagEfficiency_vs_RecoPt_Eta.at(e_eta)->Divide(BTrue_vs_RecoPt_Eta.at(e_eta));
        cout << "Divided BTagEfficiency Vector" << e_eta << endl;
        // Can you rename the histogram BTagEfficiency?
        // It now carries the name BTrue_BTag_vs_RecoPt_... instead,
        // since you copied it from that one.
   }
         cout << "Divided Histograms" << endl;

      // ---------------------------------------------------- //

   TFile* BTagEfficiencies = new TFile(outfile + ".root", "RECREATE");

   for(int e_eta=0; e_eta<12; ++e_eta){
       BTagEfficiency_vs_RecoPt_Eta.at(e_eta)->Write();
   }


   BTagEfficiencies->Write();
}


