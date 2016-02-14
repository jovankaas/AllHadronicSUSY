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

void GetBMisTag()
{
   //gROOT->ProcessLine(".L /afs/desy.de/user/c/csander/xxl-af-cms/CMSSW_7_4_6_patch1/src/macros/PlottingUtils.C+");

   //SetPlotStyle();

   // --- define output file for ps ---//
   TString outfile = "BTagEfficiency/B_Mis_TagEfficiencies_Spring15MadGraph_DeadECALTP_genPTbins";

   // The number of true b's as a function of pt and eta
   std::vector<TH1F*> BTrue_vs_RecoPt_Eta;
   // The number of true b-tagged b's as a function of reco pt and eta
   std::vector<TH1F*> BTrue_BTag_vs_RecoPt_Eta;
   // The number of true b's that are b-tagged at reco level as a
   // function of reco pt and eta
   std::vector<TH1F*> BTagEfficiency_vs_RecoPt_Eta;


   // The number of non-b's as a function of pt and eta
   std::vector<TH1F*> no_BTrue_vs_RecoPt_Eta;
   // The number of non-btagged jets as a function of pt and eta
   std::vector<TH1F*> no_BTag_vs_RecoPt_Eta;
   // The number of btags as a function of pt and eta
   std::vector<TH1F*> BTag_vs_RecoPt_Eta;
   // The number of b-tagged non-b's as a function of reco pt and eta
   std::vector<TH1F*> no_BTrue_BTag_vs_RecoPt_Eta;
   // The number of non-b's that are b-tagged at reco level as a
   // function of reco pt and eta
   std::vector<TH1F*> BMisTagEfficiency_vs_RecoPt_Eta;

   // 12 = the number Eta bins -- Change if this changes!
   BTrue_vs_RecoPt_Eta.resize(12);
   no_BTrue_vs_RecoPt_Eta.resize(12);
   no_BTag_vs_RecoPt_Eta.resize(12);
   BTag_vs_RecoPt_Eta.resize(12);
   BTrue_BTag_vs_RecoPt_Eta.resize(12);
   no_BTrue_BTag_vs_RecoPt_Eta.resize(12);
   BTagEfficiency_vs_RecoPt_Eta.resize(12);
   BMisTagEfficiency_vs_RecoPt_Eta.resize(12);

   for(int e_eta=0; e_eta<12; ++e_eta){
      char hname[100];
      // Book histograms b-(mis)tag efficiencies
      sprintf(hname, "BTrue_vs_RecoPt_Eta%i", e_eta);
      BTrue_vs_RecoPt_Eta.at(e_eta) = new TH1F(hname, hname, 200, 0., 1000.);
      sprintf(hname, "BTrue_BTag_vs_RecoPt_Eta%i", e_eta);
      BTrue_BTag_vs_RecoPt_Eta.at(e_eta) = new TH1F(hname, hname, 200, 0., 1000.);
      sprintf(hname, "BTagEfficiency_vs_RecoPt_Eta%i", e_eta);
      BTagEfficiency_vs_RecoPt_Eta.at(e_eta) = new TH1F(hname, hname, 200, 0., 1000.);
      sprintf(hname, "BTag_vs_RecoPt_Eta%i", e_eta);
      BTag_vs_RecoPt_Eta.at(e_eta) = new TH1F(hname, hname, 200, 0., 1000.);
      sprintf(hname, "no_BTrue_vs_RecoPt_Eta%i", e_eta);
      no_BTrue_vs_RecoPt_Eta.at(e_eta) = new TH1F(hname, hname, 200, 0., 1000.);
      sprintf(hname, "no_BTag_vs_RecoPt_Eta%i", e_eta);
      no_BTag_vs_RecoPt_Eta.at(e_eta) = new TH1F(hname, hname, 200, 0., 1000.);
      sprintf(hname, "no_BTrue_BTag_vs_RecoPt_Eta%i", e_eta);
      no_BTrue_BTag_vs_RecoPt_Eta.at(e_eta) = new TH1F(hname, hname, 200, 0., 1000.);
      sprintf(hname, "BMisTagEfficiency_vs_RecoPt_Eta%i", e_eta);
      BMisTagEfficiency_vs_RecoPt_Eta.at(e_eta) = new TH1F(hname, hname, 200, 0., 1000.);
   }


   string root_file;

   // madgraph
   ifstream myfile ("filelists_phys14/filelist_mc_genpTbins_deadECAL_newsoftware.txt");
   //ifstream myfile ("filelists_phys14/test.txt");
   if (myfile.is_open()) {
      while( myfile.good() ) {
         getline (myfile,root_file);
         //cout << root_file << endl;
         if(root_file == ""){
             continue;
         }
                     
         TH1F* BTrue_vs_RecoPt_temp; 
         TH1F* no_BTrue_vs_RecoPt_temp; 
         TH1F* no_BTag_vs_RecoPt_temp; 
         TH1F* BTag_vs_RecoPt_temp; 
         TH1F* BTrue_BTag_vs_RecoPt_temp; 
         TH1F* no_BTrue_BTag_vs_RecoPt_temp; 
         
         TString path = root_file;

         TFile* input_file = TFile::Open(path, "READ");


         for(int e_eta=0; e_eta<12; ++e_eta){
            char hname[100];
            // Book histograms with numbers of btrue / b-tags
            sprintf(hname, "h_trueb_RecoPt_Eta%i;1", e_eta);
            gDirectory->GetObject(hname, BTrue_vs_RecoPt_temp);
            BTrue_vs_RecoPt_Eta.at(e_eta)->Add(BTrue_vs_RecoPt_temp);
            //cout << "Added trueb Histograms" << e_eta << endl;

            sprintf(hname, "h_no_trueb_RecoPt_Eta%i;1", e_eta);
            gDirectory->GetObject(hname, no_BTrue_vs_RecoPt_temp);
            no_BTrue_vs_RecoPt_Eta.at(e_eta)->Add(no_BTrue_vs_RecoPt_temp);
            //cout << "Added no trueb Histograms" << e_eta << endl;

            sprintf(hname, "h_no_btag_RecoPt_Eta%i;1", e_eta);
            gDirectory->GetObject(hname, no_BTag_vs_RecoPt_temp);
            no_BTag_vs_RecoPt_Eta.at(e_eta)->Add(no_BTag_vs_RecoPt_temp);
            //cout << "Added no btag Histograms" << e_eta << endl;

            sprintf(hname, "h_trueb_btag_RecoPt_Eta%i;1", e_eta);
            gDirectory->GetObject(hname, BTrue_BTag_vs_RecoPt_temp);
            BTrue_BTag_vs_RecoPt_Eta.at(e_eta)->Add(BTrue_BTag_vs_RecoPt_temp);
            //cout << "Added trueb and btag Histograms" << e_eta << endl;

            sprintf(hname, "h_btag_RecoPt_Eta%i;1", e_eta);
            gDirectory->GetObject(hname, BTag_vs_RecoPt_temp);
            BTag_vs_RecoPt_Eta.at(e_eta)->Add(BTag_vs_RecoPt_temp);
            //cout << "Added btag Histograms" << e_eta << endl;

            sprintf(hname, "h_no_trueb_btag_RecoPt_Eta%i;1", e_eta);
            gDirectory->GetObject(hname, no_BTrue_BTag_vs_RecoPt_temp);
            no_BTrue_BTag_vs_RecoPt_Eta.at(e_eta)->Add(no_BTrue_BTag_vs_RecoPt_temp);
            //cout << "Added non-b and btag Histograms" << e_eta << endl;
         }

         //cout << "Added Histograms" << endl;


         input_file->Close();
      }
      myfile.close();
   }


   //cout << "Will divide Histograms" << endl;
   for(int e_eta=0; e_eta<12; ++e_eta){
        BTagEfficiency_vs_RecoPt_Eta.at(e_eta) = new TH1F(*BTrue_BTag_vs_RecoPt_Eta.at(e_eta));
        //cout << "Created BTagEfficiency Vector" << e_eta << endl;
        // Divide # btags at reco level by # true b's at generator level:
        BTagEfficiency_vs_RecoPt_Eta.at(e_eta)->Divide(BTrue_vs_RecoPt_Eta.at(e_eta));
        //cout << "Divided BTagEfficiency Vector" << e_eta << endl;
        // Can you rename the histogram BTagEfficiency?
        // It now carries the name BTrue_BTag_vs_RecoPt_... instead,
        // since you copied it from that one.

        BMisTagEfficiency_vs_RecoPt_Eta.at(e_eta) = new TH1F(*no_BTrue_BTag_vs_RecoPt_Eta.at(e_eta));
        //cout << "Created BMisTagEfficiency Vector" << e_eta << endl;
        BMisTagEfficiency_vs_RecoPt_Eta.at(e_eta)->Divide(no_BTrue_vs_RecoPt_Eta.at(e_eta));
        //cout << "Divided BTagEfficiency Vector" << e_eta << endl;
   }
         //cout << "Divided Histograms" << endl;

      // ---------------------------------------------------- //

   TFile* BTagEfficiencies = new TFile(outfile + ".root", "RECREATE");

   for(int e_eta=0; e_eta<12; ++e_eta){
       BTrue_vs_RecoPt_Eta.at(e_eta)->Write();
       BTag_vs_RecoPt_Eta.at(e_eta)->Write();
       no_BTrue_vs_RecoPt_Eta.at(e_eta)->Write();
       no_BTag_vs_RecoPt_Eta.at(e_eta)->Write();
       //no_BTrue_BTag_vs_RecoPt_Eta.at(e_eta)->Write();
       BTagEfficiency_vs_RecoPt_Eta.at(e_eta)->Write();
       BMisTagEfficiency_vs_RecoPt_Eta.at(e_eta)->Write();
   }


   BTagEfficiencies->Write();
}


