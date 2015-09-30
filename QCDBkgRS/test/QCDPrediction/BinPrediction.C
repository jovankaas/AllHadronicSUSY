//
//  BinPrediction.c
//
//
//  Created by Christian Sander on 15/06/15.
//
//

#include "BinPrediction.h"
#include "RA2bBin.h"

#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TMath.h>
#include <TString.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TGraphAsymmErrors.h>
#include <iostream>
#include <vector>
#include <string>

using namespace std;

BinPrediction::BinPrediction(TChain& QCDPrediction, TChain& RA2PreSelection, TString& Uncertainty)
{
   gROOT->ProcessLine("#include <vector>");
   
   // ------------- define all histos needed -------//
   // set histogram attributes
   int Ntries = 100;
   
   //Define search bins
   //￼Njets bins: 4-6, 7-8, 9+
   //Nb bins: 0, 1, 2, 3+
   //MHT,HT = 200,500,500,800
   //MHT,HT = 200,500,800,1200
   //MHT,HT = 200,500,1200+
   //MHT,HT = 500,750,500,1200
   //MHT,HT = 500,750,1200+
   //MHT,HT = 750+,800+
   
   std::vector<RA2bBin*> SB;
   
   RA2bBin SB1("Njet456-Nb0-MHT1-HT1", Ntries, 4, 6, 0, 0, 200., 500., 500., 800., true);
   RA2bBin SB2("Njet456-Nb0-MHT1-HT2", Ntries, 4, 6, 0, 0, 200., 500., 800., 1200., true);
   RA2bBin SB3("Njet456-Nb0-MHT1-HT3", Ntries, 4, 6, 0, 0, 200., 500., 1200., 9999., true);
   RA2bBin SB4("Njet456-Nb0-MHT2-HT1", Ntries, 4, 6, 0, 0, 500., 750., 500., 1200., true);
   RA2bBin SB5("Njet456-Nb0-MHT2-HT23", Ntries, 4, 6, 0, 0, 500., 750., 1200., 9999., true);
   RA2bBin SB6("Njet456-Nb0-MHT2-HT123", Ntries, 4, 6, 0, 0, 750., 9999., 800., 9999., true);
   
   RA2bBin SB7("Njet456-Nb1-MHT1-HT1", Ntries, 4, 6, 1, 1, 200., 500., 500., 800., true);
   RA2bBin SB8("Njet456-Nb1-MHT1-HT2", Ntries, 4, 6, 1, 1, 200., 500., 800., 1200., true);
   RA2bBin SB9("Njet456-Nb1-MHT1-HT3", Ntries, 4, 6, 1, 1, 200., 500., 1200., 9999., true);
   RA2bBin SB10("Njet456-Nb1-MHT2-HT1", Ntries, 4, 6, 1, 1, 500., 750., 500., 1200., true);
   RA2bBin SB11("Njet456-Nb1-MHT2-HT23", Ntries, 4, 6, 1, 1, 500., 750., 1200., 9999., true);
   RA2bBin SB12("Njet456-Nb1-MHT2-HT123", Ntries, 4, 6, 1, 1, 750., 9999., 800., 9999., true);
   
   RA2bBin SB13("Njet456-Nb2-MHT1-HT1", Ntries, 4, 6, 2, 2, 200., 500., 500., 800., true);
   RA2bBin SB14("Njet456-Nb2-MHT1-HT2", Ntries, 4, 6, 2, 2, 200., 500., 800., 1200., true);
   RA2bBin SB15("Njet456-Nb2-MHT1-HT3", Ntries, 4, 6, 2, 2, 200., 500., 1200., 9999., true);
   RA2bBin SB16("Njet456-Nb2-MHT2-HT1", Ntries, 4, 6, 2, 2, 500., 750., 500., 1200., true);
   RA2bBin SB17("Njet456-Nb2-MHT2-HT23", Ntries, 4, 6, 2, 2, 500., 750., 1200., 9999., true);
   RA2bBin SB18("Njet456-Nb2-MHT2-HT123", Ntries, 4, 6, 2, 2, 750., 9999., 800., 9999., true);
   
   RA2bBin SB19("Njet456-Nb3+-MHT1-HT1", Ntries, 4, 6, 3, 9, 200., 500., 500., 800., true);
   RA2bBin SB20("Njet456-Nb3+-MHT1-HT2", Ntries, 4, 6, 3, 9, 200., 500., 800., 1200., true);
   RA2bBin SB21("Njet456-Nb3+-MHT1-HT3", Ntries, 4, 6, 3, 9, 200., 500., 1200., 9999., true);
   RA2bBin SB22("Njet456-Nb3+-MHT2-HT1", Ntries, 4, 6, 3, 9, 500., 750., 500., 1200., true);
   RA2bBin SB23("Njet456-Nb3+-MHT2-HT23", Ntries, 4, 6, 3, 9, 500., 750., 1200., 9999., true);
   RA2bBin SB24("Njet456-Nb3+-MHT2-HT123", Ntries, 4, 6, 3, 9, 750., 9999., 800., 9999., true);
   
   RA2bBin SB25("Njet78-Nb0-MHT1-HT1", Ntries, 7, 8, 0, 0, 200., 500., 500., 800., true);
   RA2bBin SB26("Njet78-Nb0-MHT1-HT2", Ntries, 7, 8, 0, 0, 200., 500., 800., 1200., true);
   RA2bBin SB27("Njet78-Nb0-MHT1-HT3", Ntries, 7, 8, 0, 0, 200., 500., 1200., 9999., true);
   RA2bBin SB28("Njet78-Nb0-MHT2-HT1", Ntries, 7, 8, 0, 0, 500., 750., 500., 1200., true);
   RA2bBin SB29("Njet78-Nb0-MHT2-HT23", Ntries, 7, 8, 0, 0, 500., 750., 1200., 9999., true);
   RA2bBin SB30("Njet78-Nb0-MHT2-HT123", Ntries, 7, 8, 0, 0, 750., 9999., 800., 9999., true);
   
   RA2bBin SB31("Njet78-Nb1-MHT1-HT1", Ntries, 7, 8, 1, 1, 200., 500., 500., 800., true);
   RA2bBin SB32("Njet78-Nb1-MHT1-HT2", Ntries, 7, 8, 1, 1, 200., 500., 800., 1200., true);
   RA2bBin SB33("Njet78-Nb1-MHT1-HT3", Ntries, 7, 8, 1, 1, 200., 500., 1200., 9999., true);
   RA2bBin SB34("Njet78-Nb1-MHT2-HT1", Ntries, 7, 8, 1, 1, 500., 750., 500., 1200., true);
   RA2bBin SB35("Njet78-Nb1-MHT2-HT23", Ntries, 7, 8, 1, 1, 500., 750., 1200., 9999., true);
   RA2bBin SB36("Njet78-Nb1-MHT2-HT123", Ntries, 7, 8, 1, 1, 750., 9999., 800., 9999., true);
   
   RA2bBin SB37("Njet78-Nb2-MHT1-HT1", Ntries, 7, 8, 2, 2, 200., 500., 500., 800., true);
   RA2bBin SB38("Njet78-Nb2-MHT1-HT2", Ntries, 7, 8, 2, 2, 200., 500., 800., 1200., true);
   RA2bBin SB39("Njet78-Nb2-MHT1-HT3", Ntries, 7, 8, 2, 2, 200., 500., 1200., 9999., true);
   RA2bBin SB40("Njet78-Nb2-MHT2-HT1", Ntries, 7, 8, 2, 2, 500., 750., 500., 1200., true);
   RA2bBin SB41("Njet78-Nb2-MHT2-HT23", Ntries, 7, 8, 2, 2, 500., 750., 1200., 9999., true);
   RA2bBin SB42("Njet78-Nb2-MHT2-HT123", Ntries, 7, 8, 2, 2, 750., 9999., 800., 9999., true);
   
   RA2bBin SB43("Njet78-Nb3+-MHT1-HT1", Ntries, 7, 8, 3, 9, 200., 500., 500., 800., true);
   RA2bBin SB44("Njet78-Nb3+-MHT1-HT2", Ntries, 7, 8, 3, 9, 200., 500., 800., 1200., true);
   RA2bBin SB45("Njet78-Nb3+-MHT1-HT3", Ntries, 7, 8, 3, 9, 200., 500., 1200., 9999., true);
   RA2bBin SB46("Njet78-Nb3+-MHT2-HT1", Ntries, 7, 8, 3, 9, 500., 750., 500., 1200., true);
   RA2bBin SB47("Njet78-Nb3+-MHT2-HT23", Ntries, 7, 8, 3, 9, 500., 750., 1200., 9999., true);
   RA2bBin SB48("Njet78-Nb3+-MHT2-HT123", Ntries, 7, 8, 3, 9, 750., 9999., 800., 9999., true);
   
   RA2bBin SB49("Njet9+-Nb0-MHT1-HT1", Ntries, 9, 99, 0, 0, 200., 500., 500., 800., true);
   RA2bBin SB50("Njet9+-Nb0-MHT1-HT2", Ntries, 9, 99, 0, 0, 200., 500., 800., 1200., true);
   RA2bBin SB51("Njet9+-Nb0-MHT1-HT3", Ntries, 9, 99, 0, 0, 200., 500., 1200., 9999., true);
   RA2bBin SB52("Njet9+-Nb0-MHT2-HT1", Ntries, 9, 99, 0, 0, 500., 750., 500., 1200., true);
   RA2bBin SB53("Njet9+-Nb0-MHT2-HT23", Ntries, 9, 99, 0, 0, 500., 750., 1200., 9999., true);
   RA2bBin SB54("Njet9+-Nb0-MHT2-HT123", Ntries, 9, 99, 0, 0, 750., 9999., 800., 9999., true);
   
   RA2bBin SB55("Njet9+-Nb1-MHT1-HT1", Ntries, 9, 99, 1, 1, 200., 500., 500., 800., true);
   RA2bBin SB56("Njet9+-Nb1-MHT1-HT2", Ntries, 9, 99, 1, 1, 200., 500., 800., 1200., true);
   RA2bBin SB57("Njet9+-Nb1-MHT1-HT3", Ntries, 9, 99, 1, 1, 200., 500., 1200., 9999., true);
   RA2bBin SB58("Njet9+-Nb1-MHT2-HT1", Ntries, 9, 99, 1, 1, 500., 750., 500., 1200., true);
   RA2bBin SB59("Njet9+-Nb1-MHT2-HT23", Ntries, 9, 99, 1, 1, 500., 750., 1200., 9999., true);
   RA2bBin SB60("Njet9+-Nb1-MHT2-HT123", Ntries, 9, 99, 1, 1, 750., 9999., 800., 9999., true);
   
   RA2bBin SB61("Njet9+-Nb2-MHT1-HT1", Ntries, 9, 99, 2, 2, 200., 500., 500., 800., true);
   RA2bBin SB62("Njet9+-Nb2-MHT1-HT2", Ntries, 9, 99, 2, 2, 200., 500., 800., 1200., true);
   RA2bBin SB63("Njet9+-Nb2-MHT1-HT3", Ntries, 9, 99, 2, 2, 200., 500., 1200., 9999., true);
   RA2bBin SB64("Njet9+-Nb2-MHT2-HT1", Ntries, 9, 99, 2, 2, 500., 750., 500., 1200., true);
   RA2bBin SB65("Njet9+-Nb2-MHT2-HT23", Ntries, 9, 99, 2, 2, 500., 750., 1200., 9999., true);
   RA2bBin SB66("Njet9+-Nb2-MHT2-HT123", Ntries, 9, 99, 2, 2, 750., 9999., 800., 9999., true);
   
   RA2bBin SB67("Njet9+-Nb3+-MHT1-HT1", Ntries, 9, 99, 3, 9, 200., 500., 500., 800., true);
   RA2bBin SB68("Njet9+-Nb3+-MHT1-HT2", Ntries, 9, 99, 3, 9, 200., 500., 800., 1200., true);
   RA2bBin SB69("Njet9+-Nb3+-MHT1-HT3", Ntries, 9, 99, 3, 9, 200., 500., 1200., 9999., true);
   RA2bBin SB70("Njet9+-Nb3+-MHT2-HT1", Ntries, 9, 99, 3, 9, 500., 750., 500., 1200., true);
   RA2bBin SB71("Njet9+-Nb3+-MHT2-HT23", Ntries, 9, 99, 3, 9, 500., 750., 1200., 9999., true);
   RA2bBin SB72("Njet9+-Nb3+-MHT2-HT123", Ntries, 9, 99, 3, 9, 750., 9999., 800., 9999., true);
   
   RA2bBin SB1a("Njet456-MHT1-HT1", Ntries, 4, 6, 0, 99, 200., 500., 500., 800., true);
   RA2bBin SB2a("Njet456-MHT1-HT2", Ntries, 4, 6, 0, 99, 200., 500., 800., 1200., true);
   RA2bBin SB3a("Njet456-MHT1-HT3", Ntries, 4, 6, 0, 99, 200., 500., 1200., 9999., true);
   RA2bBin SB4a("Njet456-MHT2-HT1", Ntries, 4, 6, 0, 99, 500., 750., 500., 1200., true);
   RA2bBin SB5a("Njet456-MHT2-HT23", Ntries, 4, 6, 0, 99, 500., 750., 1200., 9999., true);
   RA2bBin SB6a("Njet456-MHT2-HT123", Ntries, 4, 6, 0, 99, 750., 9999., 800., 9999., true);
   
   RA2bBin SB7a("Njet78-MHT1-HT1", Ntries, 7, 8, 0, 99, 200., 500., 500., 800., true);
   RA2bBin SB8a("Njet78-MHT1-HT2", Ntries, 7, 8, 0, 99, 200., 500., 800., 1200., true);
   RA2bBin SB9a("Njet78-MHT1-HT3", Ntries, 7, 8, 0, 99, 200., 500., 1200., 9999., true);
   RA2bBin SB10a("Njet78-MHT2-HT1", Ntries, 7, 8, 0, 99, 500., 750., 500., 1200., true);
   RA2bBin SB11a("Njet78-MHT2-HT23", Ntries, 7, 8, 0, 99, 500., 750., 1200., 9999., true);
   RA2bBin SB12a("Njet78-MHT2-HT123", Ntries, 7, 8, 0, 99, 750., 9999., 800., 9999., true);
   
   RA2bBin SB13a("Njet9+-MHT1-HT1", Ntries, 9, 99, 0, 99, 200., 500., 500., 800., true);
   RA2bBin SB14a("Njet9+-MHT1-HT2", Ntries, 9, 99, 0, 99, 200., 500., 800., 1200., true);
   RA2bBin SB15a("Njet9+-MHT1-HT3", Ntries, 9, 99, 0, 99, 200., 500., 1200., 9999., true);
   RA2bBin SB16a("Njet9+-MHT2-HT1", Ntries, 9, 99, 0, 99, 500., 750., 500., 1200., true);
   RA2bBin SB17a("Njet9+-MHT2-HT23", Ntries, 9, 99, 0, 99, 500., 750., 1200., 9999., true);
   RA2bBin SB18a("Njet9+-MHT2-HT123", Ntries, 9, 99, 0, 99, 750., 9999., 800., 9999., true);
   
   
   SB.push_back(&SB1);
   SB.push_back(&SB2);
   SB.push_back(&SB3);
   SB.push_back(&SB4);
   SB.push_back(&SB5);
   SB.push_back(&SB6);
   SB.push_back(&SB7);
   SB.push_back(&SB8);
   SB.push_back(&SB9);
   SB.push_back(&SB10);
   SB.push_back(&SB11);
   SB.push_back(&SB12);
   SB.push_back(&SB13);
   SB.push_back(&SB14);
   SB.push_back(&SB15);
   SB.push_back(&SB16);
   SB.push_back(&SB17);
   SB.push_back(&SB18);
   SB.push_back(&SB19);
   SB.push_back(&SB20);
   SB.push_back(&SB21);
   SB.push_back(&SB22);
   SB.push_back(&SB23);
   SB.push_back(&SB24);
   SB.push_back(&SB25);
   SB.push_back(&SB26);
   SB.push_back(&SB27);
   SB.push_back(&SB28);
   SB.push_back(&SB29);
   SB.push_back(&SB30);
   SB.push_back(&SB31);
   SB.push_back(&SB32);
   SB.push_back(&SB33);
   SB.push_back(&SB34);
   SB.push_back(&SB35);
   SB.push_back(&SB36);
   SB.push_back(&SB37);
   SB.push_back(&SB38);
   SB.push_back(&SB39);
   SB.push_back(&SB40);
   SB.push_back(&SB41);
   SB.push_back(&SB42);
   SB.push_back(&SB43);
   SB.push_back(&SB44);
   SB.push_back(&SB45);
   SB.push_back(&SB46);
   SB.push_back(&SB47);
   SB.push_back(&SB48);
   SB.push_back(&SB49);
   SB.push_back(&SB50);
   SB.push_back(&SB51);
   SB.push_back(&SB52);
   SB.push_back(&SB53);
   SB.push_back(&SB54);
   SB.push_back(&SB55);
   SB.push_back(&SB56);
   SB.push_back(&SB57);
   SB.push_back(&SB58);
   SB.push_back(&SB59);
   SB.push_back(&SB60);
   SB.push_back(&SB61);
   SB.push_back(&SB62);
   SB.push_back(&SB63);
   SB.push_back(&SB64);
   SB.push_back(&SB65);
   SB.push_back(&SB66);
   SB.push_back(&SB67);
   SB.push_back(&SB68);
   SB.push_back(&SB69);
   SB.push_back(&SB70);
   SB.push_back(&SB71);
   SB.push_back(&SB72);
   
   //SB.push_back(&SB1a);
   //SB.push_back(&SB2a);
   //SB.push_back(&SB3a);
   //SB.push_back(&SB4a);
   //SB.push_back(&SB5a);
   //SB.push_back(&SB6a);
   //SB.push_back(&SB7a);
   //SB.push_back(&SB8a);
   //SB.push_back(&SB9a);
   //SB.push_back(&SB10a);
   //SB.push_back(&SB11a);
   //SB.push_back(&SB12a);
   //SB.push_back(&SB13a);
   //SB.push_back(&SB14a);
   //SB.push_back(&SB15a);
   //SB.push_back(&SB16a);
   //SB.push_back(&SB17a);
   //SB.push_back(&SB18a);

   
   // ------------------------------------------------------------------------------ //
   
   // get tree with predictions
   cout << "entries prediction tree:" << QCDPrediction.GetEntries() << endl;
   
   // variables from prediction tree
   vtxN = 0;
   NJets = 0;
   BTags = 0;
   NSmear = 0;
   weight = 0;
   HT = 0;
   MHT = 0;
   HT_seed = 0;
   Jet1Pt = 0;
   Jet2Pt = 0;
   Jet3Pt = 0;
   Jet4Pt = 0;
   Jet1Eta = 0;
   Jet2Eta = 0;
   Jet3Eta = 0;
   Jet4Eta = 0;
   DeltaPhi1 = 0;
   DeltaPhi2 = 0;
   DeltaPhi3 = 0;
   minDeltaPhiN = 0;
   
   cout << "Before: SetBranchAddress (prediction)" << endl;
   
   QCDPrediction.SetBranchAddress("NVtx",&vtxN);
   QCDPrediction.SetBranchAddress("NJets",&NJets);
   QCDPrediction.SetBranchAddress("BTags",&BTags);
   QCDPrediction.SetBranchAddress("Ntries",&NSmear);
   QCDPrediction.SetBranchAddress("Weight",&weight);
   QCDPrediction.SetBranchAddress("HT",&HT);
   QCDPrediction.SetBranchAddress("MHT",&MHT);
   QCDPrediction.SetBranchAddress("HT_seed",&HT_seed);
   QCDPrediction.SetBranchAddress("Jet1Pt",&Jet1Pt);
   QCDPrediction.SetBranchAddress("Jet2Pt",&Jet2Pt);
   QCDPrediction.SetBranchAddress("Jet3Pt",&Jet3Pt);
   QCDPrediction.SetBranchAddress("Jet4Pt",&Jet4Pt);
   QCDPrediction.SetBranchAddress("Jet1Eta",&Jet1Eta);
   QCDPrediction.SetBranchAddress("Jet2Eta",&Jet2Eta);
   QCDPrediction.SetBranchAddress("Jet3Eta",&Jet3Eta);
   QCDPrediction.SetBranchAddress("Jet4Eta",&Jet4Eta);
   QCDPrediction.SetBranchAddress("DeltaPhi1",&DeltaPhi1);
   QCDPrediction.SetBranchAddress("DeltaPhi2",&DeltaPhi2);
   QCDPrediction.SetBranchAddress("DeltaPhi3",&DeltaPhi3);
   QCDPrediction.SetBranchAddress("minDeltaPhiN",&minDeltaPhiN);
   
   cout << "After: SetBranchAddress (prediction)" << endl;
   
   int Prediction_entries = 0;
   float smear_rep = 0;
   
   // loop over entries and fill prediction histos
   ULong_t nentries = QCDPrediction.GetEntries();
   
   for ( ULong_t i = 0 ; i<nentries ; i++) {
      
      QCDPrediction.GetEntry(i);
      
      if( i%1000000 == 0 ) std::cout << "event (prediction): " << i << '\n';
      
      bool mdp = DeltaPhiCut_prediction();

      if( HT > 500. && MHT > 200. && NJets >= 4 && vtxN >= 0 && mdp) {
         
         for (int it = 0; it < SB.size(); ++it){
            if (SB.at(it)->FillPred((int)(NSmear-1), (int)NJets, (int)BTags, (double)MHT, (double)HT, mdp, (double)weight)){
               //cout << "HT: " << HT;
               //cout << ", MHT: " << MHT;
               //cout << ", NJets: " << NJets;
               //cout << ", BTags: " << BTags;
               //cout << ", MDP: " << mdp;
               //cout << ", weight: " << weight << endl;
               //break;
            };
         }
         
      }
   }
   cout << "Prediction entries final: " <<  Prediction_entries << endl;
   
   // ------------------------------------------------------------- //
   
   ///////////////////////////////////////////////////////////////////////
   // get event selections
   cout << "entries selection tree:" << RA2PreSelection.GetEntries() << endl;
   
   // variables from selection tree
   vtxN_RA2 = 0;
   NLeptons_RA2 = 0;
   NJets_RA2 = 0;
   BTags_RA2 = 0;
   weight_RA2 = 0;
   HT_RA2 = 0;
   MHT_RA2 = 0;
   Jet1Pt_RA2 = 0;
   Jet2Pt_RA2 = 0;
   Jet3Pt_RA2 = 0;
   Jet4Pt_RA2 = 0;
   Jet1Eta_RA2 = 0;
   Jet2Eta_RA2 = 0;
   Jet3Eta_RA2 = 0;
   Jet4Eta_RA2 = 0;
   DeltaPhi1_RA2 = 0;
   DeltaPhi2_RA2 = 0;
   DeltaPhi3_RA2 = 0;
   minDeltaPhiN_RA2 = 0;
   
   cout << "Before: SetBranchAddress (expectation)" << endl;
   
   RA2PreSelection.SetBranchAddress("NVtx",&vtxN_RA2);
   RA2PreSelection.SetBranchAddress("GoodLeptons",&NLeptons_RA2);
   RA2PreSelection.SetBranchAddress("NJets",&NJets_RA2);
   RA2PreSelection.SetBranchAddress("BTags",&BTags_RA2);
   RA2PreSelection.SetBranchAddress("Weight",&weight_RA2);
   RA2PreSelection.SetBranchAddress("HT",&HT_RA2);
   RA2PreSelection.SetBranchAddress("MHT",&MHT_RA2);
   RA2PreSelection.SetBranchAddress("Jet1Pt",&Jet1Pt_RA2);
   RA2PreSelection.SetBranchAddress("Jet2Pt",&Jet2Pt_RA2);
   RA2PreSelection.SetBranchAddress("Jet3Pt",&Jet3Pt_RA2);
   RA2PreSelection.SetBranchAddress("Jet4Pt",&Jet4Pt_RA2);
   RA2PreSelection.SetBranchAddress("Jet1Eta",&Jet1Eta_RA2);
   RA2PreSelection.SetBranchAddress("Jet2Eta",&Jet2Eta_RA2);
   RA2PreSelection.SetBranchAddress("Jet3Eta",&Jet3Eta_RA2);
   RA2PreSelection.SetBranchAddress("Jet4Eta",&Jet4Eta_RA2);
   RA2PreSelection.SetBranchAddress("DeltaPhi1",&DeltaPhi1_RA2);
   RA2PreSelection.SetBranchAddress("DeltaPhi2",&DeltaPhi2_RA2);
   RA2PreSelection.SetBranchAddress("DeltaPhi3",&DeltaPhi3_RA2);
   RA2PreSelection.SetBranchAddress("minDeltaPhiN",&minDeltaPhiN_RA2);
   
   cout << "After: SetBranchAddress (expectation)" << endl;
   
   // loop over entries and fill selection histos
   Int_t nentries2 = RA2PreSelection.GetEntries();
   
   for ( Int_t i = 0 ; i<nentries2 ; i++) {
      
      RA2PreSelection.GetEntry(i);
      
      if( i%100000 == 0 ) std::cout << "event (selection): " << i << '\n';
      
      bool mdp_RA2 = DeltaPhiCut_selection();
      if( HT_RA2 > 500. && MHT_RA2 > 200. && NJets_RA2 >= 4 && NLeptons_RA2 == 0 && vtxN_RA2 >= 0 && mdp_RA2) {
         
         for (int it = 0; it < SB.size(); ++it){
            if (SB.at(it)->FillExp((int) NJets_RA2, (int)BTags_RA2, (double)MHT_RA2, (double)HT_RA2, mdp_RA2, (double)weight_RA2));
         }
         
      }
   }
   
   //----------------------------------------------------------//
   
   const int NBins = SB.size();
   std::vector<std::string> binname;
   binname.resize(NBins);
   std::vector<double> exp;
   exp.resize(NBins);
   std::vector<double> exp_err;
   exp_err.resize(NBins);
   std::vector<double> pred;
   pred.resize(NBins);
   std::vector<double> pred_err;
   pred_err.resize(NBins);
   
   cout << "Variation: " << Uncertainty << endl;
   cout << "No." << ", " << "BinName" << ", " << "mean_pred" << ", " << "var_pred"  << ", " << "mean_exp" << ", " << "var_exp" << endl;
   for (int it = 0; it < SB.size(); ++it){
      double mean_pred, mean_exp;
      double var_pred, var_exp;
      if ( SB.at(it)->GetPred(mean_pred, var_pred) ){
         cout << it << ", " << SB.at(it)->BinName << ", " << mean_pred << ", " << var_pred;
         pred.at(it) = mean_pred;
         pred_err.at(it) = var_pred;
         binname.at(it) = SB.at(it)->BinName;
      }
      if ( SB.at(it)->GetExp(mean_exp, var_exp) ){
         cout << ", " << mean_exp << ", " << var_exp << endl;
         exp.at(it) = mean_exp;
         exp_err.at(it) = var_exp;
      }
   }
   
   plotClosure(NBins, binname, exp, exp_err, pred, pred_err);
   
}

////////////////////////////////////////////////////////////////////////////////////////
bool BinPrediction::DeltaPhiCut_prediction()
{
   bool deltaPhiCut = true;
   if( NJets == 2 ) {
      if( DeltaPhi1 < 0.5 ||
         DeltaPhi2 < 0.5 ) deltaPhiCut = false;
   }
   if( NJets >= 3 ) {
      if( DeltaPhi1 < 0.5 ||
         DeltaPhi2 < 0.5 ||
         DeltaPhi3 < 0.3 ) deltaPhiCut = false;
   }
   //   if (minDeltaPhiN < 6) deltaPhiCut = false;
   
   return deltaPhiCut;
}
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
bool BinPrediction::DeltaPhiCut_selection()
{
   bool deltaPhiCut = true;
   if( NJets_RA2 == 2 ) {
      if( DeltaPhi1_RA2 < 0.5 ||
         DeltaPhi2_RA2 < 0.5 ) deltaPhiCut = false;
   }
   if( NJets_RA2 >= 3 ) {
      if( DeltaPhi1_RA2 < 0.5 ||
         DeltaPhi2_RA2 < 0.5 ||
         DeltaPhi3_RA2 < 0.3 ) deltaPhiCut = false;
   }
   //   if (minDeltaPhiN_RA2 < 6) deltaPhiCut = false;
   
   return deltaPhiCut;
}
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
void BinPrediction::plotClosure(int NBins, std::vector<std::string> binname, std::vector<double> exp, std::vector<double> exp_err, std::vector<double> pred, std::vector<double> pred_err) {
   
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
   
   enum {cqcd=kRed+3, cltau=kYellow, cllep=kRed+1, czinv=kGreen+1, cdata=kBlack, cpred=kBlue+2};
   
   TH1F *h_pred  = new TH1F("h_pred",  "h_pred",  NBins, 0, NBins);
   TH1F *h_exp = new TH1F("h_exp", "h_exp", NBins, 0, NBins);
   
   h_pred->SetTitle("");
   h_exp->SetTitle("");
   
   h_pred  ->SetLineColor(cqcd );
   h_exp ->SetLineColor(cdata);
   
   h_pred  ->SetFillColor(cqcd );
   h_exp ->SetFillColor(cdata);
   
   h_pred  ->SetMarkerColor(cqcd );
   h_exp ->SetMarkerColor(cdata);
   
   h_exp ->SetMarkerStyle(20);
   h_exp ->SetMarkerSize (1.2);
   h_exp ->SetLineWidth(2);
   
   h_exp->GetYaxis()->SetTitleSize(0.08);
   h_exp->GetYaxis()->SetLabelSize(0.06);
   
   h_pred->GetYaxis()->SetTitleSize(0.08);
   h_pred->GetYaxis()->SetLabelSize(0.06);
   
   //gStyle->SetHatchesSpacing(0.5);
   gStyle->SetHatchesLineWidth(1);
   
   for (int ibin=1; ibin<=h_pred->GetNbinsX(); ibin++){
      h_pred ->SetBinContent(ibin, pred.at(ibin-1));
      h_pred->SetBinError  (ibin, pred_err.at(ibin-1));
      h_pred->GetXaxis()->SetBinLabel(ibin, binname.at(ibin-1).c_str());
      h_exp->SetBinContent(ibin, exp.at(ibin-1));
      h_exp->SetBinError(ibin, exp_err.at(ibin-1));
      h_exp->GetXaxis()->SetBinLabel(ibin, binname.at(ibin-1).c_str());
   }
   
   TCanvas *c = new TCanvas("c", "c", 1200, 800);
   c->cd();
   TPad *pad1 = new TPad("pad1", "pad1", 0, 0.45, 1, 1);
   pad1->SetFillStyle(4000);
   pad1->Draw();
   pad1->SetLogy();
   pad1->SetTopMargin(0.1);
   pad1->SetBottomMargin(0);
   pad1->SetRightMargin(0.05);
   pad1->SetLeftMargin(0.15);
   pad1->cd();
   
   TH1F *vFrame1 = pad1->DrawFrame(0.0, 0.1, (double) NBins, 20000.0);
   vFrame1->GetYaxis()->SetTitle("Events");
   vFrame1->GetYaxis()->SetTitleOffset(-1.3);
   vFrame1->GetYaxis()->SetTitleSize(0.15);
   vFrame1->GetYaxis()->SetTitleFont(42);
   vFrame1->GetYaxis()->SetLabelOffset(-0.04);
   vFrame1->GetYaxis()->SetLabelSize(0.05);
   vFrame1->GetYaxis()->SetLabelFont(42);
   vFrame1->GetXaxis()->SetLabelOffset(1.0);
   vFrame1->GetYaxis()->SetTickLength(0.02);
   vFrame1->GetYaxis()->SetTicks("+");
   vFrame1->SetFillStyle(4000);
   h_pred->Draw("same histe");
   h_exp->Draw("same pe");
   
   TLatex *title = new TLatex(0., 30000, "Simulation, L = 10 fb^{ -1}, #sqrt{s} = 13 TeV");
   title->SetNDC(0);
   title->SetTextFont(42);
   title->SetTextSize(0.06);
   title->Draw("same");
   
   pad1->cd();
   c->cd();
   TLegend *leg1 = new TLegend(0.7,0.85,0.95,0.93,NULL,"NDC");
   leg1->SetLineColor(0);
   leg1->SetLineStyle(1);
   leg1->SetLineWidth(1);
   leg1->SetFillColor(0);
   leg1->SetFillStyle(4000);
   leg1->SetTextSize(0.025);
   leg1->SetTextFont(42);
   leg1->AddEntry(h_exp, "MC Expectation", "lp");
   leg1->AddEntry(h_pred,  "R+S Prediction", "lf");
   leg1->Draw("same");
   
   double x[NBins];
   double y[NBins];
   double ey_UP[NBins];
   double ey_DN[NBins];
   double r[NBins];
   double r_e[NBins];
   double bins_e[NBins];
   for(int ii=0; ii<NBins; ii++){
      bins_e[ii] = 0.5;
      x[ii] = 0.5+ii;
      y[ii] = 0.;
      ey_UP[ii] = pred_err[ii]/pred[ii];
      ey_DN[ii] = pred_err[ii]/pred[ii];
      r[ii] = (pred[ii]-exp[ii])/exp[ii];
      r_e[ii] = sqrt((pred_err[ii]/exp[ii])*(pred_err[ii]/exp[ii])+(exp_err[ii]*pred[ii]/exp[ii]/exp[ii])*(exp_err[ii]*pred[ii]/exp[ii]/exp[ii]));
   }
   
   TGraphAsymmErrors *ratio = new TGraphAsymmErrors(NBins, x, r, bins_e, bins_e, r_e, r_e);
   
   ratio->SetMarkerColor(cdata);
   ratio->SetMarkerStyle(20);
   ratio->SetMarkerSize (1.2);
   ratio->SetLineColor(cdata);
   ratio->SetLineWidth(2);
   
   TH1F *h_ratio = new TH1F("h_ratio_mc", "h_ratio_mc", NBins, 0, NBins);
   
   h_ratio->SetMarkerColor(cdata);
   
   h_ratio->SetMarkerStyle(20);
   h_ratio->SetMarkerSize (1.2);
   h_ratio->SetLineWidth(1);
   h_ratio->SetMinimum(-1.2);
   h_ratio->SetMaximum(2.2);
   
   for (int ibin=1; ibin<=h_exp->GetNbinsX(); ibin++){
      h_ratio->SetBinContent(ibin, 0);
      h_ratio->SetBinError(ibin, 0);
      h_ratio->GetXaxis()->SetBinLabel(ibin, binname.at(ibin-1).c_str());
   }
   
   c->cd();
   TPad *pad2 = new TPad("pad2", "pad2", 0., 0., 1, 0.5);
   pad2->SetTopMargin(0.0);
   pad2->SetRightMargin(0.05);
   pad2->SetLeftMargin(0.15);
   pad2->SetBottomMargin(0.4);
   pad2->Draw("same");
   pad2->cd();
   
   TH1F *vFrame2 = pad2->DrawFrame(0.0, 0.45, (double) NBins, 1.5);
   vFrame2->GetYaxis()->SetTitle("(Pred.-Exp.)/Exp.");
   vFrame2->GetYaxis()->SetTitleOffset(2.6);
   vFrame2->GetYaxis()->SetTitleSize(0.08);
   vFrame2->GetYaxis()->SetTitleFont(42);
   vFrame2->GetYaxis()->SetLabelOffset(0.08);
   vFrame2->GetYaxis()->SetLabelSize(0.06);
   vFrame2->GetYaxis()->SetLabelFont(42);
   vFrame2->GetXaxis()->SetLabelOffset(1.0);
   vFrame2->GetYaxis()->SetTickLength(-0.02);
   vFrame2->GetYaxis()->SetTicks("+");
   
   
   h_ratio->SetTitle("");
   h_ratio->GetYaxis()->SetTitle("(Pred.-Exp.)/Exp.");
   h_ratio->GetYaxis()->SetTitleOffset(0.6);
   h_ratio->GetYaxis()->SetLabelOffset(0.02);
   h_ratio->GetYaxis()->SetTickLength(0.02);
   h_ratio->GetYaxis()->SetTitleSize(0.08);
   h_ratio->GetYaxis()->SetLabelSize(0.06);
   h_ratio->GetXaxis()->SetLabelOffset(0.01);
   h_ratio->Draw("e5");
   ratio->Draw("P0Z");
   
   c->SaveAs("BinByBinClosure_MGMLM_bestMatching_withCleverRBcorr_pt10_angResNew_withNeutrinos.pdf");
   //c->SaveAs("BinByBinClosure_test.pdf");
   
}



