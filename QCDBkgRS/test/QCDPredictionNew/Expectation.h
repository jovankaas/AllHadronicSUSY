//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 21 16:50:14 2015 by ROOT version 6.02/05
// from TTree PreSelection/PreSelection
// found on file: ../QCDSmearingClosure_OnMC.root
//////////////////////////////////////////////////////////

#ifndef Expectation_h
#define Expectation_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1F.h>
#include <TH2F.h>
#include "RA2bBin.h"

// Header file for the classes stored in the TTree if any.

class Expectation : public TSelector {
   public :
   TTree *fChain;   //!pointer to the analyzed TTree or TChain
   
   // Fixed size dimensions of array or collections stored in the TTree if any.
   
   // Declaration of leaf types
   UInt_t RunNum;
   UInt_t LumiBlockNum;
   UInt_t EvtNum;
   Int_t NVtx;
   Int_t GoodLeptons;
   Int_t NJets;
   Int_t BTags;
   Float_t HT;
   Float_t minDeltaPhiN;
   Float_t MHT;
   Float_t Jet1Pt;
   Float_t Jet2Pt;
   Float_t Jet3Pt;
   Float_t Jet4Pt;
   Float_t Jet1Eta;
   Float_t Jet2Eta;
   Float_t Jet3Eta;
   Float_t Jet4Eta;
   Float_t DeltaPhi1;
   Float_t DeltaPhi2;
   Float_t DeltaPhi3;
   Float_t Weight;
   
   // List of branches
   TBranch *b_RunNum;   //!
   TBranch *b_LumiBlockNum;   //!
   TBranch *b_EvtNum;   //!
   TBranch *b_NVtx;   //!
   TBranch *b_GoodLeptons;   //!
   TBranch *b_NJets;   //!
   TBranch *b_BTags;   //!
   TBranch *b_HT;   //!
   TBranch *b_minDeltaPhiN;   //!
   TBranch *b_MHT;   //!
   TBranch *b_Jet1Pt;   //!
   TBranch *b_Jet2Pt;   //!
   TBranch *b_Jet3Pt;   //!
   TBranch *b_Jet4Pt;   //!
   TBranch *b_Jet1Eta;   //!
   TBranch *b_Jet2Eta;   //!
   TBranch *b_Jet3Eta;   //!
   TBranch *b_Jet4Eta;   //!
   TBranch *b_DeltaPhi1;   //!
   TBranch *b_DeltaPhi2;   //!
   TBranch *b_DeltaPhi3;   //!
   TBranch *b_Weight;   //!
   
   Expectation(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~Expectation() { }
   virtual Int_t Version() const { return 2; }
   virtual void Begin(TTree *tree);
   virtual void SlaveBegin(TTree *tree);
   virtual void Init(TTree *tree);
   virtual Bool_t Notify();
   virtual Bool_t Process(Long64_t entry);
   virtual Int_t GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void SetOption(const char *option) { fOption = option; }
   virtual void SetObject(TObject *obj) { fObject = obj; }
   virtual void SetInputList(TList *input) { fInput = input; }
   virtual TList *GetOutputList() const { return fOutput; }
   virtual void SlaveTerminate();
   virtual void Terminate();
   virtual bool DeltaPhiCut();
   
   std::vector<RA2bBin*> SB;
   TH1F* yields;
   int TotEvents;

   ClassDef(Expectation,0);
};

#endif

#ifdef Expectation_cxx
void Expectation::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
   
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);
   
   fChain->SetBranchAddress("RunNum", &RunNum, &b_RunNum);
   fChain->SetBranchAddress("LumiBlockNum", &LumiBlockNum, &b_LumiBlockNum);
   fChain->SetBranchAddress("EvtNum", &EvtNum, &b_EvtNum);
   fChain->SetBranchAddress("NVtx", &NVtx, &b_NVtx);
   fChain->SetBranchAddress("GoodLeptons", &GoodLeptons, &b_GoodLeptons);
   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("BTags", &BTags, &b_BTags);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("minDeltaPhiN", &minDeltaPhiN, &b_minDeltaPhiN);
   fChain->SetBranchAddress("MHT", &MHT, &b_MHT);
   fChain->SetBranchAddress("Jet1Pt", &Jet1Pt, &b_Jet1Pt);
   fChain->SetBranchAddress("Jet2Pt", &Jet2Pt, &b_Jet2Pt);
   fChain->SetBranchAddress("Jet3Pt", &Jet3Pt, &b_Jet3Pt);
   fChain->SetBranchAddress("Jet4Pt", &Jet4Pt, &b_Jet4Pt);
   fChain->SetBranchAddress("Jet1Eta", &Jet1Eta, &b_Jet1Eta);
   fChain->SetBranchAddress("Jet2Eta", &Jet2Eta, &b_Jet2Eta);
   fChain->SetBranchAddress("Jet3Eta", &Jet3Eta, &b_Jet3Eta);
   fChain->SetBranchAddress("Jet4Eta", &Jet4Eta, &b_Jet4Eta);
   fChain->SetBranchAddress("DeltaPhi1", &DeltaPhi1, &b_DeltaPhi1);
   fChain->SetBranchAddress("DeltaPhi2", &DeltaPhi2, &b_DeltaPhi2);
   fChain->SetBranchAddress("DeltaPhi3", &DeltaPhi3, &b_DeltaPhi3);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
}

Bool_t Expectation::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.
   
   return kTRUE;
}

#endif // #ifdef Expectation_cxx
