#define Expectation_cxx
// The class definition in Expectation.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("Expectation.C")
// root> T->Process("Expectation.C","some options")
// root> T->Process("Expectation.C+")
//

#include "Expectation.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>


void Expectation::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).
   HT_ = NULL;
	 MHT_ = NULL;
	 NJets_= NULL;
	 BTags_ = NULL;
	 TH1Prediction_=NULL;

   TString option = GetOption();

}

void Expectation::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
   SearchBin_= new SearchBins();
	 std::cout<<"Searchregions: "<<SearchBin_->GetAmountOfBins()<<"endl";
	 HT_ = new TH1F("HT","HT",500,0,5000);
	 GetOutputList()->Add(HT_);
	 MHT_ = new TH1F("MHT","MHT",300,0,3000);
	 GetOutputList()->Add(MHT_);
	 NJets_ = new TH1F("NJets","NJets",15, 0, 15);
	 GetOutputList()->Add(NJets_);
	 BTags_ = new TH1F("BTags","BTags",15, 0, 15);
	 GetOutputList()->Add(BTags_);
	 TH1Prediction_ = new TH1F("Bin","Bin",SearchBin_->GetAmountOfBins(),0.5,SearchBin_->GetAmountOfBins()+0.5);
	 GetOutputList()->Add(TH1Prediction_);

   TString option = GetOption();

}

Bool_t Expectation::Process(Long64_t entry)
{
	fChain->GetTree()->GetEntry(entry);
	 	//if(DeltaPhi1 < deltaPhi1_ || DeltaPhi2 < deltaPhi2_ || DeltaPhi3 < deltaPhi3_ )return kTRUE;
	//if( HT < HTCut_ || MHT < MHTCut_ || NJets < NJetsCut_) return kTRUE;
	if( Weight > 30000. ) return kTRUE;
	
	HT_->Fill(HT,Weight);
	MHT_->Fill(MHT,Weight);
	NJets_->Fill(NJets,Weight);
	BTags_->Fill(BTags,Weight);
	unsigned int bin = SearchBin_->GetBinNumber(HT,MHT,NJets,BTags);
	TH1Prediction_->Fill(bin,Weight);
	
	
	
	return kTRUE;
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either Expectation::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
}

void Expectation::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void Expectation::Terminate()
{
	GetOutputList()->Print();
	HT_ = (TH1F*)GetOutputList()->FindObject("HT");
	MHT_ = (TH1F*)GetOutputList()->FindObject("MHT");
	NJets_ = (TH1F*)GetOutputList()->FindObject("NJets");
	BTags_ = (TH1F*)GetOutputList()->FindObject("BTags");
	TH1Prediction_ = (TH1F*)GetOutputList()->FindObject("Bin");
	TFile *outPutFile = new TFile("Expectation.root","RECREATE"); ;
	outPutFile->cd();
	HT_->Write();
	DoRebinning(HT_,0);
	MHT_->Write();
	DoRebinning(MHT_,0);
	NJets_->Write();
	DoRebinning(NJets_,0);
	BTags_->Write();
	DoRebinning(BTags_,0);
	TH1Prediction_->Write();
	DoRebinning(TH1Prediction_,0);
	TTree *result  =  new TTree("QCDExpectation","Expectation");
	Float_t prediction;
	Float_t uncertainty;
	UShort_t bin;
	result->Branch("Bin",&bin,"Bin/s");
	result->Branch("Uncertainty",&uncertainty,"Uncertainty/F");
	result->Branch("Integral",&prediction,"Integral/F");
	for(int i=0; i< ((TH1F*)GetOutputList()->FindObject("Bin"))->GetXaxis()->GetNbins();i++)
	{
		bin=i+1;
		prediction= ((TH1F*)GetOutputList()->FindObject("Bin"))->GetBinContent(i+1);
		uncertainty = ((TH1F*)GetOutputList()->FindObject("Bin"))->GetBinError(i+1);
		result->Fill();
	}
	result->Write();
	outPutFile->Close();
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}



SearchBins::SearchBins()
{
	
	// HTmin,HTmax,MHTmin,MHTmax,NJetsmin,NJetsmax,BTagsmin,BTagsmax
	// NJets 4,6 BTags=0
	// fixed ht Njets btags all MHT bins
	
	bins_.push_back( Bin(500,800,200,500,4,6,-1,0) );
	bins_.push_back( Bin(800,1200,200,500,4,6,-1,0) );
	bins_.push_back( Bin(1200,99999,200,500,4,6,-1,0) );
	
	bins_.push_back( Bin(500,1200,500,750,4,6,-1,0) );
	bins_.push_back( Bin(1200,99999,500,750,4,6,-1,0) );
	
	bins_.push_back( Bin(800,99999,750,9999,4,6,-1,0) );
	
	// NJets 4,6 BTags=1
	// fixed ht Njets btags all MHT bins
	bins_.push_back( Bin(500,800,200,500,4,6,1,1) );
	bins_.push_back( Bin(800,1200,200,500,4,6,1,1) );
	bins_.push_back( Bin(1200,99999,200,500,4,6,1,1) );
	
	bins_.push_back( Bin(500,1200,500,750,4,6,1,1) );
	bins_.push_back( Bin(1200,99999,500,750,4,6,1,1) );
	
	bins_.push_back( Bin(800,99999,750,9999,4,6,1,1) );
	
	// NJets 4,6 BTags=2
	// fixed ht Njets btags all MHT bins
	bins_.push_back( Bin(500,800,200,500,4,6,2,2) );
	bins_.push_back( Bin(800,1200,200,500,4,6,2,2) );
	bins_.push_back( Bin(1200,99999,200,500,4,6,2,2) );
	
	bins_.push_back( Bin(500,1200,500,750,4,6,2,2) );
	bins_.push_back( Bin(1200,99999,500,750,4,6,2,2) );
	
	bins_.push_back( Bin(800,99999,750,9999,4,6,2,2) );
	
	// NJets 4,6 BTags=>3
	// fixed ht Njets btags all MHT bins
	bins_.push_back( Bin(500,800,200,500,4,6,3,9999) );
	bins_.push_back( Bin(800,1200,200,500,4,6,3,9999) );
	bins_.push_back( Bin(1200,99999,200,500,4,6,3,9999) );
	
	bins_.push_back( Bin(500,1200,500,750,4,6,3,9999) );
	bins_.push_back( Bin(1200,99999,500,750,4,6,3,9999) );
	
	bins_.push_back( Bin(800,99999,750,9999,4,6,3,9999) );
	
	// NJewts 7,8 BTags=0
	bins_.push_back( Bin(500,800,200,500,7,8,-1,0) );
	bins_.push_back( Bin(800,1200,200,500,7,8,-1,0) );
	bins_.push_back( Bin(1200,99999,200,500,7,8,-1,0) );
	
	bins_.push_back( Bin(500,1200,500,750,7,8,-1,0) );
	bins_.push_back( Bin(1200,99999,500,750,7,8,-1,0) );
	
	bins_.push_back( Bin(800,99999,750,9999,7,8,-1,0) );
	
	// NJewts 7,8 BTags=1
	bins_.push_back( Bin(500,800,200,500,7,8,1,1) );
	bins_.push_back( Bin(800,1200,200,500,7,8,1,1) );
	bins_.push_back( Bin(1200,99999,200,500,7,8,1,1) );
	
	bins_.push_back( Bin(500,1200,500,750,7,8,1,1) );
	bins_.push_back( Bin(1200,99999,500,750,7,8,1,1) );
	
	bins_.push_back( Bin(800,99999,750,9999,7,8,1,1) );
	
	// NJewts 7,8 BTags=2
	bins_.push_back( Bin(500,800,200,500,7,8,2,2) );
	bins_.push_back( Bin(800,1200,200,500,7,8,2,2) );
	bins_.push_back( Bin(1200,99999,200,500,7,8,2,2) );
	
	bins_.push_back( Bin(500,1200,500,750,7,8,2,2) );
	bins_.push_back( Bin(1200,99999,500,750,7,8,2,2) );
	
	bins_.push_back( Bin(800,99999,750,9999,7,8,2,2) );
	
	// NJewts 7,8 BTags=>3
	bins_.push_back( Bin(500,800,200,500,7,8,3,9999) );
	bins_.push_back( Bin(800,1200,200,500,7,8,3,9999) );
	bins_.push_back( Bin(1200,99999,200,500,7,8,3,9999) );
	
	bins_.push_back( Bin(500,1200,500,750,7,8,3,9999) );
	bins_.push_back( Bin(1200,99999,500,750,7,8,3,9999) );
	
	bins_.push_back( Bin(800,99999,750,9999,7,8,3,9999) );
	
	
	// NJewts 9,9999 BTags=0
	bins_.push_back( Bin(500,800,200,500,9,9999,-1,0) );
	bins_.push_back( Bin(800,1200,200,500,9,9999,-1,0) );
	bins_.push_back( Bin(1200,99999,200,500,9,9999,-1,0) );
	
	bins_.push_back( Bin(500,1200,500,750,9,9999,-1,0) );
	bins_.push_back( Bin(1200,99999,500,750,9,9999,-1,0) );
	
	bins_.push_back( Bin(800,99999,750,9999,9,9999,-1,0) );
	
	
	// NJewts 9,9999 BTags=1
	bins_.push_back( Bin(500,800,200,500,9,9999,1,1) );
	bins_.push_back( Bin(800,1200,200,500,9,9999,1,1) );
	bins_.push_back( Bin(1200,99999,200,500,9,9999,1,1) );
	
	bins_.push_back( Bin(500,1200,500,750,9,9999,1,1) );
	bins_.push_back( Bin(1200,99999,500,750,9,9999,1,1) );
	
	bins_.push_back( Bin(800,99999,750,9999,9,9999,1,1) );
	
	
	// NJewts 9,9999 BTags=2
	bins_.push_back( Bin(500,800,200,500,9,9999,2,2) );
	bins_.push_back( Bin(800,1200,200,500,9,9999,2,2) );
	bins_.push_back( Bin(1200,99999,200,500,9,9999,2,2) );
	
	bins_.push_back( Bin(500,1200,500,750,9,9999,2,2) );
	bins_.push_back( Bin(1200,99999,500,750,9,9999,2,2) );
	
	bins_.push_back( Bin(800,99999,750,9999,9,9999,2,2) );
	
	
	// NJewts 9,9999 BTags=>3
	bins_.push_back( Bin(500,800,200,500,9,9999,3,9999) );
	bins_.push_back( Bin(800,1200,200,500,9,9999,3,9999) );
	bins_.push_back( Bin(1200,99999,200,500,9,9999,3,9999) );
	
	bins_.push_back( Bin(500,1200,500,750,9,9999,3,9999) );
	bins_.push_back( Bin(1200,99999,500,750,9,9999,3,9999) );
	
	bins_.push_back( Bin(800,99999,750,9999,9,9999,3,9999) );
	
	for(unsigned int i=0; i<bins_.size();i++)
	{
		usedBin_.push_back(0); 
	}
}

unsigned int SearchBins::GetBinNumber(double HT, double MHT, int NJets, int BTags)
{
	unsigned int result =999;
	int match =-1;
	for(unsigned int i=0; i<bins_.size();i++)
	{
		//              std::cout<<"Bin["<<i<<"]: HT["<<bins_[i].HTmin_<<","<<bins_[i].HTmax_<<"] MHT["<<bins_[i].MHTmin_<<","<<bins_[i].MHTmax_<<"] NJets["<<bins_[i].NJetsmin_<<","<<bins_[i].NJetsmax_<<"] BTags["<<bins_[i].BTagsmin_<<","<<bins_[i].BTagsmax_<<"]\n";
		if(HT>bins_[i].HTmin_ && 
			HT<bins_[i].HTmax_ &&
			MHT>bins_[i].MHTmin_ && 
			MHT<bins_[i].MHTmax_ &&
			NJets+0.1>bins_[i].NJetsmin_ && 
			NJets-0.1<bins_[i].NJetsmax_ &&
			BTags+0.1>bins_[i].BTagsmin_ && 
			BTags-0.1<bins_[i].BTagsmax_
		)
		{
			result=i;
			match++;
			usedBin_[i]=usedBin_[i]+1;
		}
	}
	if(match==-1)
	{
		//     std::cout<<"Error event fits in no bin!!! HT: "<<HT<<", MHT: "<<MHT<<", NJets: "<<NJets<<", BTags: "<<BTags<<std::endl;
		result=999;
	}
	if(match>0)
	{
		std::cout<<"Error event fits in more than one bin!!!! HT: "<<HT<<", MHT: "<<MHT<<", NJets: "<<NJets<<", BTags: "<<BTags<<std::endl;
	}
	return result+1; // to not use the 0 bin but start counting at 1
}
void SearchBins::PrintUsed()
{
	for(unsigned int i=0; i< usedBin_.size();i++) std::cout<<"Bin["<<i<<"] has been used: "<<usedBin_[i]<<std::endl; 
}

void Expectation::DoRebinning(TH1F* selection_raw, int Nbins)
{
	//do some non-equidistant binning
	if (Nbins < 0) {
		int nbins = 19;
		// HT binning
		double vbins[20] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1300, 1500, 1700, 2000, 2500, 3000, 4000, 5000};
		
		// MHT binning
		if (Nbins == -2) {
			nbins = 17;
			vbins[0]  = 0.;
			vbins[1]  = 10.;
			vbins[2]  = 20.;
			vbins[3]  = 30.;
			vbins[4]  = 40.;
			vbins[5]  = 50.;
			vbins[6]  = 60.;
			vbins[7]  = 70.;
			vbins[8]  = 80.;
			vbins[9] = 90.;
			vbins[10] = 100.;
			vbins[11] = 120.;
			vbins[12] = 150.;
			vbins[13] = 200.;
			vbins[14] = 300.;
			vbins[15] = 500.;
			vbins[16] = 800.;
			vbins[17] = 1500.;
		}
		
		TH1F* temp2 = (TH1F*) selection_raw->Clone();
		selection_raw->GetXaxis()->Set(nbins, &vbins[0]);
		for (int i = 0; i <= selection_raw->GetXaxis()->GetNbins() + 1; ++i) {
			selection_raw->SetBinContent(i, 0);
			selection_raw->SetBinError(i, 0);
		}
		
		int bin = 0;
		double sum2 = 0., content = 0.;
		for (int i = 0; i <= temp2->GetXaxis()->GetNbins() + 1; ++i) {
			int this_bin = selection_raw->GetXaxis()->FindBin(temp2->GetXaxis()->GetBinCenter(i));
			if (this_bin > bin) {
				double binWidth = selection_raw->GetXaxis()->GetBinWidth(bin);
				selection_raw->SetBinContent(bin, content/binWidth);
				selection_raw->SetBinError(bin, sqrt(sum2)/binWidth);
				bin = this_bin;
				sum2 = content = 0.;
			}
			sum2 += pow(temp2->GetBinError(i), 2);
			content += temp2->GetBinContent(i);
		}
		
	} else if(Nbins>0) {
		selection_raw->Rebin(Nbins);
	}
}

