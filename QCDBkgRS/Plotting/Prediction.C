#define Prediction_cxx
// The class definition in Prediction.h has been generated automatically
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
// root> T->Process("Prediction.C")
// root> T->Process("Prediction.C","some options")
// root> T->Process("Prediction.C+")
//

#include "Prediction.h"
#include <TH2.h>
#include <TStyle.h>
#include <TH2.h>
#include <TStyle.h>
#include <iostream>

void Prediction::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).
   HT_ = NULL;
	 MHT_ = NULL;
	 NJets_= NULL;
	 BTags_ = NULL;
	 TH2Prediction_=NULL;

   TString option = GetOption();

}

void Prediction::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
   SearchBin_= new SearchBins();
	 std::cout<<"Searchregions: "<<SearchBin_->GetAmountOfBins()<<" Ntries: "<< Ntries_<<"endl";
	 HT_ = new TH2F("HT","HT",500,0,5000,Ntries_, 0.5, Ntries_ + 0.5);
	 GetOutputList()->Add(HT_);
	 MHT_ = new TH2F("MHT","MHT",300,0,3000,Ntries_, 0.5, Ntries_ + 0.5);
	 GetOutputList()->Add(MHT_);
	 NJets_ = new TH2F("NJets","NJets",15, 0, 15,Ntries_, 0.5, Ntries_ + 0.5);
	 GetOutputList()->Add(NJets_);
	 BTags_ = new TH2F("BTags","BTags",15, 0, 15,Ntries_, 0.5, Ntries_ + 0.5);
	 GetOutputList()->Add(BTags_);
	 TH2Prediction_ = new TH2F("Bin","Bin",SearchBin_->GetAmountOfBins(),0.5,SearchBin_->GetAmountOfBins()+0.5,Ntries_, 0.5, Ntries_ + 0.5);
	 GetOutputList()->Add(TH2Prediction_);
// 	 GetOutputList()->Add(prediction_);
   TString option = GetOption();

}

Bool_t Prediction::Process(Long64_t entry)
{
	fChain->GetTree()->GetEntry(entry);
 	//if(DeltaPhi1 < deltaPhi1_ || DeltaPhi2 < deltaPhi2_ || DeltaPhi3 < deltaPhi3_ )return kTRUE;
	//if( HT < HTCut_ || MHT < MHTCut_ || NJets < NJetsCut_) return kTRUE;
	if( Weight > 30000. ) return kTRUE;
		
	HT_->Fill(HT,Ntries,Weight);
	MHT_->Fill(MHT,Ntries,Weight);
	NJets_->Fill(NJets,Ntries,Weight);
	BTags_->Fill(BTags,Ntries,Weight);
	unsigned int bin = SearchBin_->GetBinNumber(HT,MHT,NJets,BTags);
	TH2Prediction_->Fill(bin,Ntries,Weight);
	


   return kTRUE;
}

void Prediction::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void Prediction::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   GetOutputList()->Print();
   HT_ = (TH2F*)GetOutputList()->FindObject("HT");
	 MHT_ = (TH2F*)GetOutputList()->FindObject("MHT");
	 NJets_ = (TH2F*)GetOutputList()->FindObject("NJets");
	 BTags_ = (TH2F*)GetOutputList()->FindObject("BTags");
	 TH2Prediction_ = (TH2F*)GetOutputList()->FindObject("Bin");
   TFile *outPutFile = new TFile("Prediction.root","RECREATE");
	 outPutFile->cd();
	 DoRebinning(HT_,0);
	 CalcPrediction(HT_)->Write();
	 DoRebinning(MHT_,0);
	 CalcPrediction(MHT_)->Write();
	 
	 DoRebinning(NJets_,0);
	 CalcPrediction(NJets_)->Write();
	 DoRebinning(BTags_,0);
	 CalcPrediction(BTags_)->Write();
	 
	 TH2Prediction_->Write();
	 DoRebinning(TH2Prediction_,0);
	 TH1F* TH1Prediction = CalcPrediction(TH2Prediction_);
	 TH1Prediction->Write();
	 TTree *result  =  new TTree("QCDPrediction","Prediction");
	 Float_t prediction;
	 Float_t uncertainty;
	 UShort_t bin;
	 result->Branch("Bin",&bin,"Bin/s");
	 result->Branch("Uncertainty",&uncertainty,"Uncertainty/F");
	 result->Branch("Prediction",&prediction,"Prediction/F");
	 for(int i=0; i< TH1Prediction->GetXaxis()->GetNbins();i++)
	 {
		 bin=i+1;
		 prediction= TH1Prediction->GetBinContent(i+1);
		 uncertainty = TH1Prediction->GetBinError(i+1);
		 result->Fill();
	 }
	 result->Write();
	 outPutFile->Close();

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
void Prediction::DoRebinning(TH2F* prediction_raw,  int Nbins)
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
		
		TH2F* temp = (TH2F*) prediction_raw->Clone();
		prediction_raw->GetXaxis()->Set(nbins, &vbins[0]);
		for (int j = 0; j <= prediction_raw->GetYaxis()->GetNbins() + 1; ++j) {
			for (int i = 0; i <= prediction_raw->GetXaxis()->GetNbins() + 1; ++i) {
				prediction_raw->SetBinContent(i, j, 0);
				prediction_raw->SetBinError(i, j, 0);
			}
		}
		
		//loop over y-axis
		for (int j = 1; j < temp->GetYaxis()->GetNbins() + 1; ++j) {
			int bin = 0;
			double sum2 = 0., content = 0.;
			for (int i = 0; i <= temp->GetXaxis()->GetNbins() + 1; ++i) {
				int this_bin = prediction_raw->GetXaxis()->FindBin(temp->GetXaxis()->GetBinCenter(i));
				if (this_bin > bin) {
					double binWidth = prediction_raw->GetXaxis()->GetBinWidth(bin);
					prediction_raw->SetBinContent(bin, j, content/binWidth);
					prediction_raw->SetBinError(bin, j, sqrt(sum2)/binWidth);
					bin = this_bin;
					sum2 = content = 0.;
				}
				sum2 += pow(temp->GetBinError(i, j), 2);
				content += temp->GetBinContent(i, j);
			}
			
		}
		
		
		
	} else if(Nbins>0) {
		prediction_raw->Rebin2D(Nbins, 1);
	}
}

TH1F* Prediction::CalcPrediction(TH2F* prediction_raw) {
	TH1F* prediction = new TH1F();
	TString Name = prediction_raw->GetName();
	prediction = (TH1F*) prediction_raw->ProjectionX();
	prediction->Reset();
	for (int i = 0; i <= prediction_raw->GetXaxis()->GetNbins() + 1; ++i) {
		TH1F h = *((TH1F*) prediction_raw->ProjectionY("py", i, i));
		
		double summ = 0;
		double sumv = 0;
		int N = 0;
		
		//// Calculate mean
		for (int j = 1; j <= h.GetNbinsX(); ++j) {
			//for (int j = 501; j <= 1000; ++j) {
				summ += h.GetBinContent(j);
				++N;
		}
		double mean = summ / N;
		
		//// Calculated variance
		for (int j = 1; j <= h.GetNbinsX(); ++j) {
			//for (int j = 501; j <= 1000; ++j) {
				sumv += pow(mean - h.GetBinContent(j), 2);
		}
		double variance = sqrt(sumv / N);
		//cout << "i, mean, variance: " << i << " " << mean << " " << variance << endl;
		
		prediction->SetBinContent(i, mean);
		prediction->SetBinError(i, variance);
	}
	prediction->SetName(Name);
	return prediction;
}
