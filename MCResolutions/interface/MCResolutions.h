#ifndef MCResolutions_H
#define MCResolutions_H

// system include files
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TFile.h"

// b-Tagging
#include "DataFormats/BTauReco/interface/JetTag.h"

//
// class declaration
//

class MCResolutions: public edm::EDAnalyzer {
   public:
   explicit MCResolutions(const edm::ParameterSet&);
   ~MCResolutions();

   private:
   virtual void beginJob();
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   virtual void endJob();

   // ----------member data ---------------------------
   edm::InputTag _jetTag;
   edm::InputTag _leptonTag;
   std::string _btagTag;
   edm::InputTag _genJetTag;
   edm::InputTag _weightName;

   double _btagCut;
   double _deltaRMatch;
   double _deltaRMatchVeto;
   double _absPtVeto;
   double _relPtVeto;
   double _GenJetPtCut;
   std::string _fileName;

   double weight;

   // JetResponse in Pt and eta bins
   int PtBin(const double& pt);
   int EtaBin(const double& eta);

   // Resize histo vectors
   void ResizeHistoVector(std::vector<std::vector<TH1F*> > &histoVector);

   // total
   std::vector<std::vector<TH1F*> > h_tot_JetAll_JetResPt_Pt;
   std::vector<std::vector<TH1F*> > h_tot_JetAll_JetResPhi_Pt;
   std::vector<std::vector<TH1F*> > h_tot_JetAll_JetResEta_Pt;
   // with btag
   std::vector<std::vector<TH1F*> > h_b_JetAll_JetResPt_Pt;
   std::vector<std::vector<TH1F*> > h_b_JetAll_JetResPhi_Pt;
   std::vector<std::vector<TH1F*> > h_b_JetAll_JetResEta_Pt;
   // without btag
   std::vector<std::vector<TH1F*> > h_nob_JetAll_JetResPt_Pt;
   std::vector<std::vector<TH1F*> > h_nob_JetAll_JetResPhi_Pt;
   std::vector<std::vector<TH1F*> > h_nob_JetAll_JetResEta_Pt;
   //
   // Btag efficiencies
   std::vector<TH1F*> h_trueb_RecoPt;
   std::vector<TH1F*> h_no_trueb_RecoPt;
   std::vector<TH1F*> h_trueb_btag_RecoPt;
   std::vector<TH1F*> h_no_trueb_btag_RecoPt;
   std::vector<TH1F*> h_no_trueb_no_btag_RecoPt;
   std::vector<TH1F*> h_trueb_no_btag_RecoPt;
   std::vector<TH1F*> h_btag_RecoPt;
   std::vector<TH1F*> h_no_btag_RecoPt;

   std::vector<double> PtBinEdges;
   std::vector<double> EtaBinEdges;

   TFile* hfile;

};

#endif
