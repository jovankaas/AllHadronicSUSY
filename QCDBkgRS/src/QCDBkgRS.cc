// -*- C++ -*-
//
// Package:    QCDBkgRS
// Class:      QCDBkgRS
//
/**\class QCDBkgRS QCDBkgRS.cc RA2Classic/QCDBkgRS/src/QCDBkgRS.cc

 Description: [one line class summary]

 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  Kristin Heine,,,DESY
//         Created:  Tue Aug  7 15:55:02 CEST 2012
// $Id: QCDBkgRS.cc,v 1.6 2013/01/08 11:14:22 kheine Exp $
//
//

#include "AllHadronicSUSY/QCDBkgRS/interface/QCDBkgRS.h"


#include <TROOT.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPostScript.h>


//--------------------------------------------------------------------------
QCDBkgRS::QCDBkgRS(const edm::ParameterSet& iConfig)
{
   rebalancedJetPt_ = iConfig.getParameter<double> ("RebalanceJetPt");
   rebalanceMode_ = iConfig.getParameter<std::string> ("RebalanceMode");
   nSmearedJets_ = iConfig.getParameter<int> ("NSmearedJets");
   smearedJetPt_ = iConfig.getParameter<double> ("SmearedJetPt");
   PtBinEdges_scaling_ = iConfig.getParameter<std::vector<double> > ("PtBinEdges_scaling");
   EtaBinEdges_scaling_ = iConfig.getParameter<std::vector<double> > ("EtaBinEdges_scaling");
   AdditionalSmearing_ = iConfig.getParameter<std::vector<double> > ("AdditionalSmearing");
   LowerTailScaling_ = iConfig.getParameter<std::vector<double> > ("LowerTailScaling");
   UpperTailScaling_ = iConfig.getParameter<std::vector<double> > ("UpperTailScaling");
   AdditionalSmearing_variation_ = iConfig.getParameter<double> ("AdditionalSmearing_variation");
   LowerTailScaling_variation_ = iConfig.getParameter<double> ("LowerTailScaling_variation");
   UpperTailScaling_variation_ = iConfig.getParameter<double> ("UpperTailScaling_variation");
   smearCollection_ = iConfig.getParameter<std::string> ("SmearCollection");
   vertices_ = iConfig.getParameter<edm::InputTag>("VertexCollection");
   genjets_ = iConfig.getParameter<edm::InputTag> ("genjetCollection");
   jets_ = iConfig.getParameter<edm::InputTag> ("jetCollection");
   jets_reb_ = iConfig.getParameter<std::string> ("jetCollection_reb");
   jets_smeared_ = iConfig.getParameter<std::string> ("jetCollection_smeared");
   genjets_smeared_ = iConfig.getParameter<std::string> ("genjetCollection_smeared");
   btagTag_ = iConfig.getParameter<std::string> ("btagTag");
   btagCut_ = iConfig.getParameter<double> ("btagCut");
   leptonTag_ = iConfig.getParameter<edm::InputTag> ("leptonTag");
   uncertaintyName_ = iConfig.getParameter<std::string> ("uncertaintyName");
   inputhistPtHF_ = iConfig.getParameter<std::string> ("InputHistoPt_HF");
   inputhistEtaHF_ = iConfig.getParameter<std::string> ("InputHistoEta_HF");
   inputhistPhiHF_ = iConfig.getParameter<std::string> ("InputHistoPhi_HF");
   inputhistPtNoHF_ = iConfig.getParameter<std::string> ("InputHistoPt_NoHF");
   inputhistEtaNoHF_ = iConfig.getParameter<std::string> ("InputHistoEta_NoHF");
   inputhistPhiNoHF_ = iConfig.getParameter<std::string> ("InputHistoPhi_NoHF");
   RebalanceCorrectionFile_ = iConfig.getParameter<std::string> ("RebalanceCorrectionFile");
   BTagEfficiencyFile_ = iConfig.getParameter<std::string> ("BTagEfficiencyFile");
   controlPlots_ = iConfig.getParameter<bool> ("ControlPlots");
   isData_ = iConfig.getParameter<bool> ("IsData");
   isMadgraph_ = iConfig.getParameter<bool> ("IsMadgraph");
   smearingfile_ = iConfig.getParameter<std::string> ("SmearingFile");
   outputfile_ = iConfig.getParameter<std::string> ("OutputFile");
   NRebin_ = iConfig.getParameter<int> ("NRebin");
   weightName_ = iConfig.getParameter<edm::InputTag> ("weightName");
   absoluteTailScaling_ = iConfig.getParameter<bool> ("absoluteTailScaling");
   cleverPrescaleTreating_ = iConfig.getParameter<bool> ("cleverPrescaleTreating");
   useRebalanceCorrectionFactors_ = iConfig.getParameter<bool> ("useRebalanceCorrectionFactors");
   useBTagEfficiencyFactors_ = iConfig.getParameter<bool> ("useBTagEfficiencyFactors");
   useCleverRebalanceCorrectionFactors_ = iConfig.getParameter<bool> ("useCleverRebalanceCorrectionFactors");
   A0RMS_ = iConfig.getParameter<double> ("A0RMS");
   A1RMS_ = iConfig.getParameter<double> ("A1RMS");
   probExtreme_ = iConfig.getParameter<double> ("probExtreme");
   HTSeedTag_ = iConfig.getParameter<edm::InputTag> ("HTSeedTag");
   HTSeedMin_ = iConfig.getParameter<double> ("HTSeedMin");
   NJetsSeedTag_ = iConfig.getParameter<edm::InputTag> ("NJetsSeedTag");
   NJetsSeedMin_ = iConfig.getParameter<int> ("NJetsSeedMin");
   MHTmin_ = iConfig.getParameter<double> ("MHTmin");
   MHTmax_ = iConfig.getParameter<double> ("MHTmax");
   HTmin_ = iConfig.getParameter<double> ("HTmin");
   HTmax_ = iConfig.getParameter<double> ("HTmax");
   NbinsMHT_ = iConfig.getParameter<int> ("NbinsMHT");
   NbinsHT_ = iConfig.getParameter<int> ("NbinsHT");
   Ntries_ = iConfig.getParameter<int> ("Ntries");
   NJetsSave_ = iConfig.getParameter<int> ("NJetsSave");
   HTSave_ = iConfig.getParameter<double> ("HTSave");
   MHTSave_ = iConfig.getParameter<double> ("MHTSave");
   JetsHTPt_ = iConfig.getParameter<double> ("JetsHTPt");
   JetsHTEta_ = iConfig.getParameter<double> ("JetsHTEta");
   JetsMHTPt_ = iConfig.getParameter<double> ("JetsMHTPt");
   JetsMHTEta_ = iConfig.getParameter<double> ("JetsMHTEta");
   JetDeltaMin_ = iConfig.getParameter<std::vector<double> > ("JetDeltaMin");
   MHTcut_low_ = iConfig.getParameter<double> ("MHTcut_low");
   MHTcut_medium_ = iConfig.getParameter<double> ("MHTcut_medium");
   MHTcut_high_ = iConfig.getParameter<double> ("MHTcut_high");
   HTcut_low_ = iConfig.getParameter<double> ("HTcut_low");
   HTcut_medium_ = iConfig.getParameter<double> ("HTcut_medium");
   HTcut_high_ = iConfig.getParameter<double> ("HTcut_high");
   HTcut_veryhigh_ = iConfig.getParameter<double> ("HTcut_veryhigh");
   HTcut_extremehigh_ = iConfig.getParameter<double> ("HTcut_extremehigh");

   PtBinEdges_ = iConfig.getParameter<std::vector<double> > ("PtBinEdges");
   EtaBinEdges_ = iConfig.getParameter<std::vector<double> > ("EtaBinEdges");

   unsigned int needed_dim = (PtBinEdges_scaling_.size() - 1) * (EtaBinEdges_scaling_.size() - 1);
   if (AdditionalSmearing_.size() != needed_dim) {
      throw edm::Exception(edm::errors::Configuration, "AdditionalSmearing has not correct dimension");
   }
   if (LowerTailScaling_.size() != needed_dim) {
      throw edm::Exception(edm::errors::Configuration, "LowerTailScaling has not correct dimension");
   }
   if (UpperTailScaling_.size() != needed_dim) {
      throw edm::Exception(edm::errors::Configuration, "UpperTailScaling has not correct dimension");
   }

   // Different seed per initialization
   gRandom->SetSeed(0);
   rand_ = new TRandom3(0);

   // get object of class SmearFunction
   smearFunc_ = new SmearFunction(iConfig);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
QCDBkgRS::~QCDBkgRS()
{
   PtBinEdges_.clear();
   EtaBinEdges_.clear();

   delete smearFunc_;

   if (rand_)
      delete rand_;
}
//--------------------------------------------------------------------------

//
// member functions
//

//--------------------------------------------------------------------------
std::string QCDBkgRS::GetName(const std::string plot, const std::string uncert, const std::string ptbin) const {
   std::string result(plot);
   if (uncert != "" && uncert != " ")
      result += "_" + uncert;
   if (ptbin != "" && ptbin != " ")
      result += "_" + ptbin;
   return result;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
int QCDBkgRS::GetIndex(const double& x, const std::vector<double>* vec) {
   int index = -1;
   // this is a check
   //int index = 0;
   for (std::vector<double>::const_iterator it = vec->begin(); it != vec->end(); ++it) {
      if ((*it) > fabs(x))
         break;
      ++index;
   }
   if (index < 0)
      index = 0;
   if (index > (int) vec->size() - 2)
      index = vec->size() - 2;

   return index;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// pt resolution for KinFitter
double QCDBkgRS::JetResolution_Pt2(const double& pt, const double& eta, const int& i) {
   int i_jet;
   i < 1 ? i_jet = i : i_jet = 0;
   int i_eta = GetIndex(eta, &EtaBinEdges_);
   //return pow(pt, 2) * pow(smearFunc_->getSigmaPtScaledForRebalancing(i_jet, i_eta)->Eval(pt), 2);
   return pow(pt, 2) * pow(smearFunc_->getSigmaPtForRebalancing(i_jet, i_eta)->Eval(pt), 2);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// relative pt resolution for KinFitter
double QCDBkgRS::JetResolution_Ptrel(const double& pt, const double& eta, const int& i) {
   int i_jet;
   i < 1 ? i_jet = i : i_jet = 0;
   int i_eta = GetIndex(eta, &EtaBinEdges_);
   return smearFunc_->getSigmaPtScaledForRebalancing(i_jet, i_eta)->Eval(pt);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// eta resolution for KinFitter
double QCDBkgRS::JetResolution_Eta(const double& pt, const double& eta, const int& i, const int& i_flav) {
   int i_jet;
   i < 1 ? i_jet = i : i_jet = 0;
   int i_eta = GetIndex(eta, &EtaBinEdges_);
   int i_Pt = GetIndex(pt, &PtBinEdges_);
   return smearFunc_->SigmaEtaHist.at(i_flav).at(i_jet).at(i_eta).at(i_Pt);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// phi resolution for KinFitter
double QCDBkgRS::JetResolution_Phi(const double& pt, const double& eta, const int& i, const int& i_flav) {
   int i_jet;
   i < 1 ? i_jet = i : i_jet = 0;
   int i_eta = GetIndex(eta, &EtaBinEdges_);
   int i_Pt = GetIndex(pt, &PtBinEdges_);
   return smearFunc_->SigmaPhiHist.at(i_flav).at(i_jet).at(i_eta).at(i_Pt);
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// pt resolution for smearing
double QCDBkgRS::JetResolutionHist_Pt_Smear(const double& pt, const double& eta, const int& i, const double& HT, const int& NJets, const bool btag) {
   int i_jet;
   i < 1 ? i_jet = i : i_jet = 0;
   int i_Pt = GetIndex(pt, &PtBinEdges_);
   int i_eta = GetIndex(eta, &EtaBinEdges_);

   double res = 1.0;
   if( btag ){
      // get heavy flavor smear function
      res = smearFunc_->getSmearFunc(1, i_jet, i_eta, i_Pt)->GetRandom();
      //cout << "B-jet: response " << res << endl;
   }
   else {
      // get no heavy flavor smear function
      res = smearFunc_->getSmearFunc(0, i_jet, i_eta, i_Pt)->GetRandom();
      //cout << "noB-jet: response " << res << endl;
   }

   return res;
}
//--------------------------------------------------------------------------
// Get Rebalance corrections:
//--------------------------------------------------------------------------
double QCDBkgRS::GetRebalanceCorrection(double jet_pt, bool btag)
{
   if( jet_pt > 1000. ) jet_pt = 999.;

   if ( jet_pt > rebalancedJetPt_ ) {
      int i_bin = h_RebCorrectionFactor->FindBin(jet_pt);

      if (!useCleverRebalanceCorrectionFactors_){
         double result = 0;
         if (btag) {
            result = h_RebCorrectionFactor_b->GetBinContent(i_bin);
         } else {
            result = h_RebCorrectionFactor->GetBinContent(i_bin);
         }
         if (result < 0.01) result = 1.;
         return result;
      } else {
         double result = 0;
         if (btag) {
            result = h_2DRebCorrectionFactor_b_py.at(i_bin-1)->GetRandom();;
         } else {
            result = h_2DRebCorrectionFactor_py.at(i_bin-1)->GetRandom();;
         }
         if (result < 0.01) result = 1.;
         return result;
      }
   }

   else return 1.;

}
//--------------------------------------------------------------------------

double QCDBkgRS::ReadXYHist(TH1F* hist, double x)
{
   // Given x, return y for a 2d histogram.
   //cout << "Called ReadXYHist..." << endl;
   int the_bin = hist->GetXaxis()->FindBin(x);
   if ( hist->IsBinOverflow(the_bin) || hist->IsBinUnderflow(the_bin) ){
      // overflow/underflow can be empty
      cout << "Warning: Reading from under-/overflow; x=" << x << "bin=" << endl;
      // Not optimal; I would like to return the python equiv of None.
      // How to do that?
      return 0;
   }
   double y = hist->GetBinContent(the_bin);
   return y;
}


//--------------------------------------------------------------------------
// Get BTag Efficiencies:
//--------------------------------------------------------------------------
double QCDBkgRS::GetBTagEfficiency(double pt, double eta)
{
    //cout << "Called GetBTagEfficiency..." << endl;
    // Set to unrealistic value:
    double btageff = 0.;
    int eta_bin = GetIndex(eta, &EtaBinEdges_);
    h_BTagEfficiencyFactor = BTagEfficiencyFactors[eta_bin];
    //btageff = h_BTagEfficiencyFactor->GetBinContent(h_BTagEfficiencyFactor->GetXaxis()->FindBin(pt));
    btageff = ReadXYHist(h_BTagEfficiencyFactor, pt);
    //cout << "Found efficiency " << btageff << " for pt " << pt << " and eta " << eta << endl;

    if( btageff == 0.0 ){
        // Reset to unrealistic value.
        btageff = 0.;
    }

    return btageff;

}
//--------------------------------------------------------------------------
// Get BMisTag Efficiencies:
//--------------------------------------------------------------------------
double QCDBkgRS::GetBMisTagEfficiency(double pt, double eta)
{
    //cout << "Called GetBmistagEfficiency..." << endl;
    double bmistageff = 1.;
    int eta_bin = GetIndex(eta, &EtaBinEdges_);
    h_BMisTagEfficiencyFactor = BMisTagEfficiencyFactors[eta_bin];
    //bmistageff = h_BMisTagEfficiencyFactor->GetBinContent(h_BMisTagEfficiencyFactor->GetXaxis()->FindBin(pt));
    bmistageff = ReadXYHist(h_BMisTagEfficiencyFactor, pt);
    //cout << "Found mistag efficiency " << bmistageff << " for pt " << pt << " and eta " << eta << endl;

    return bmistageff;
}

//
//--------------------------------------------------------------------------
// Get NBTag:
//--------------------------------------------------------------------------
double QCDBkgRS::GetNBTag(double pt, double eta)
{
    //cout << "Called GetNBTag..." << endl;
    double nbtag = 0.;
    int eta_bin = GetIndex(eta, &EtaBinEdges_);
    h_NBTag = NBTags[eta_bin];
    //btageff = h_BTagEfficiencyFactor->GetBinContent(h_BTagEfficiencyFactor->GetXaxis()->FindBin(pt));
    nbtag = ReadXYHist(h_NBTag, pt);
    //cout << "Found NBTag " << nbtag << " for pt " << pt << " and eta " << eta << endl;


    return nbtag;

}

//
//--------------------------------------------------------------------------
// Get NnoBTag:
//--------------------------------------------------------------------------
double QCDBkgRS::GetNnoBTag(double pt, double eta)
{
    //cout << "Called GetNnoBTag..." << endl;
    double n_no_btag = 0.;
    int eta_bin = GetIndex(eta, &EtaBinEdges_);
    h_NnoBTag = NnoBTags[eta_bin];
    //btageff = h_BTagEfficiencyFactor->GetBinContent(h_BTagEfficiencyFactor->GetXaxis()->FindBin(pt));
    n_no_btag = ReadXYHist(h_NnoBTag, pt);
    //cout << "Found NnoBTag " << n_no_btag << " for pt " << pt << " and eta " << eta << endl;


    return n_no_btag;

}

//
//--------------------------------------------------------------------------
// Get NBTrue:
//--------------------------------------------------------------------------
double QCDBkgRS::GetNBTrue(double pt, double eta)
{
    //cout << "Called GetNBTrue..." << endl;
    double nbtrue = 0.;
    int eta_bin = GetIndex(eta, &EtaBinEdges_);
    h_NBTrue = NBTrues[eta_bin];
    //btageff = h_BTagEfficiencyFactor->GetBinContent(h_BTagEfficiencyFactor->GetXaxis()->FindBin(pt));
    nbtrue = ReadXYHist(h_NBTrue, pt);
    //cout << "Found NBTrue ratio " << nbtrue << " for pt " << pt << " and eta " << eta << endl;


    return nbtrue;

}



//--------------------------------------------------------------------------
// rebalance the events using a kinematic fit and transverse momentum balance
bool QCDBkgRS::RebalanceJets_KinFitter(edm::View<pat::Jet>* Jets_rec, std::vector<pat::Jet> &Jets_reb) {

   bool result = true;

   //// Interface to KinFitter
   TKinFitter* myFit = new TKinFitter();

   std::vector<TLorentzVector*> lvec_m;

   std::vector<TMatrixD*> covMat_m;

   std::vector<TFitParticleEtEtaPhi*> fitted;
   std::vector<TFitParticleEtEtaPhi*> measured;
   std::map<int, const pat::Jet*> JetMap;
   double dPx = 0;
   double dPy = 0;
   double HTreco = 0;
   double HTreb = 0;
   double MHTx_low = 0;
   double MHTy_low = 0;

   //// Fill measured particles to vector
   int i = 0;
   for (edm::View<pat::Jet>::const_iterator it = Jets_rec-> begin(); it != Jets_rec->end(); ++it) {

      //if (it->pt() < rebalancedJetPt_ || abs(it->pt()) > 3.0) {
      //if (it->pt() < rebalancedJetPt_ || it->chargedEmEnergyFraction()>0.9 || it->muonEnergyFraction()>0.9) {
      if (it->pt() < rebalancedJetPt_) {
         if (rebalanceMode_ == "MHTall") {
            MHTx_low -= it->px();
            MHTy_low -= it->py();
            pat::Jet rebalancedJet(*((pat::Jet*) &(*it)));
            Jets_reb.push_back(rebalancedJet);
         }
      } else {

         JetMap[i] = &(*it);

         // The particles before fitting
         double tmppx, tmppy, tmppz, tmpe;

         if( useRebalanceCorrectionFactors_ ) {
            bool btag = (it->bDiscriminator(btagTag_) > btagCut_);
            tmppx = it->px()/GetRebalanceCorrection( it->pt(), btag );
            tmppy = it->py()/GetRebalanceCorrection( it->pt(), btag );
            tmppz = it->pz()/GetRebalanceCorrection( it->pt(), btag );
            tmpe = it->energy()/GetRebalanceCorrection( it->pt(), btag );
         }
         else {
            tmppx = it->px();
            tmppy = it->py();
            tmppz = it->pz();
            tmpe = it->energy();
         }

         TLorentzVector* lv = new TLorentzVector(tmppx, tmppy, tmppz, tmpe);
         lvec_m.push_back(lv);
         TMatrixD* cM = new TMatrixD(3, 3);
         (*cM)(0, 0) = JetResolution_Pt2(it->pt(), it->eta(), i);
         (*cM)(1, 1) = pow(JetResolution_Eta(it->pt(), it->eta(), 0, 0), 2);
         (*cM)(2, 2) = pow(JetResolution_Phi(it->pt(), it->eta(), 0, 0), 2);
         covMat_m.push_back(cM);
         char name[10];
         sprintf(name, "jet%i", i);
         TFitParticleEtEtaPhi* jet1 = new TFitParticleEtEtaPhi(name, name, lvec_m.back(), covMat_m.back());
         measured.push_back(jet1);
         TFitParticleEtEtaPhi* jet2 = new TFitParticleEtEtaPhi(name, name, lvec_m.back(), covMat_m.back());
         fitted.push_back(jet2);
         myFit->addMeasParticle(fitted.back());
         ++i;
      }
   }

   //// Add momentum constraints
   double MET_constraint_x = 0.;
   double MET_constraint_y = 0.;
   if (rebalanceMode_ == "MHTall") {
      //// rebalance MHT of all jets
      MET_constraint_x = MHTx_low;
      MET_constraint_y = MHTy_low;
   } else if (rebalanceMode_ == "MHThigh") {
      //// rebalance MHT of fitted jets
      MET_constraint_x = 0.;
      MET_constraint_y = 0.;
   } else {
      //// default: rebalance MHT of fitted jets
      MET_constraint_x = 0.;
      MET_constraint_y = 0.;
   }
   TFitConstraintEp* momentumConstr1 = new TFitConstraintEp("px", "px", 0, TFitConstraintEp::pX, MET_constraint_x);
   TFitConstraintEp* momentumConstr2 = new TFitConstraintEp("py", "py", 0, TFitConstraintEp::pY, MET_constraint_y);
   for (unsigned int i = 0; i < fitted.size(); ++i) {
      momentumConstr1->addParticle(fitted.at(i));
      momentumConstr2->addParticle(fitted.at(i));
   }
   myFit->addConstraint(momentumConstr1);
   myFit->addConstraint(momentumConstr2);

   //// Set fit parameters
   myFit->setVerbosity(0);
   myFit->setMaxNbIter(100);
   myFit->setMaxF(0.01 * 2);
   myFit->setMaxDeltaS(1.e-3);
   myFit->fit();

   //cout << "KinFitter: " << myFit->getStatus() << endl;
   int status = myFit->getStatus();

   double chi2 = 0;
   //double F = 0;
   double prob = 0;
   //if (status == 0 || status == 1) {
   if (status == 0) {
      chi2 = myFit->getS();
      //F = myFit->getF();
      int dof = myFit->getNDF();
      prob = TMath::Prob(chi2, dof);
      //if (prob < 1.e-8) result = false;
      //    if (prob < 0.01) result = false;
      //cout << "chi2, prop, F = " << chi2 << " " << prob << " " << F << endl;
   } else {
      chi2 = 99999;
      prob = 0;
      //F = 99999;
      result = false;
      //cout << "chi2, prop, F = " << chi2 << " " << prob << " " << F << endl;
   }
   if (controlPlots_)
      h_fitProb->Fill(prob);

   //// Get the output of KinFitter
   for (unsigned int i = 0; i < measured.size(); ++i) {
      // create new rebalanced Jet
      pat::Jet::LorentzVector newP4(fitted.at(i)->getCurr4Vec()->Px(), fitted.at(i)->getCurr4Vec()->Py(),
                                    fitted.at(i)->getCurr4Vec()->Pz(), fitted.at(i)->getCurr4Vec()->E());
      pat::Jet rebalancedJet(*((pat::Jet*) JetMap[i]));
      HTreco += rebalancedJet.pt();
      //cout << "RECO: " << i << "th: pt = " << rebalancedJet.pt() << " phi = " << rebalancedJet.phi() << endl;
      rebalancedJet.setP4(newP4);
      //cout << "REB: " << i << "th: pt = " << rebalancedJet.pt() << " phi = " << rebalancedJet.phi() << endl;
      Jets_reb.push_back(rebalancedJet);
      dPx -= newP4.Px() - measured.at(i)->getCurr4Vec()->Px();
      dPy -= newP4.Py() - measured.at(i)->getCurr4Vec()->Py();
      HTreb += rebalancedJet.pt();
   }
   //cout << "HT reco = " << HTreco << endl;
   //cout << "HT reb  = " << HTreb << endl;

   delete myFit;
   for (unsigned int i = 0; i < measured.size(); ++i) {
      delete lvec_m.at(i);
      delete covMat_m.at(i);
      delete measured.at(i);
      delete fitted.at(i);
   }
   delete momentumConstr1;
   delete momentumConstr2;

   return result;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void QCDBkgRS::SmearingJets(const std::vector<pat::Jet> &Jets_reb, std::vector<pat::Jet> &Jets_smeared) {

   double dPx = 0;
   double dPy = 0;
   double HT = calcHT(Jets_reb);
   int NJets_reb = calcNJets(Jets_reb);

   TRandom3 *r;
   // Make new map where we will store our own decision of
   // whether a jet is b-tagged (rather than take the given
   // b-tag information) which will then later only be used
   // to count the number of b-jets in an event:
   std::map <const reco::Jet*, bool> dynamic_jet_btag_map;
   for (int i = 1; i <= Ntries_; ++i) {
      int Ntries2 = 1;
      double w = weight_;
      if (cleverPrescaleTreating_ == true && weight_ > 1) {
         Ntries2 = (int) weight_;
         if (Ntries2 > 100) Ntries2 = 100;
         w = weight_ / Ntries2;
      }
      for (int j = 1; j <= Ntries2; ++j) {
         double btag_correction = 1.;
         double p_btrue = 0.;
         double p_btag = 0.;
         Jets_smeared.clear();
         int i_jet = 0;

         // Clear map:
         dynamic_jet_btag_map.clear();

         // Iterate over jets:
         for (std::vector<pat::Jet>::const_iterator it = Jets_reb.begin(); it != Jets_reb.end(); ++it) {

            if (it->pt() > smearedJetPt_) {
               bool btag = (it->bDiscriminator(btagTag_) > btagCut_);

               // Smearing
               //double scale = JetResolutionHist_Pt_Smear(it->pt(), it->eta(), i_jet, HT, NJets_reb, btag);

               //-------------------------------------------------------
               // Compute the probability that a jet was a true b-jet:
               //-------------------------------------------------------
               // Here you actually want the reco pt, not the rebalanced pt:
               double BTagEff = GetBTagEfficiency(it->pt(), it->eta());
               double NBTrue = GetNBTrue(it->pt(), it->eta());
               double NBTag = GetNBTag(it->pt(), it->eta());
               double NnoBTag = GetNnoBTag(it->pt(), it->eta());
               //cout << "Found BTagEff " << BTagEff << endl;
               //cout << "Found NBTrue " << NBTrue << endl;
               //cout << "Found NBTag " << NBTag << endl;
               //cout << "Found NnoBTag " << NnoBTag << endl;
               if(btag){
                   p_btrue =  BTagEff * NBTrue/NBTag;
               } else {
                   p_btrue = (1 - BTagEff) * NBTrue/NnoBTag;
               }
               //cout << "P(BTrue): " << p_btrue << endl;
               // Pick a random number between 0 and 1:

               r = new TRandom3(0); // generate a number in interval ]0,1] (0 is excluded)
               double random_number = r->Rndm();
               delete r;

               // The particle is a true b with probability p_btrue;
               // decide with the random number:
               bool btrue = p_btrue > random_number;
               int i_flav = 0;
               //if (btag){
               // Decide which jet resolution template to use based on whether the
               // particle was a true b (not based on its b-tag):
               if (btrue){
                  i_flav = 1;
                  //cout << "There was a true b jet" << endl;
               }
               double scale = JetResolutionHist_Pt_Smear(it->pt(), it->eta(), i_jet, HT, NJets_reb, btrue);
               //double scale = JetResolutionHist_Pt_Smear(it->pt(), it->eta(), i_jet, HT, NJets_reb, i_flav);
               double newE = it->energy() * scale;
               double newMass = it->mass() * scale;
               double newEta = rand_->Gaus(it->eta(), JetResolution_Eta(it->pt(), it->eta(), i_jet, i_flav));
               double newPhi = rand_->Gaus(it->phi(), JetResolution_Phi(it->pt(), it->eta(), i_jet, i_flav));
               double newPt = sqrt(newE*newE-newMass*newMass)/cosh(newEta);



               // No smearing:
               //double newEta = it->eta();
               //double newPhi = it->phi();
               //
               pat::Jet::PolarLorentzVector newP4(newPt, newEta, newPhi, it->mass());
               pat::Jet smearedJet(*it);
               smearedJet.setP4(newP4);
               Jets_smeared.push_back(smearedJet);
               dPx -= newP4.Px() - it->px();
               dPy -= newP4.Py() - it->py();


               //-------------------------------------------------------
               // Compute the probability that a jet was b-tagged:
               //-------------------------------------------------------
               double newBTagEff = GetBTagEfficiency(newPt, newEta);
               double newBMisTagEff = GetBMisTagEfficiency(newPt, newEta);
               //cout << "Found NewBTagEff " << newBTagEff << endl;
               //cout << "Found NewBMisTagEff " << newBMisTagEff << endl;
               // Btag correction factors
               // This is deprecated
               double oldBTagEff = BTagEff;
               //cout << "NewBtageff: "<< newBTagEff << " old btageff: " << oldBTagEff << " correction: " << newBTagEff/oldBTagEff << endl;
               // Take care not to divide by zero:
               if(btag && (oldBTagEff > 0)){

                    // BTag efficiency should be less than 1!
                    // If btag efficiency is 1, this is because
                    // it could not be read from the file. Ignore
                    // and keep correction factor 1:
                    if (newBTagEff < 1.0 && oldBTagEff < 1.0){
                        //cout << "NewBtageff: "<< newBTagEff << "old btageff: " << oldBTagEff << " correction: " << btag_correction << endl;
                        btag_correction *= newBTagEff/oldBTagEff;

                    }
               }
               //cout << "Btag correction factor:" << btag_correction << endl;

               // Compute probability that the particle is b-tagged.
               // Ignore the given b-tag information here:
               p_btag = p_btrue * newBTagEff + (1 - p_btrue) * newBMisTagEff;
               //cout << "p(BTag): " << p_btag << endl;
               r = new TRandom3(0); // generate a number in interval ]0,1] (0 is excluded)
               double random_btag = r->Rndm();
               delete r;
               bool btagged = p_btag > random_btag;
               if(btagged) {
                   //cout << "There was a b-tagged jet" << endl;
               }

               //if(p_btag > 1){
               //    cout << " " << endl;
               //    cout << "Found p(BTag) > 1!" << endl;
               //    cout << "p(BTag): " << p_btag << endl;
               //    cout << "p(BTrue): " << p_btrue << endl;
               //    cout << "eff(BTag,new): " << newBTagEff << endl;
               //    cout << "eff(BMisTag,new): " << newBMisTagEff << endl;
               //    cout << "Random number: " << random_btag << endl;
               //    cout << "Considered b-tagged:" << btagged << endl;
               //    cout << " " << endl;
               //}

               // Add btag to map.
               // This is the information used to determine the number
               // of b's in an event.
               dynamic_jet_btag_map[&(Jets_smeared.back())] = btagged;
               //-------------------------------------------------------

               h_pBTrue_smear->Fill(p_btrue);
               h_pBTag_smear->Fill(p_btag);
               h_random_BTrue_smear->Fill(random_number);
               h_random_BTag_smear->Fill(random_btag);
               ++i_jet;
            } else {
               pat::Jet smearedJet(*it);
               Jets_smeared.push_back(smearedJet);
            }
         }
         h_BTagCorrectionFactor_smear->Fill(btag_correction);
         GreaterByPt<reco::Candidate> ptComparator_;
         std::sort(Jets_smeared.begin(), Jets_smeared.end(), ptComparator_);

         //Fill HT and MHT prediction histos for i-th iteration of smearing
         int NJets = calcNJets(Jets_smeared);
         if (NJets >= NJetsSave_) {
            //FillPredictions(Jets_smeared, i, w*btag_correction);
            FillPredictions(Jets_smeared, i, w, dynamic_jet_btag_map);

            if( HT_pred > HTSave_ && MHT_pred > MHTSave_){
               PredictionTree->Fill();
            }

            // clean variables in tree
            weight = 0.;
            Ntries_pred = 0.;
            Njets_pred = 0;
            BTags_pred = 0;
            HT_pred = 0.;
            MHT_pred = 0.;
            Jet1Pt_pred = 0.;
            Jet2Pt_pred = 0.;
            Jet3Pt_pred = 0.;
            Jet4Pt_pred = 0.;
            Jet1Eta_pred = 0.;
            Jet2Eta_pred = 0.;
            Jet3Eta_pred = 0.;
            Jet4Eta_pred = 0.;
            DeltaPhi1_pred = 0.;
            DeltaPhi2_pred = 0.;
            DeltaPhi3_pred = 0.;
            DeltaPhiMinN_pred = 0.;
         }
      }
   }

   return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void QCDBkgRS::SmearingGenJets(edm::View<reco::GenJet>* Jets_gen, edm::View<pat::PackedGenParticle>* genParticles, std::vector<reco::GenJet> &GenJets_smeared) {

   double dPx = 0;
   double dPy = 0;

   // calculate quantities needed for smearing
   double HT;
   int NJets_gen = 0;
   for (edm::View<reco::GenJet>::const_iterator it = Jets_gen->begin(); it != Jets_gen->end(); ++it) {
      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         ++NJets_gen;
         HT += it->pt();
      }
   }

   TRandom3 *r;
   // Make new map where we will store our own decision of
   // whether a jet is b-tagged (rather than take the given
   // b-tag information) which will then later only be used
   // to count the number of b-jets in an event:
   std::map <const reco::GenJet*, bool> dynamic_genjet_btag_map;
   for (int i = 1; i <= Ntries_; ++i) {
      int Ntries2 = 1;
      double w = weight_;
      if (cleverPrescaleTreating_ == true && weight_ > 1) {
         Ntries2 = (int) weight_;
         if (Ntries2 > 100) Ntries2 = 100;
         w = weight_ / Ntries2;
      }

      for (int j = 1; j <= Ntries2; ++j) {

         double btag_correction = 1.;

         dynamic_genjet_btag_map.clear();
         GenJets_smeared.clear();
         std::map <const reco::GenJet*, bool> genJet2_btag;

         int i_jet = 0;

         for (edm::View<reco::GenJet>::const_iterator it = Jets_gen->begin(); it != Jets_gen->end(); ++it) {

            //// Find additional neutrinos and add them back since these are not included in the genJets
            math::PtEtaPhiMLorentzVector neutrinos(0., 0., 0., 0.);
            for(edm::View<pat::PackedGenParticle>::const_iterator cand = genParticles->begin(); cand!=genParticles->end(); ++cand)
            {
               if ( cand->status()==1 && (abs(cand->pdgId())==12 || abs(cand->pdgId())==14 || abs(cand->pdgId())==16)){
                  double dR = deltaR(*it, *cand);
                  if (dR < 0.4) neutrinos += cand->p4();
               }
            }

            if ((it->p4()+neutrinos).pt() > smearedJetPt_) {

               // Remember that all jets are matched exactly.
               // genJet_btag tells whether the matched jet is b-tagged:
               // Originally, this information was used to choose
               // between b and light-flavor templates.
               bool btag = genJet_btag[&(*it)];


               //-------------------------------------------------------
               // Compute the probability that a jet was a true b-jet:
               //-------------------------------------------------------
               double BTagEff = GetBTagEfficiency(it->pt(), it->eta());
               double NBTrue = GetNBTrue(it->pt(), it->eta());
               double p_btrue = 0.;
               double p_btag = 0.;
               if(btag){

                   double NBTag = GetNBTag(it->pt(), it->eta());
                   p_btrue =  BTagEff * NBTrue/NBTag;
               } else {
                   double NnoBTag = GetNnoBTag(it->pt(), it->eta());
                   p_btrue = (1 - BTagEff) * NBTrue/NnoBTag;
               }
               //cout << "P(BTrue) (gen): " << p_btrue << endl;
               // Pick a random number between 0 and 1:
               r = new TRandom3(0); // generate a number in interval ]0,1] (0 is excluded)
               double random_number = r->Rndm();
               // The particle is a true b with probability p_btrue;
               // decide with the random number:

               // TEST: Take btrue to be the real btrue:
               bool btrue_random = p_btrue > random_number; // If the random number is < pbtrue, it is a true b
               //
               //
               //
               bool btrue = false;
                 //--------- Decide on btrue based on b in genlevel constituents ---------

            //   // From GenJet::getGenConstituents
            //   // and from GenJet::getGenConstituent:
            //   // Decide on btrue based on whether there was a b inside:
            //   for (reco::Candidate::const_iterator daugh = it->begin(); daugh != it->end (); daugh++) {
            //        if ( (&*daugh != NULL) && (daugh != it->end ())) { // in range

            //            const reco::Candidate* constituent = &*daugh; // deref
            //            //double constituent_pt = constituent->pt();
            //            int constituent_pid = std::abs(constituent->pdgId());
            //            constituent_pid = constituent_pid % 10000;

            //            // Check if pid is that of a b meson:
            //            if (((500 <= constituent_pid) && (constituent_pid < 600)) || ((5000 <= constituent_pid) && (constituent_pid < 6000))){
            //                double conesize = 0.5;
            //                // Check if it was inside a cone of 0.5:
            //                // (Should you check for pt?)
            //                if((deltaR(*it, *constituent) < conesize)){
            //                    btrue = true;
            //                }
            //            }
            //        }
            //   }

               // Decide btrue based on probability (see above):
               btrue = btrue_random;

               // For example: pbtrue = 0.1 (10%): then the random number is less likely to be
               // less than 0.1, more likely to be greater than 0.1 -> true b if less than 0.1.
               delete r;
               int i_flav = 0;
               //if (btag) i_flav = 1;
               // Decide which jet resolution template to use based on whether the
               // particle was a true b (not based on its b-tag):
               if (btrue){
                  i_flav = 1;
               }
               //-------------------------------------------------------

               double invScale = (it->p4()+neutrinos).pt() / it->pt();
               //double scale = JetResolutionHist_Pt_Smear(it->pt(), it->eta(), i_jet, HT, NJets_gen, btag);
               // Use templates based on whether it was a true b-jet:
               double scale = JetResolutionHist_Pt_Smear(it->pt(), it->eta(), i_jet, HT, NJets_gen, btrue);
               double newE = it->energy() * scale * invScale;
               double newMass = it->mass() * scale;
               double newEta = rand_->Gaus(it->eta(), JetResolution_Eta(it->pt(), it->eta(), i_jet, i_flav));
               double newPhi = rand_->Gaus(it->phi(), JetResolution_Phi(it->pt(), it->eta(), i_jet, i_flav));
               double newPt = sqrt(newE*newE-newMass*newMass)/cosh(newEta);

               // No smearing:
               //double newEta = it->eta();
               //double newPhi = it->phi();

               reco::GenJet::PolarLorentzVector newP4(newPt, newEta, newPhi, it->mass());
               reco::GenJet smearedJet(*it);
               smearedJet.setP4(newP4);
               GenJets_smeared.push_back(smearedJet);
               dPx -= newP4.Px() - it->px();
               dPy -= newP4.Py() - it->py();
               genJet2_btag[&(GenJets_smeared.back())] = btag;

               //-------------------------------------------------------
               // Compute the probability that a jet was b-tagged:
               //-------------------------------------------------------
               double newBTagEff = GetBTagEfficiency(newPt, newEta);
               double newBMisTagEff = GetBMisTagEfficiency(newPt, newEta);
               // Btag correction factors
               // This is deprecated
               double oldBTagEff = GetBTagEfficiency(it->pt(), it->eta());
               double oldBMisTagEff = GetBMisTagEfficiency(it->pt(), it->eta());
               if(btag){
                    // BTag efficiency should be less than 1!
                    // If btag efficiency is 1, this is because
                    // it could not be read from the file. Ignore
                    // and keep correction factor 1:
                    if (newBTagEff < 1.0 && oldBTagEff < 1.0){
                        btag_correction *= newBTagEff/oldBTagEff;
                    }
               }

               // Compute probability that the particle is b-tagged.
               // Ignore the given b-tag information here:
               p_btag = p_btrue * newBTagEff + (1 - p_btrue) * newBMisTagEff;
               //cout << "P(Btag) (gen): " << p_btag << endl;
               r = new TRandom3(0); // generate a number in interval ]0,1] (0 is excluded)
               double random_btag = r->Rndm();
               bool btagged = p_btag > random_btag; // If the random number is < pbtag, it is btagged
               // For example: pbtag = 0.1 (10%): then the random number is less likely to be
               // less than 0.1, more likely to be greater than 0.1 -> btagged if less than 0.1.
               delete r;
               //if(btagged) {
               //    cout << "There was a b-tagged genjet" << endl;
               //}
               dynamic_genjet_btag_map[&(GenJets_smeared.back())] = btagged;
               //cout << "This jet has " << ( btagged ? "been" : "not been") << " b-tagged" << endl;
               //cout << "Originally this event was " << (btag ? "b-tagged" : "not b-tagged") <<endl;
               //dynamic_genjet_btag_map[*it] = btagged;
               //-------------------------------------------------------

               h_pBTrue_smear->Fill(p_btrue);
               h_pBTag_smear->Fill(p_btag);
               h_random_BTrue_smear->Fill(random_number);
               h_random_BTag_smear->Fill(random_btag);

               //int i_eta = GetIndex(newEta, &EtaBinEdges);
               //if(btrue){
               //     h_trueb_SmearedPt.at(i_eta)->Fill(newPt, w);
               //     if(btagged){
               //         h_trueb_btag_SmearedPt.at(i_eta)->Fill(newPt, w);
               //     }
               //} else {
               //     h_no_trueb_SmearedPt.at(i_eta)->Fill(newPt, w);
               //     if(btagged){
               //         h_no_trueb_btag_SmearedPt.at(i_eta)->Fill(newPt, w);
               //     }
               //}
               //if(btagged){
               //     h_btag_SmearedPt.at(i_eta)->Fill(newPt, w);
               //}
               //if(btag){
               //     h_btag_old_SmearedPt.at(i_eta)->Fill(newPt, w);
               //}

               h_btageff_SmearPt_Response->Fill(newPt, scale, w*newBTagEff);
               h_bmistageff_SmearPt_Response->Fill(newPt, scale, w*newBMisTagEff);
               h_btageff_genPt_Response->Fill(it->pt(), scale, w*oldBTagEff);
               h_bmistageff_genPt_Response->Fill(it->pt(), scale, w*oldBMisTagEff);
               double oldPt = it->pt();
               double oldEta = it->eta();

               // Save TH2F instead:
               if(btrue_random){
                    h_trueb_random_SmearedPt_Eta->Fill(newPt, newEta, w);
                    h_trueb_genPt_Eta->Fill(oldPt, oldEta, w);
               }
               if(btrue){
                    h_trueb_SmearedPt_Eta->Fill(newPt, newEta, w);
                    if(btagged){
                        h_trueb_btag_SmearedPt_Eta->Fill(newPt, newEta, w);
                        h_trueb_btag_genPt_Eta->Fill(oldPt, oldEta, w);
                    }
               } else {
                    h_no_trueb_SmearedPt_Eta->Fill(newPt, newEta, w);
                    h_no_trueb_genPt_Eta->Fill(oldPt, oldEta, w);
                    if(btagged){
                        h_no_trueb_btag_SmearedPt_Eta->Fill(newPt, newEta, w);
                        h_no_trueb_btag_genPt_Eta->Fill(oldPt, oldEta, w);
                    }
               }
               if(btagged){
                    h_btag_SmearedPt_Eta->Fill(newPt, newEta, w);
               }
               if(btag){
                    h_btag_old_SmearedPt_Eta->Fill(newPt, newEta, w);
               }

               ++i_jet;
            } else {
               reco::GenJet smearedJet(*it);
               GenJets_smeared.push_back(smearedJet);
            }
         }
         h_BTagCorrectionFactor_smear->Fill(btag_correction);
         GreaterByPt<reco::Candidate> ptComparator_;
         std::sort(GenJets_smeared.begin(), GenJets_smeared.end(), ptComparator_);

         //Fill HT and MHT prediction histos for i-th iteration of smearing
         int NJets = calcNJets_gen(GenJets_smeared);
         if (NJets >= NJetsSave_) {
            //FillPredictions_gen(GenJets_smeared, i, w*btag_correction, genJet2_btag);
            // for each iteration, use a newly computed map genjet <-> btag.
            FillPredictions_gen(GenJets_smeared, i, w, dynamic_genjet_btag_map);
            int NB_Old = calcNBJets_gen(GenJets_smeared, genJet2_btag);
            h_NB_Old->Fill(NB_Old, w);
            //cout << "Number of old btagged jets: " << NB_Old << endl;
            int NB_New = calcNBJets_gen(GenJets_smeared, dynamic_genjet_btag_map);
            h_NB_New->Fill(NB_New, w);
            //cout << "Number of new btagged jets: " << NB_New << endl;

            if( HT_pred > HTSave_ && MHT_pred > MHTSave_){
               PredictionTree->Fill();
            }

            // clean variables in tree
            weight = 0.;
            Ntries_pred = 0.;
            Njets_pred = 0;
            BTags_pred = 0;
            HT_pred = 0.;
            MHT_pred = 0.;
            Jet1Pt_pred = 0.;
            Jet2Pt_pred = 0.;
            Jet3Pt_pred = 0.;
            Jet4Pt_pred = 0.;
            Jet1Eta_pred = 0.;
            Jet2Eta_pred = 0.;
            Jet3Eta_pred = 0.;
            Jet4Eta_pred = 0.;
            DeltaPhi1_pred = 0.;
            DeltaPhi2_pred = 0.;
            DeltaPhi3_pred = 0.;
            DeltaPhiMinN_pred = 0.;
         }
      }
   }

   return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
int QCDBkgRS::calcNJets(const std::vector<pat::Jet>& Jets_smeared) {
   int NJets = 0;
   for (vector<pat::Jet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         ++NJets;
      }
   }
   return NJets;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
int QCDBkgRS::calcNBJets(const std::vector<pat::Jet>& Jets_smeared, std::map <const reco::Jet*, bool>& btag_map) {
   int NBJets = 0;
   for (vector<pat::Jet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      //if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_ && it->bDiscriminator(btagTag_) > btagCut_) {
      //
      if ((it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) && btag_map[&(*it)]) {
         ++NBJets;
      }
   }
   return NBJets;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
double QCDBkgRS::calcHT(const std::vector<pat::Jet>& Jets_smeared) {
   double HT = 0;
   for (vector<pat::Jet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         HT += it->pt();
      }
   }
   return HT;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
bool QCDBkgRS::calcMinDeltaPhi(const std::vector<pat::Jet>& Jets_smeared, math::PtEtaPhiMLorentzVector& MHT) {
   bool result = true;
   unsigned int i = 0;
   for (vector<pat::Jet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         if (i < JetDeltaMin_.size()) {
            if (std::abs(deltaPhi(MHT, *it)) < JetDeltaMin_.at(i))
               result = false;
            ++i;
         } else {
            break;
         }
      }
   }
   return result;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
math::PtEtaPhiMLorentzVector QCDBkgRS::calcMHT(const std::vector<pat::Jet>& Jets_smeared, double ptcut) {
   math::PtEtaPhiMLorentzVector MHT(0, 0, 0, 0);
   for (vector<pat::Jet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > ptcut && std::abs(it->eta()) < JetsMHTEta_) {
         MHT -= it->p4();
      }
   }
   return MHT;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void QCDBkgRS::FillLeadingJetPredictions(const std::vector<pat::Jet>& Jets_smeared) {
   int NJets = 0;
   for (vector<pat::Jet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         ++NJets;

         if( NJets == 1 ) {
            Jet1Pt_pred = it->pt();
            Jet1Eta_pred = it->eta();
         }
         if( NJets == 2 ) {
            Jet2Pt_pred = it->pt();
            Jet2Eta_pred = it->eta();
         }
         if( NJets == 3 ) {
            Jet3Pt_pred = it->pt();
            Jet3Eta_pred = it->eta();
         }
         if( NJets == 4 ) {
            Jet4Pt_pred = it->pt();
            Jet4Eta_pred = it->eta();
            break;
         }
      }
   }
   if( NJets < 3 ) {
      Jet3Pt_pred = -1.;
      Jet3Eta_pred = 9999.;
   }
   if( NJets < 4 ) {
      Jet4Pt_pred = -1.;
      Jet4Eta_pred = 9999.;
   }

   return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void QCDBkgRS::FillDeltaPhiPredictions(const std::vector<pat::Jet>& Jets_smeared, math::PtEtaPhiMLorentzVector& vMHT) {

   int NJets = 0;
   double dpnhat_min = 9999;

   for (vector<pat::Jet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {

      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         ++NJets;

         double deltaT = 0.0;
         float jres = 0.1;
         double sum = 0.0;
         double dphi = 0.0;
         double dpnhat = 0;

         if (NJets < 4) {
            dphi = std::abs(deltaPhi(vMHT, *it));
            for (vector<pat::Jet>::const_iterator jt = Jets_smeared.begin(); jt != Jets_smeared.end(); ++jt) {
               if (&(*it) == &(*jt)) continue;
               sum += (jt->px()*it->py()-it->px()*jt->py()) * (jt->px()*it->py()-it->px()*jt->py());
            }
            deltaT = jres*sqrt(sum)/it->pt();
         }

         if (deltaT/vMHT.pt() >= 1.0) dpnhat = dphi/(TMath::Pi()/2.0);
         else dpnhat=dphi/asin(deltaT/vMHT.pt());
         if (dpnhat < dpnhat_min) dpnhat_min = dpnhat;

         if( NJets == 1 ) {
            DeltaPhi1_pred = dphi;
         }
         if( NJets == 2 ) {
            DeltaPhi2_pred = dphi;
         }
         if( NJets == 3 ) {
            DeltaPhi3_pred = dphi;
            break;
         }
      }
   }

   DeltaPhiMinN_pred = dpnhat_min;

   if( NJets == 2 ) {
      DeltaPhi3_pred = 9999.;
   }

   return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
int QCDBkgRS::calcNJets_gen(const std::vector<reco::GenJet>& Jets_smeared) {
   int NJets = 0;
   for (vector<reco::GenJet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         ++NJets;
      }
   }
   return NJets;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
int QCDBkgRS::calcNBJets_gen(const std::vector<reco::GenJet>& Jets_smeared, std::map <const reco::GenJet*, bool>& b_map) {
   int NBJets = 0;
   for (vector<reco::GenJet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         if (b_map[&(*it)]){
            ++NBJets;
            //cout << "Found a bjet: NBJets is now " << NBJets << endl;
         }
      }
   }
   return NBJets;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
double QCDBkgRS::calcHT_gen(const std::vector<reco::GenJet>& Jets_smeared) {
   double HT = 0;
   for (vector<reco::GenJet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         HT += it->pt();
      }
   }
   return HT;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
bool QCDBkgRS::calcMinDeltaPhi_gen(const std::vector<reco::GenJet>& Jets_smeared, math::PtEtaPhiMLorentzVector& MHT) {
   bool result = true;
   unsigned int i = 0;
   for (vector<reco::GenJet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         if (i < JetDeltaMin_.size()) {
            if (std::abs(deltaPhi(MHT, *it)) < JetDeltaMin_.at(i))
               result = false;
            ++i;
         } else {
            break;
         }
      }
   }
   return result;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
math::PtEtaPhiMLorentzVector QCDBkgRS::calcMHT_gen(const std::vector<reco::GenJet>& Jets_smeared, double ptcut) {
   math::PtEtaPhiMLorentzVector MHT(0, 0, 0, 0);
   for (vector<reco::GenJet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > ptcut && std::abs(it->eta()) < JetsMHTEta_) {
         MHT -= it->p4();
      }
   }
   return MHT;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void QCDBkgRS::FillLeadingJetPredictions_gen(const std::vector<reco::GenJet>& Jets_smeared) {
   int NJets = 0;
   for (vector<reco::GenJet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {
      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         ++NJets;

         if( NJets == 1 ) {
            Jet1Pt_pred = it->pt();
            Jet1Eta_pred = it->eta();

         }
         if( NJets == 2 ) {
            Jet2Pt_pred = it->pt();
            Jet2Eta_pred = it->eta();
         }
         if( NJets == 3 ) {
            Jet3Pt_pred = it->pt();
            Jet3Eta_pred = it->eta();
         }
         if( NJets == 4 ) {
            Jet4Pt_pred = it->pt();
            Jet4Eta_pred = it->eta();
            break;
         }
      }
   }
   if( NJets < 3 ) {
      Jet3Pt_pred = -1.;
      Jet3Eta_pred = 9999.;
   }
   if( NJets < 4 ) {
      Jet4Pt_pred = -1.;
      Jet4Eta_pred = 9999.;
   }

   return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void QCDBkgRS::FillDeltaPhiPredictions_gen(const std::vector<reco::GenJet>& Jets_smeared, math::PtEtaPhiMLorentzVector& vMHT) {

   int NJets = 0;
   double dpnhat_min = 9999;

   for (vector<reco::GenJet>::const_iterator it = Jets_smeared.begin(); it != Jets_smeared.end(); ++it) {

      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         ++NJets;

         double deltaT = 0.0;
         float jres = 0.1;
         double sum = 0.0;
         double dphi = 0;
         double dpnhat = 0;

         if (NJets < 4) {
            dphi = std::abs(deltaPhi(vMHT, *it));
            for (vector<reco::GenJet>::const_iterator jt = Jets_smeared.begin(); jt != Jets_smeared.end(); ++jt) {
               if (&(*it) == &(*jt)) continue;
               sum += (jt->px()*it->py()-it->px()*jt->py()) * (jt->px()*it->py()-it->px()*jt->py());
            }
            deltaT = jres*sqrt(sum)/it->pt();
         }

         if (deltaT/vMHT.pt() >= 1.0) dpnhat = dphi/(TMath::Pi()/2.0);
         else dpnhat=dphi/asin(deltaT/vMHT.pt());
         if (dpnhat < dpnhat_min) dpnhat_min = dpnhat;

         if( NJets == 1 ) {
            DeltaPhi1_pred = dphi;
         }
         if( NJets == 2 ) {
            DeltaPhi2_pred = dphi;
         }
         if( NJets == 3 ) {
            DeltaPhi3_pred = dphi;
            break;
         }
      }
   }

   DeltaPhiMinN_pred = dpnhat_min;

   if( NJets == 2 ) {
      DeltaPhi3_pred = 9999.;
   }

   return;
}
//--------------------------------------------------------------------------


//--------------------------------------------------------------------------
void QCDBkgRS::FillPredictions(const std::vector<pat::Jet>& Jets_smeared, const int& i, const double& w, std::map <const reco::Jet*, bool>& dynamic_jet_btag_map) {

   int NJets = calcNJets(Jets_smeared);
   int NBJets = calcNBJets(Jets_smeared, dynamic_jet_btag_map);
   double HT = calcHT(Jets_smeared);
   math::PtEtaPhiMLorentzVector vMHT = calcMHT(Jets_smeared, JetsMHTPt_);
   math::PtEtaPhiMLorentzVector vMET = calcMHT(Jets_smeared, 30.);
   double MHT = vMHT.pt();
   //   bool minDeltaPhi = calcMinDeltaPhi(Jets_smeared, vMHT);

   weight = w;
   Ntries_pred = i;
   Njets_pred = NJets;
   BTags_pred = NBJets;
   HT_pred = HT;
   MHT_pred = MHT;
   FillDeltaPhiPredictions(Jets_smeared, vMET);
   FillLeadingJetPredictions(Jets_smeared);

   return;
}
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
void QCDBkgRS::FillPredictions_gen(const std::vector<reco::GenJet>& Jets_smeared, const int& i,
                                   const double& w, std::map <const reco::GenJet*, bool>& b_map ) {

   int NJets = calcNJets_gen(Jets_smeared);
   int NBJets = calcNBJets_gen(Jets_smeared, b_map);
   double HT = calcHT_gen(Jets_smeared);
   math::PtEtaPhiMLorentzVector vMHT = calcMHT_gen(Jets_smeared, JetsMHTPt_);
   math::PtEtaPhiMLorentzVector vMET = calcMHT_gen(Jets_smeared, 30.);
   double MHT = vMHT.pt();
   //  bool minDeltaPhi = calcMinDeltaPhi_gen(Jets_smeared, vMHT);

   weight = w;
   Ntries_pred = i;
   Njets_pred = NJets;
   BTags_pred = NBJets;
   HT_pred = HT;
   MHT_pred = MHT;
   FillDeltaPhiPredictions_gen(Jets_smeared, vMET);
   FillLeadingJetPredictions_gen(Jets_smeared);

   return;
}
//--------------------------------------------------------------------------

// ------------ method called for each event  ------------
void QCDBkgRS::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //LeptonVeto
   edm::Handle<int> NLeptons;
   iEvent.getByLabel(leptonTag_, NLeptons);
   if ((*NLeptons) != 0)
      return;

   //NJets
   edm::Handle<int> NJetsSeed;
   iEvent.getByLabel(NJetsSeedTag_, NJetsSeed);
   if ((*NJetsSeed) < NJetsSeedMin_)
      return;

   //HT
   edm::Handle<double> HTSeed;
   iEvent.getByLabel(HTSeedTag_, HTSeed);
   if ((*HTSeed) < HTSeedMin_)
      return;

   //Weight
   edm::Handle<double> event_weight;
   iEvent.getByLabel(weightName_, event_weight);
   weight_ = (event_weight.isValid() ? (*event_weight) : 1.0);
   //if (!event_weight.isValid()) cout << "weight not found" << endl;

   if (controlPlots_) {
      h_weight->Fill(log10(weight_));
      h_weightedWeight->Fill(log10(weight_), weight_);
   }

   // Number of vertices
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByLabel(vertices_,vertices);
   if( vertices.isValid() ) {
      vtxN = vertices->size();
   }

   //GenJets
   edm::Handle<edm::View<reco::GenJet> > gj;
   edm::View<reco::GenJet> Jets_gen;
   bool gj_present = iEvent.getByLabel(genjets_, gj);
   if (gj_present) {
      Jets_gen = *gj;
   }

   //PATJets
   edm::Handle<edm::View<pat::Jet> > Jets;
   iEvent.getByLabel(jets_, Jets);
   edm::View<pat::Jet> Jets_rec = *Jets;

   //GenParticles
   edm::Handle<edm::View<pat::PackedGenParticle> > gp;
   iEvent.getByLabel("packedGenParticles", gp);
   edm::View<pat::PackedGenParticle> genParticles = *gp;

   // collection of rebalanced jets
   std::auto_ptr<vector<pat::Jet> > Jets_reb(new vector<pat::Jet> );

   // collection of smeared jets
   std::auto_ptr<vector<pat::Jet> > Jets_smeared(new vector<pat::Jet> );

   // collection of smeared gen jets
   std::auto_ptr<vector<reco::GenJet> > GenJets_smeared(new vector<reco::GenJet> );

   // ------------------------------------------------------------------------ //
   // plots for reco jets
   HT_seed = 0;
   double HT_rec = 0;
   double HTall_rec = 0;
   double HTlow_rec = 0;
   double HThigh_rec = 0;
   const reco::Jet* Jet1_rec = 0;
   const reco::Jet* Jet2_rec = 0;
   math::PtEtaPhiMLorentzVector vMHTall_rec(0., 0., 0., 0.);
   math::PtEtaPhiMLorentzVector vMHTlow_rec(0., 0., 0., 0.);
   math::PtEtaPhiMLorentzVector vMHThigh_rec(0., 0., 0., 0.);
   int NJets_reco = 0;
   int JetCounter = 0;
   //// Fill measured particles to vector
   for (edm::View<pat::Jet>::const_iterator it = Jets_rec.begin(); it != Jets_rec.end(); ++it) {
      JetCounter++;
      h_JetPt_reco->Fill(it->pt(), weight_);

      if( JetCounter == 1 ) {
         Jet1_rec = &(*it);
      }
      if( JetCounter == 2 ) {
         Jet2_rec = &(*it);
      }

      // all HT cuts
      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         HT_seed += it->pt();
         NJets_reco++;
         HT_rec += it->pt();
      }

      // HT cuts without eta
      if (it->pt() < JetsHTPt_) {
         HTall_rec += it->pt();
         HTlow_rec += it->pt();
      } else {
         HTall_rec += it->pt();
         HThigh_rec += it->pt();
      }

      // MHT cuts without eta
      if (it->pt() < JetsMHTPt_) {
         vMHTall_rec -= it->p4();
         vMHTlow_rec -= it->p4();
      } else {
         vMHTall_rec -= it->p4();
         vMHThigh_rec -= it->p4();
      }
   }

   ////////// checks on high MHT events
   if (vMHThigh_rec.pt()>200 && debug > 0){
      int NJ = 0;
      cout << "-------------------------------" << endl;
      cout << "RecoJet summary: " << endl;
      for (edm::View<pat::Jet>::const_iterator it = Jets_rec.begin(); it != Jets_rec.end(); ++it) {
         if (NJ < 5) {
            cout << "Jet number: " << NJ << endl;
            cout << "pt, eta, phi: " << it->pt() << ", " << it->eta() << ", " << it->phi() << endl;
            cout << "neutral hadron fraction: " << it->neutralHadronEnergyFraction() << endl;
            cout << "charged hadron fraction: " << it->chargedHadronEnergyFraction() << endl;
            cout << "charged EM fraction: " << it->chargedEmEnergyFraction() << endl;
            cout << "photon fraction: " << it->neutralEmEnergyFraction() << endl;
            cout << "muon fraction: " << it->muonEnergyFraction() << endl;
            cout << "N constituents: " << it->nConstituents() << endl;
            cout << "charged hadron multiplicity: " << it->chargedHadronMultiplicity() << endl;
            double dRmin = 100;
            const reco::GenJet* matchedGenJet = 0;
            for (edm::View<reco::GenJet>::const_iterator jt = Jets_gen.begin(); jt != Jets_gen.end(); ++jt) {
               double dR = deltaR(*jt, *it);
               if (dR < dRmin) {
                  dRmin = dR;
                  matchedGenJet = &(*jt);
               }
            }
            if (dRmin < 0.5) {
               cout << "Matched genJet dR: " << dRmin << endl;
               cout << "pt, eta, phi: " << matchedGenJet->pt() << ", " << matchedGenJet->eta() << ", " << matchedGenJet->phi() << endl;
            } else {
               cout << "NOT MATCHED!" << endl;
            }
            ++NJ;
         }
      }
      cout << " " << endl;
      cout << "GenJet summary: " << endl;
      int NGJ = 0;
      for (edm::View<reco::GenJet>::const_iterator jt = Jets_gen.begin(); jt != Jets_gen.end(); ++jt) {
         if (NGJ < 5) {
            cout << "Jet number: " << NGJ << endl;
            cout << "pt, eta, phi: " << jt->pt() << ", " << jt->eta() << ", " << jt->phi() << endl;
            double dRmin = 100;
            const reco::Jet* matchedJet = 0;
            for (edm::View<pat::Jet>::const_iterator it = Jets_rec.begin(); it != Jets_rec.end(); ++it) {
               double dR = deltaR(*jt, *it);
               if (dR < dRmin) {
                  dRmin = dR;
                  matchedJet = &(*it);
               }
            }
            if (dRmin < 0.5) {
               cout << "Matched recoJet dR: " << dRmin << endl;
               cout << "pt, eta, phi: " << matchedJet->pt() << ", " << matchedJet->eta() << ", " << matchedJet->phi() << endl;
            } else {
               cout << "NOT MATCHED!" << endl;
            }
            ++NGJ;
         }
      }
   }

   // fill control plots for reasonable event weights
   if (controlPlots_ && weight_ < 30000) {

      // reco HT for different NJet bins
      if( NJets_reco == 2 ){
         h_HT_JetBin1_rec->Fill(HT_rec, weight_);
      }
      if( NJets_reco == 3 ||  NJets_reco == 4 || NJets_reco == 5){
         h_HT_JetBin2_rec->Fill(HT_rec, weight_);
      }
      if( NJets_reco == 6 ||  NJets_reco == 7){
         h_HT_JetBin3_rec->Fill(HT_rec, weight_);
      }
      if( NJets_reco >= 8){
         h_HT_JetBin4_rec->Fill(HT_rec, weight_);
      }

      // delta Phi between first two reco jets
      if( JetCounter >= 2 ){
         h_deltaPhiJet1Jet2_rec->Fill(std::abs(deltaPhi(Jet1_rec->phi(), Jet2_rec->phi())), weight_);
      }

      // control plots for reco jets
      h_nJets_reco->Fill(NJets_reco, weight_);
      h_HT_rec->Fill(HT_rec, weight_);
      h_HTall_rec->Fill(HTall_rec, weight_);
      h_HThigh_rec->Fill(HThigh_rec, weight_);
      h_MHTall_rec->Fill(vMHTall_rec.pt(), weight_);
      h_MHThigh_rec->Fill(vMHThigh_rec.pt(), weight_);
   }
   // ------------------------------------------------------------------------ //

   // ------------------------------------------------------------------------ //
   // count seed events for PU Uncertainty
   if (controlPlots_ && weight_ < 30000) {
      if( NJets_reco == 2 ) h_SeedEvents_HT_NJet2->Fill(HT_rec, vtxN, weight_);
      if( NJets_reco == 3 ) h_SeedEvents_HT_NJet3->Fill(HT_rec, vtxN, weight_);
      if( NJets_reco == 4 ) h_SeedEvents_HT_NJet4->Fill(HT_rec, vtxN, weight_);
      if( NJets_reco == 5 ) h_SeedEvents_HT_NJet5->Fill(HT_rec, vtxN, weight_);
      if( NJets_reco == 6 ) h_SeedEvents_HT_NJet6->Fill(HT_rec, vtxN, weight_);
      if( NJets_reco == 7 ) h_SeedEvents_HT_NJet7->Fill(HT_rec, vtxN, weight_);
      if( NJets_reco >= 8 ) h_SeedEvents_HT_NJet8->Fill(HT_rec, vtxN, weight_);
   }
   // ------------------------------------------------------------------------ //

   // ------------------------------------------------------------------------ //
   // plots for gen jets
   double HT_gen = 0;
   double HTall_gen = 0;
   double HTlow_gen = 0;
   double HThigh_gen = 0;
   math::PtEtaPhiMLorentzVector vMHTall_gen(0., 0., 0., 0.);
   math::PtEtaPhiMLorentzVector vMHTlow_gen(0., 0., 0., 0.);
   math::PtEtaPhiMLorentzVector vMHThigh_gen(0., 0., 0., 0.);
   int NJets_gen = 0;
   const reco::GenJet* Jet1_gen = 0;
   const reco::GenJet* Jet2_gen = 0;
   JetCounter = 0;
   //// Fill measured particles to vector
   for (edm::View<reco::GenJet>::const_iterator it = Jets_gen.begin(); it != Jets_gen.end(); ++it) {
      JetCounter++;
      h_JetPt_gen->Fill(it->pt(), weight_);

      if( JetCounter == 1 ) {
         Jet1_gen = &(*it);
      }
      if( JetCounter == 2 ) {
         Jet2_gen = &(*it);
      }

      // all HT cuts
      if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
         NJets_gen++;
         HT_gen += it->pt();
      }

      // HT cuts without eta
      if (it->pt() < JetsHTPt_) {
         HTall_gen += it->pt();
         HTlow_gen += it->pt();
      } else {
         HTall_gen += it->pt();
         HThigh_gen += it->pt();
      }

      // MHT cuts without eta
      if (it->pt() < JetsMHTPt_) {
         vMHTall_gen -= it->p4();
         vMHTlow_gen -= it->p4();
      } else {
         vMHTall_gen -= it->p4();
         vMHThigh_gen -= it->p4();
      }
   }

   // fill control plots for reasonable event weights
   if (controlPlots_ && weight_ < 30000) {

      // gen HT for different NJet bins
      if( NJets_gen == 2 ){
         h_HT_JetBin1_gen->Fill(HT_gen, weight_);
      }
      if( NJets_gen == 3 ||  NJets_gen == 4 || NJets_gen == 5){
         h_HT_JetBin2_gen->Fill(HT_gen, weight_);
      }
      if( NJets_gen == 6 ||  NJets_gen == 7){
         h_HT_JetBin3_gen->Fill(HT_gen, weight_);
      }
      if( NJets_gen >= 8){
         h_HT_JetBin4_gen->Fill(HT_gen, weight_);
      }

      // delta Phi between first two gen jets
      if( JetCounter >= 2 ){
         h_deltaPhiJet1Jet2_gen->Fill(std::abs(deltaPhi(Jet1_gen->phi(), Jet2_gen->phi())), weight_);
      }

      // control plots for gen jets
      h_nJets_gen->Fill(NJets_gen, weight_);
      h_HT_gen->Fill(HT_gen, weight_);
      h_HTall_gen->Fill(HTall_gen, weight_);
      h_HThigh_gen->Fill(HThigh_gen, weight_);
      h_MHTall_gen->Fill(vMHTall_gen.pt(), weight_);
      h_MHThigh_gen->Fill(vMHThigh_gen.pt(), weight_);
   }
   //cout << "HT gen  = " << HThigh_gen << endl;
   // ------------------------------------------------------------------------ //

   // ------------------------------------------------------------------------ //
   // plots for matched/not-matched reco - gen jets
   if (gj_present) {
      int NNN = 0;
      for (edm::View<reco::GenJet>::const_iterator it = Jets_gen.begin(); it != Jets_gen.end(); ++it) {
         double dRmin = 100;
         const reco::Jet* matchedJet = 0;
         double addActivity = 0;
         for (edm::View<pat::Jet>::const_iterator jt = Jets_rec.begin(); jt != Jets_rec.end(); ++jt) {
            double dR = deltaR(*jt, *it);
            if (dR < 0.7) {
               addActivity += jt->pt();
            }
            if (dR < dRmin) {
               dRmin = dR;
               matchedJet = &(*jt);
            }
         }
         addActivity -= matchedJet->pt();
         if (controlPlots_ && weight_ < 30000) {
            if (dRmin < 0.5) {
               // pt spectrum of all matched reco jets
               h_RecJetMatched_Pt->Fill(matchedJet->pt(), weight_);
               if (vMHThigh_rec.pt() > 100.){
                  if( NNN == 0 ) {
                     h_DeltaPhiRecoGenJet1_GenHT_HighMHT->Fill(dRmin, HT_gen, weight_);
                     h_AddRelJetActivity07GenJet1_GenHT_HighMHT->Fill(addActivity/matchedJet->pt(), HT_gen, weight_);
                  }
                  if( NNN == 1 ){
                     h_DeltaPhiRecoGenJet2_GenHT_HighMHT->Fill(dRmin, HT_gen, weight_);
                     h_AddRelJetActivity07GenJet2_GenHT_HighMHT->Fill(addActivity/matchedJet->pt(), HT_gen, weight_);
                  }
                  if( NNN == 2 ){
                     h_DeltaPhiRecoGenJet3_GenHT_HighMHT->Fill(dRmin, HT_gen, weight_);
                     h_AddRelJetActivity07GenJet3_GenHT_HighMHT->Fill(addActivity/matchedJet->pt(), HT_gen, weight_);
                  }
               } else {
                  if( NNN == 0 ) {
                     h_DeltaPhiRecoGenJet1_GenHT_LowMHT->Fill(dRmin, HT_gen, weight_);
                     h_AddRelJetActivity07GenJet1_GenHT_LowMHT->Fill(addActivity/matchedJet->pt(), HT_gen, weight_);
                  }
                  if( NNN == 1 ){
                     h_DeltaPhiRecoGenJet2_GenHT_LowMHT->Fill(dRmin, HT_gen, weight_);
                     h_AddRelJetActivity07GenJet2_GenHT_LowMHT->Fill(addActivity/matchedJet->pt(), HT_gen, weight_);
                  }
                  if( NNN == 2 ){
                     h_DeltaPhiRecoGenJet3_GenHT_LowMHT->Fill(dRmin, HT_gen, weight_);
                     h_AddRelJetActivity07GenJet3_GenHT_LowMHT->Fill(addActivity/matchedJet->pt(), HT_gen, weight_);
                  }
               }
               // pt spectrum of matched reco jets in NJet bins
               if( NJets_gen == 2 ) {
                  h_RecJetMatched_JetBin1_Pt->Fill(matchedJet->pt(), weight_);
               }
               if( NJets_gen == 3 ||  NJets_gen == 4 || NJets_gen == 5){
                  h_RecJetMatched_JetBin2_Pt->Fill(matchedJet->pt(), weight_);
               }
               if( NJets_gen == 6 ||  NJets_gen == 7){
                  h_RecJetMatched_JetBin3_Pt->Fill(matchedJet->pt(), weight_);
               }
               if( NJets_gen >= 8){
                  h_RecJetMatched_JetBin4_Pt->Fill(matchedJet->pt(), weight_);
               }

               // resolution of reco jets compared to gen jets
               if (fabs(it->eta()) < 1.5)
                  h_RecJetRes_Pt->Fill(it->pt(), matchedJet->pt() / it->pt(), weight_);
               if (it->pt() > 100.)
                  h_RecJetRes_Eta->Fill(it->eta(), matchedJet->pt() / it->pt(), weight_);
            }
            else if (matchedJet != 0) {
               // pt spectrum of not-matched reco jets
               h_RecJetNotMatched_Pt->Fill(matchedJet->pt(), weight_);

               // pt spectrum of not-matched reco jets in NJet bins
               if( NJets_gen == 2 ) {
                  h_RecJetNotMatched_JetBin1_Pt->Fill(matchedJet->pt(), weight_);
               }
               if( NJets_gen == 3 ||  NJets_gen == 4 || NJets_gen == 5){
                  h_RecJetNotMatched_JetBin2_Pt->Fill(matchedJet->pt(), weight_);
               }
               if( NJets_gen == 6 ||  NJets_gen == 7){
                  h_RecJetNotMatched_JetBin3_Pt->Fill(matchedJet->pt(), weight_);
               }
               if( NJets_gen >= 8){
                  h_RecJetNotMatched_JetBin4_Pt->Fill(matchedJet->pt(), weight_);
               }
            }
         }
         ++NNN;
      }
   }
   // ------------------------------------------------------------------------ //

   Jets_reb->reserve(Jets_rec.size());
   Jets_smeared->reserve(Jets_rec.size());

   //
   // Rebalance multi jet system
   //
   bool isRebalanced = false;

   if (smearCollection_ == "Reco") {
      if (Jets_rec.size() > 2) {
         isRebalanced = RebalanceJets_KinFitter(&Jets_rec, *(Jets_reb.get()));
      }

      if (!isRebalanced) {
         //cout << "Bad event: Not possible to rebalance!" << endl;
         weight_ = 0;
      }

      // sort rebalanced jets
      GreaterByPt<pat::Jet> ptComparator_;
      std::sort(Jets_reb->begin(), Jets_reb->end(), ptComparator_);

      // plots for reb jets
      double HT_reb = 0;
      double HThigh_reb = 0;
      math::PtEtaPhiMLorentzVector vMHThigh_reb(0., 0., 0., 0.);
      double HTlow_reb = 0;
      math::PtEtaPhiMLorentzVector vMHTlow_reb(0., 0., 0., 0.);
      double HTall_reb = 0;
      math::PtEtaPhiMLorentzVector vMHTall_reb(0., 0., 0., 0.);
      int NJets_reb = 0;
      const pat::Jet* Jet1_reb = 0;
      const pat::Jet* Jet2_reb = 0;
      JetCounter = 0;
      for (vector<pat::Jet>::const_iterator it = Jets_reb-> begin(); it != Jets_reb->end(); ++it) {
         JetCounter++;
         h_JetPt_reb->Fill(it->pt(), weight_);

         if( JetCounter == 1 ) {
            Jet1_reb = &(*it);
         }
         if( JetCounter == 2 ) {
            Jet2_reb = &(*it);
         }

         // all HT cuts
         if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
            NJets_reb++;
            HT_reb += it->pt();
         }

         // HT cuts without eta
         if (it->pt() < JetsHTPt_) {
            HTall_reb += it->pt();
            HTlow_reb += it->pt();
         } else {
            HTall_reb += it->pt();
            HThigh_reb += it->pt();
         }

         // MHT cuts without eta
         if (it->pt() < JetsMHTPt_) {
            vMHTall_reb -= it->p4();
            vMHTlow_reb -= it->p4();
         } else {
            vMHTall_reb -= it->p4();
            vMHThigh_reb -= it->p4();
         }
      }

      if (!isRebalanced) {
         cout << "Bad event: Can't be rebalanced!!!" << endl;
         cout << "Reconstructed: HT, MHT = " << HThigh_rec << ", " << vMHThigh_rec.pt() << endl;
         cout << "Rebalanced: HT, MHT = " << HTall_reb << ", " << vMHTall_reb.pt() << endl;
      }

      // fill control plots for reasonable event weights
      if (controlPlots_ && weight_ < 30000) {

         // reb HT for different NJet bins
         if( NJets_reb == 2 ){
            h_HT_JetBin1_reb->Fill(HT_reb, weight_);
         }
         if( NJets_reb == 3 ||  NJets_reb == 4 || NJets_reb == 5){
            h_HT_JetBin2_reb->Fill(HT_reb, weight_);
         }
         if( NJets_reb == 6 ||  NJets_reb == 7){
            h_HT_JetBin3_reb->Fill(HT_reb, weight_);
         }
         if( NJets_reb >= 8){
            h_HT_JetBin4_reb->Fill(HT_reb, weight_);
         }

         // delta Phi between first two reb jets
         if( JetCounter >= 2 ){
            h_deltaPhiJet1Jet2_reb->Fill(std::abs(deltaPhi(Jet1_reb->phi(), Jet2_reb->phi())), weight_);
         }

         // control plots for reb jets
         h_nJets_reb->Fill(NJets_reb, weight_);
         h_HT_reb->Fill(HT_reb, weight_);
         h_HTall_reb->Fill(HTall_reb, weight_);
         h_HThigh_reb->Fill(HThigh_reb, weight_);
         //         if (abs(HThigh_reb - HThigh_rec) > 100) {
         //            cout << "WARNING!!!!!!!!!!!!!!!!!!!" << endl;
         //            cout << uncertaintyName_ << ": HTall reb = " << HTall_reb << " HThigh reb = " << HThigh_reb
         //                  << " HTall rec = " << HTall_rec << " HThigh rec = " << HThigh_rec << " weight = " << weight_ << endl;
         //         }
         h_MHTall_reb->Fill(vMHTall_reb.pt(), weight_);
         h_MHThigh_reb->Fill(vMHThigh_reb.pt(), weight_);
      }
      // ------------------------------------------------------------------------ //

      // ------------------------------------------------------------------------ //
      //// plots for matched/not-matched reb - gen jets
      if (smearCollection_ == "Reco") {
         if (gj_present) {
            for (edm::View<reco::GenJet>::const_iterator it = Jets_gen.begin(); it != Jets_gen.end(); ++it) {
               double dRmin = 100;
               const pat::Jet* matchedJet = 0;
               // match reb jets to gen jets
               for (vector<pat::Jet>::const_iterator jt = Jets_reb-> begin(); jt != Jets_reb->end(); ++jt) {
                  double dR = deltaR(*jt, *it);
                  if (dR < dRmin) {
                     dRmin = dR;
                     matchedJet = &(*jt);
                  }
               }
               // match reco jets to gen jets
               double dRmin_reco = 100;
               const reco::Jet* matchedJet_reco = 0;
               for (edm::View<pat::Jet>::const_iterator jt = Jets_rec.begin(); jt != Jets_rec.end(); ++jt) {
                  double dR_reco = deltaR(*jt, *it);
                  if (dR_reco < dRmin_reco) {
                     dRmin_reco = dR_reco;
                     matchedJet_reco = &(*jt);
                  }
               }

               if (controlPlots_ && weight_ < 30000) {
                  if (dRmin < 0.15) {
                     //// Find additional neutrinos and add them back since these are not included in the genJets
                     math::PtEtaPhiMLorentzVector neutrinos(0., 0., 0., 0.);
                     for(edm::View<pat::PackedGenParticle>::const_iterator cand = genParticles.begin(); cand!=genParticles.end(); ++cand)
                     {
                        if ( cand->status()==1 && (abs(cand->pdgId())==12 || abs(cand->pdgId())==14 || abs(cand->pdgId())==16)){
                           double dR = deltaR(*it, *cand);
                           if (dR < 0.4) neutrinos += cand->p4();
                        }
                     }

                     // resolution of reb jets compared to gen jets
                     if (fabs(it->eta()) < 1.5)
                        h_RebJetRes_Pt->Fill((it->p4()+neutrinos).pt(), matchedJet->pt() / (it->p4()+neutrinos).pt(), weight_);
                     if (it->pt() > 100.)
                        h_RebJetRes_Eta->Fill(it->eta(), matchedJet->pt() / (it->p4()+neutrinos).pt(), weight_);

                     // get correction factor for rebalancing vs. reco jet pt
                     if(dRmin_reco < 0.15) {
                        if (matchedJet->bDiscriminator(btagTag_) > btagCut_){
                           h_RebCorrection_vsReco_b->Fill(matchedJet_reco->pt(), matchedJet->pt() / (it->p4()+neutrinos).pt(), weight_);
                        }
                        else {
                           h_RebCorrection_vsReco->Fill(matchedJet_reco->pt(), matchedJet->pt() / (it->p4()+neutrinos).pt(), weight_);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   // ------------------------------------------------------------------------ //


   //
   // Smear rebalanced multi jet system
   //
   if (smearCollection_ == "Reco") {
      SmearingJets(*(Jets_reb.get()), *(Jets_smeared.get()));
   } else if (smearCollection_ == "Gen") {
      // Fill the map if genJet can be associated to reco b-tag
      if (gj_present) {
         for (edm::View<reco::GenJet>::const_iterator it = Jets_gen.begin(); it != Jets_gen.end(); ++it) {
            double dRmin = 100;
            const pat::Jet* matchedJet = 0;
            // match reb jets to gen jets
            for (edm::View<pat::Jet>::const_iterator jt = Jets_rec.begin(); jt != Jets_rec.end(); ++jt) {
               double dR = deltaR(*jt, *it);
               if (dR < dRmin) {
                  dRmin = dR;
                  matchedJet = &(*jt);
               }
            }
            if (dRmin < 0.15){
               if (matchedJet->bDiscriminator(btagTag_) > btagCut_){
                  genJet_btag[&(*it)] = true;
               } else {
                  genJet_btag[&(*it)] = false;
               }
            }
         }
      }
      SmearingGenJets(&Jets_gen, &genParticles, *(GenJets_smeared.get()));
   }

   // plots for smeared jets
   double HT_smeared = 0;
   double HThigh_smeared = 0;
   math::PtEtaPhiMLorentzVector vMHThigh_smeared(0., 0., 0., 0.);
   int NJets_smeared = 0;
   const pat::Jet* Jet1_smear = 0;
   const pat::Jet* Jet2_smear = 0;
   JetCounter = 0;
   if (smearCollection_ == "Reco") {
      for (vector<pat::Jet>::const_iterator it = Jets_smeared-> begin(); it != Jets_smeared->end(); ++it) {
         JetCounter++;
         h_JetPt_smear->Fill(it->pt(), weight_);

         if( JetCounter == 1 ) {
            Jet1_smear = &(*it);
         }
         if( JetCounter == 2 ) {
            Jet2_smear = &(*it);
         }

         // all HT cuts
         if (it->pt() > JetsHTPt_ && std::abs(it->eta()) < JetsHTEta_) {
            NJets_smeared++;
            HT_smeared += it->pt();
         }

         vMHThigh_smeared -= it->p4();
         HThigh_smeared += it->pt();
      }
   } else if (smearCollection_ == "Gen") {
      for (vector<reco::GenJet>::const_iterator it = GenJets_smeared-> begin(); it != GenJets_smeared->end(); ++it) {
         vMHThigh_smeared -= it->p4();
         HThigh_smeared += it->pt();
      }
   }

   // fill control plots for reasonable event weights
   if (controlPlots_ && weight_ < 30000) {

      // smeared HT for different NJet bins
      if( NJets_smeared == 2 ){
         h_HT_JetBin1_smeared->Fill(HT_smeared, weight_);
      }
      if( NJets_smeared == 3 ||  NJets_smeared == 4 || NJets_smeared == 5){
         h_HT_JetBin2_smeared->Fill(HT_smeared, weight_);
      }
      if( NJets_smeared == 6 ||  NJets_smeared == 7){
         h_HT_JetBin3_smeared->Fill(HT_smeared, weight_);
      }
      if( NJets_smeared >= 8){
         h_HT_JetBin4_smeared->Fill(HT_smeared, weight_);
      }

      // delta Phi between first two smeared jets
      if( JetCounter >= 2 ){
         h_deltaPhiJet1Jet2_smeared->Fill(std::abs(deltaPhi(Jet1_smear->phi(), Jet2_smear->phi())), weight_);
      }

      // control plots for smeared jets
      h_nJets_smear->Fill(NJets_smeared, weight_);
      h_HT_smeared->Fill(HT_smeared, weight_);
      h_HTall_smeared->Fill(HTlow_rec + HThigh_smeared, weight_);
      h_HThigh_smeared->Fill(HThigh_smeared, weight_);
      h_MHTall_smeared->Fill((vMHTlow_rec + vMHThigh_smeared).pt(), weight_);
      h_MHThigh_smeared->Fill(vMHThigh_smeared.pt(), weight_);
   }
   // ------------------------------------------------------------------------ //

   // ------------------------------------------------------------------------ //
   //// plots for matched/not-matched smeared gen - gen jets
   if (smearCollection_ == "Gen") {
      if (gj_present) {
         for (edm::View<reco::GenJet>::const_iterator it = Jets_gen.begin(); it != Jets_gen.end(); ++it) {
            double dRmin = 100;
            const reco::GenJet* matchedJet = 0;
            for (vector<reco::GenJet>::const_iterator jt = GenJets_smeared-> begin(); jt != GenJets_smeared->end(); ++jt) {
               double dR = deltaR(*jt, *it);
               if (dR < dRmin) {
                  dRmin = dR;
                  matchedJet = &(*jt);
               }
            }
            if (controlPlots_ && weight_ < 30000) {
               if (dRmin < 0.15) {
                  // resolution of smeared gen jets compared to gen jets
                  if (fabs(it->eta()) < 1.5)
                     h_SmearedJetRes_Pt->Fill(it->pt(), matchedJet->pt() / it->pt(), weight_);
                  if (it->pt() > 100.)
                     h_SmearedJetRes_Eta->Fill(it->eta(), matchedJet->pt() / it->pt(), weight_);
               }
            }
         }
      }
   }
   //// plots for matched/not-matched smeared - gen jets
   else if (smearCollection_ == "Reco") {
      if (gj_present) {
         for (edm::View<reco::GenJet>::const_iterator it = Jets_gen.begin(); it != Jets_gen.end(); ++it) {
            double dRmin = 100;
            const pat::Jet* matchedJet = 0;
            for (vector<pat::Jet>::const_iterator jt = Jets_smeared-> begin(); jt != Jets_smeared->end(); ++jt) {
               double dR = deltaR(*jt, *it);
               if (dR < dRmin) {
                  dRmin = dR;
                  matchedJet = &(*jt);
               }
            }
            if (controlPlots_ && weight_ < 30000) {
               // resolution of smeared jets compared to gen jets
               if (dRmin < 0.15) {
                  if (fabs(it->eta()) < 1.5)
                     h_SmearedJetRes_Pt->Fill(it->pt(), matchedJet->pt() / it->pt(), weight_);
                  if (it->pt() > 100.)
                     h_SmearedJetRes_Eta->Fill(it->eta(), matchedJet->pt() / it->pt(), weight_);
               }
            }
         }
      }
   }
}
//--------------------------------------------------------------------------


// ------------ method called once each job just before starting event loop  ------------
void QCDBkgRS::beginJob()
{

   debug = 0;

   edm::Service<TFileService> fs;
   if (!fs) {
      throw edm::Exception(edm::errors::Configuration, "TFile Service is not registered in cfg file");
   }

   if (controlPlots_) {
      h_SeedEvents_HT_NJet2 = fs->make<TH2F> ("SeedEvents_HT_NJet2", "Seed Events", NbinsHT_, HTmin_, HTmax_, 60, 0, 60 );
      h_SeedEvents_HT_NJet2->Sumw2();
      h_SeedEvents_HT_NJet3 = fs->make<TH2F> ("SeedEvents_HT_NJet3", "Seed Events", NbinsHT_, HTmin_, HTmax_, 60, 0, 60 );
      h_SeedEvents_HT_NJet3->Sumw2();
      h_SeedEvents_HT_NJet4 = fs->make<TH2F> ("SeedEvents_HT_NJet4", "Seed Events", NbinsHT_, HTmin_, HTmax_, 60, 0, 60 );
      h_SeedEvents_HT_NJet4->Sumw2();
      h_SeedEvents_HT_NJet5 = fs->make<TH2F> ("SeedEvents_HT_NJet5", "Seed Events", NbinsHT_, HTmin_, HTmax_, 60, 0, 60 );
      h_SeedEvents_HT_NJet5->Sumw2();
      h_SeedEvents_HT_NJet6 = fs->make<TH2F> ("SeedEvents_HT_NJet6", "Seed Events", NbinsHT_, HTmin_, HTmax_, 60, 0, 60 );
      h_SeedEvents_HT_NJet6->Sumw2();
      h_SeedEvents_HT_NJet7 = fs->make<TH2F> ("SeedEvents_HT_NJet7", "Seed Events", NbinsHT_, HTmin_, HTmax_, 60, 0, 60 );
      h_SeedEvents_HT_NJet7->Sumw2();
      h_SeedEvents_HT_NJet8 = fs->make<TH2F> ("SeedEvents_HT_NJet8", "Seed Events", NbinsHT_, HTmin_, HTmax_, 60, 0, 60 );
      h_SeedEvents_HT_NJet8->Sumw2();

      h_DeltaPhiRecoGenJet1_GenHT_LowMHT = fs->make<TH2F> ("h_DeltaPhiRecoGenJet1_GenHT_LowMHT", "h_DeltaPhiRecoGenJet1_GenHT_LowMHT", 100, 0, 0.5, 80, 0, 4000 );
      h_DeltaPhiRecoGenJet1_GenHT_LowMHT->Sumw2();
      h_DeltaPhiRecoGenJet2_GenHT_LowMHT = fs->make<TH2F> ("h_DeltaPhiRecoGenJet2_GenHT_LowMHT", "h_DeltaPhiRecoGenJet2_GenHT_LowMHT", 100, 0, 0.5, 80, 0, 4000 );
      h_DeltaPhiRecoGenJet2_GenHT_LowMHT->Sumw2();
      h_DeltaPhiRecoGenJet3_GenHT_LowMHT = fs->make<TH2F> ("h_DeltaPhiRecoGenJet3_GenHT_LowMHT", "h_DeltaPhiRecoGenJet3_GenHT_LowMHT", 100, 0, 0.5, 80, 0, 4000 );
      h_DeltaPhiRecoGenJet3_GenHT_LowMHT->Sumw2();

      h_DeltaPhiRecoGenJet1_GenHT_HighMHT = fs->make<TH2F> ("h_DeltaPhiRecoGenJet1_GenHT_HighMHT", "h_DeltaPhiRecoGenJet1_GenHT_HighMHT", 100, 0, 0.5, 80, 0, 4000 );
      h_DeltaPhiRecoGenJet1_GenHT_HighMHT->Sumw2();
      h_DeltaPhiRecoGenJet2_GenHT_HighMHT = fs->make<TH2F> ("h_DeltaPhiRecoGenJet2_GenHT_HighMHT", "h_DeltaPhiRecoGenJet2_GenHT_HighMHT", 100, 0, 0.5, 80, 0, 4000 );
      h_DeltaPhiRecoGenJet2_GenHT_HighMHT->Sumw2();
      h_DeltaPhiRecoGenJet3_GenHT_HighMHT = fs->make<TH2F> ("h_DeltaPhiRecoGenJet3_GenHT_HighMHT", "h_DeltaPhiRecoGenJet3_GenHT_HighMHT", 100, 0, 0.5, 80, 0, 4000 );
      h_DeltaPhiRecoGenJet3_GenHT_HighMHT->Sumw2();

      h_AddRelJetActivity07GenJet1_GenHT_LowMHT = fs->make<TH2F> ("h_AddRelJetActivity07GenJet1_GenHT_LowMHT", "h_AddRelJetActivity07GenJet1_GenHT_LowMHT", 100, 0, 1., 80, 0, 4000 );
      h_AddRelJetActivity07GenJet1_GenHT_LowMHT->Sumw2();
      h_AddRelJetActivity07GenJet2_GenHT_LowMHT = fs->make<TH2F> ("h_AddRelJetActivity07GenJet2_GenHT_LowMHT", "h_AddRelJetActivity07GenJet2_GenHT_LowMHT", 100, 0, 1., 80, 0, 4000 );
      h_AddRelJetActivity07GenJet2_GenHT_LowMHT->Sumw2();
      h_AddRelJetActivity07GenJet3_GenHT_LowMHT = fs->make<TH2F> ("h_AddRelJetActivity07GenJet3_GenHT_LowMHT", "h_AddRelJetActivity07GenJet3_GenHT_LowMHT", 100, 0, 1., 80, 0, 4000 );
      h_AddRelJetActivity07GenJet3_GenHT_LowMHT->Sumw2();

      h_AddRelJetActivity07GenJet1_GenHT_HighMHT = fs->make<TH2F> ("h_AddRelJetActivity07GenJet1_GenHT_HighMHT", "h_AddRelJetActivity07GenJet1_GenHT_HighMHT", 100, 0, 1., 80, 0, 4000 );
      h_AddRelJetActivity07GenJet1_GenHT_HighMHT->Sumw2();
      h_AddRelJetActivity07GenJet2_GenHT_HighMHT = fs->make<TH2F> ("h_AddRelJetActivity07GenJet2_GenHT_HighMHT", "h_AddRelJetActivity07GenJet2_GenHT_HighMHT", 100, 0, 1., 80, 0, 4000 );
      h_AddRelJetActivity07GenJet2_GenHT_HighMHT->Sumw2();
      h_AddRelJetActivity07GenJet3_GenHT_HighMHT = fs->make<TH2F> ("h_AddRelJetActivity07GenJet3_GenHT_HighMHT", "h_AddRelJetActivity07GenJet3_GenHT_HighMHT", 100, 0, 1., 80, 0, 4000 );
      h_AddRelJetActivity07GenJet3_GenHT_HighMHT->Sumw2();

      h_nJets_gen = fs->make<TH1F> ("NJets_gen", "NJets", 15, 0., 15);
      h_nJets_gen->Sumw2();
      h_nJets_reco = fs->make<TH1F> ("NJets_reco", "NJets", 15, 0., 15);
      h_nJets_reco->Sumw2();
      h_nJets_reb = fs->make<TH1F> ("NJets_reb", "NJets", 15, 0., 15);
      h_nJets_reb->Sumw2();
      h_nJets_smear = fs->make<TH1F> ("NJets_smear", "NJets", 15, 0., 15);
      h_nJets_smear->Sumw2();

      // Control plots for btag corrections
      h_pBTrue_smear = fs->make<TH1F> ("pBTrue_smear", "pBTrue", 1000, 0., 100.);
      h_pBTrue_smear->Sumw2();
      h_pBTag_smear = fs->make<TH1F> ("pBTag_smear", "pBTag", 1000, 0., 100.);
      h_pBTag_smear->Sumw2();
      h_BTagCorrectionFactor_smear = fs->make<TH1F> ("BTagCorrectionFactor_smear", "BTagCorr", 100000, 0., 1000.);
      h_BTagCorrectionFactor_smear->Sumw2();
      h_NB_Old = fs->make<TH1F> ("NB_Old", "NBOld", 15, 0., 15.);
      h_NB_Old->Sumw2();
      h_NB_New = fs->make<TH1F> ("NB_New", "NBNew", 15, 0., 15.);
      h_NB_New->Sumw2();
      h_random_BTag_smear = fs->make<TH1F> ("Random_Number_BTag_smear", "BTagCorr", 1000, 0., 1.);
      h_random_BTag_smear->Sumw2();
      h_random_BTrue_smear = fs->make<TH1F> ("Random_Number_BTrue_smear", "BTagCorr", 1000, 0., 1.);
      h_random_BTrue_smear->Sumw2();
      //EtaBinEdges.push_back(0.0); // Like in MCResolutions
      //EtaBinEdges.push_back(0.3);
      //EtaBinEdges.push_back(0.5);
      //EtaBinEdges.push_back(0.8);
      //EtaBinEdges.push_back(1.1);
      //EtaBinEdges.push_back(1.4);
      //EtaBinEdges.push_back(1.7);
      //EtaBinEdges.push_back(2.0);
      //EtaBinEdges.push_back(2.3);
      //EtaBinEdges.push_back(2.8);
      //EtaBinEdges.push_back(3.2);
      //EtaBinEdges.push_back(4.1);
      //EtaBinEdges.push_back(5.0);
      //for (unsigned int i_eta = 0; i_eta < EtaBinEdges.size() - 1; ++i_eta) {
      //      char hname[100];
      //      cout << "Book histograms for btag effs" << endl;
      //      sprintf(hname, "h_trueb_SmearedPt_Eta%i", i_eta);
      //      h_trueb_SmearedPt_temp = fs->make<TH1F> (hname, hname, 200, 0., 1000.);
      //      h_trueb_SmearedPt.at(i_eta) = *h_trueb_SmearedPt_temp;
      //      h_trueb_SmearedPt.at(i_eta)->Sumw2();
      //      sprintf(hname, "h_no_trueb_SmearedPt_Eta%i", i_eta);
      //      h_no_trueb_SmearedPt.at(i_eta) = fs->make<TH1F> (hname, hname, 200, 0., 1000.);
      //      h_no_trueb_SmearedPt.at(i_eta)->Sumw2();
      //      sprintf(hname, "h_trueb_btag_SmearedPt_Eta%i", i_eta);
      //      h_trueb_btag_SmearedPt.at(i_eta) = fs->make<TH1F> (hname, hname, 200, 0., 1000.);
      //      h_trueb_btag_SmearedPt.at(i_eta)->Sumw2();
      //      sprintf(hname, "h_no_trueb_btag_SmearedPt_Eta%i", i_eta);
      //      h_no_trueb_btag_SmearedPt.at(i_eta) = fs->make<TH1F> (hname, hname, 200, 0., 1000.);
      //      h_no_trueb_btag_SmearedPt.at(i_eta)->Sumw2();
      //      sprintf(hname, "h_btag_SmearedPt_Eta%i", i_eta);
      //      h_btag_SmearedPt.at(i_eta) = fs->make<TH1F> (hname, hname, 200, 0., 1000.);
      //      h_btag_SmearedPt.at(i_eta)->Sumw2();
      //      sprintf(hname, "h_btag_old_SmearedPt_Eta%i", i_eta);
      //      h_btag_old_SmearedPt.at(i_eta) = fs->make<TH1F> (hname, hname, 200, 0., 1000.);
      //      h_btag_old_SmearedPt.at(i_eta)->Sumw2();
      //}

      cout << "Book histograms for btag effs" << endl;
      int etabins = 12; //EtaBinEdges_.size() - 1; This was 12 in MC Resolutions
      int etamax = EtaBinEdges_.back();
      int etamin = EtaBinEdges_.front();
      cout << "# eta bins: " << etabins << " and etamin: " << etamin << " and etamax: " << etamax << endl;
      h_btageff_genPt_Response = fs->make<TH2F> ("h_btageff_genPt_Response", "h_btageff_genPt_Response", 200, 0., 1000., 10, 0., 5.);
      h_btageff_genPt_Response->Sumw2();
      h_bmistageff_genPt_Response = fs->make<TH2F> ("h_bmistageff_genPt_Response", "h_bmistageff_genPt_Response", 200, 0., 1000., 10, 0, 5.);
      h_bmistageff_genPt_Response->Sumw2();
      h_btageff_SmearPt_Response = fs->make<TH2F> ("h_btageff_SmearPt_Response", "h_btageff_SmearPt_Response", 200, 0., 1000., 10, 0., 5.);
      h_btageff_SmearPt_Response->Sumw2();
      h_bmistageff_SmearPt_Response = fs->make<TH2F> ("h_bmistageff_SmearPt_Response", "h_bmistageff_SmearPt_Response", 200, 0., 1000., 10, 0., 5.);
      h_bmistageff_SmearPt_Response->Sumw2();
      h_trueb_SmearedPt_Eta = fs->make<TH2F> ("h_trueb_SmearedPt_Eta", "h_trueb_SmearedPt_Eta", 200, 0., 1000., etabins, etamin, etamax);
      h_trueb_SmearedPt_Eta->Sumw2();
      h_trueb_random_SmearedPt_Eta = fs->make<TH2F> ("h_trueb_random_SmearedPt_Eta", "h_trueb_random_SmearedPt_Eta", 200, 0., 1000., etabins, etamin, etamax);
      h_trueb_random_SmearedPt_Eta->Sumw2();
      h_trueb_genPt_Eta = fs->make<TH2F> ("h_trueb_genPt_Eta", "h_trueb_genPt_Eta", 200, 0., 1000., etabins, etamin, etamax);
      h_trueb_genPt_Eta->Sumw2();
      h_no_trueb_genPt_Eta = fs->make<TH2F> ("h_no_trueb_genPt_Eta", "h_no_trueb_genPt_Eta", 200, 0., 1000., etabins, etamin, etamax);
      h_no_trueb_genPt_Eta->Sumw2();
      h_no_trueb_SmearedPt_Eta = fs->make<TH2F> ("h_no_trueb_SmearedPt_Eta", "h_no_trueb_SmearedPt_Eta", 200, 0., 1000., etabins, etamin, etamax);
      h_no_trueb_SmearedPt_Eta->Sumw2();
      h_trueb_btag_genPt_Eta = fs->make<TH2F> ("h_trueb_btag_genPt_Eta", "h_trueb_btag_genPt_Eta", 200, 0., 1000., etabins, etamin, etamax);
      h_trueb_btag_genPt_Eta->Sumw2();
      h_trueb_btag_SmearedPt_Eta = fs->make<TH2F> ("h_trueb_btag_SmearedPt_Eta", "h_trueb_btag_SmearedPt_Eta", 200, 0., 1000., etabins, etamin, etamax);
      h_trueb_btag_SmearedPt_Eta->Sumw2();
      h_no_trueb_btag_genPt_Eta = fs->make<TH2F> ("h_no_trueb_btag_genPt_Eta", "h_no_trueb_btag_genPt_Eta", 200, 0., 1000., etabins, etamin, etamax);
      h_no_trueb_btag_genPt_Eta->Sumw2();
      h_no_trueb_btag_SmearedPt_Eta = fs->make<TH2F> ("h_no_trueb_btag_SmearedPt_Eta", "h_no_trueb_btag_SmearedPt_Eta", 200, 0., 1000., etabins, etamin, etamax);
      h_no_trueb_btag_SmearedPt_Eta->Sumw2();
      h_btag_SmearedPt_Eta = fs->make<TH2F> ("h_btag_SmearedPt_Eta", "h_btag_SmearedPt_Eta", 200, 0., 1000., etabins, etamin, etamax);
      h_btag_SmearedPt_Eta->Sumw2();
      h_btag_old_SmearedPt_Eta = fs->make<TH2F> ("h_btag_old_SmearedPt_Eta", "h_btag_old_SmearedPt_Eta", 200, 0., 1000., etabins, etamin, etamax);
      h_btag_old_SmearedPt_Eta->Sumw2();

      h_JetPt_gen = fs->make<TH1F> ("JetPt_gen", "Jet pt", 1000, 0., 1000.);
      h_JetPt_gen->Sumw2();
      h_JetPt_reco = fs->make<TH1F> ("JetPt_reco", "Jet pt", 1000, 0., 1000.);
      h_JetPt_reco->Sumw2();
      h_JetPt_reb = fs->make<TH1F> ("JetPt_reb", "Jet pt", 1000, 0., 1000.);
      h_JetPt_reb->Sumw2();
      h_JetPt_smear = fs->make<TH1F> ("JetPt_smear", "Jet pt", 1000, 0., 1000.);
      h_JetPt_smear->Sumw2();

      h_deltaR_rebCorr = fs->make<TH1F> ("deltaR_rebCorr", "deltaR", 400, 0., 2.);
      h_deltaR_rebCorr->Sumw2();

      h_RebCorrection_vsReco = fs->make<TH2F> ("RebCorrection_vsReco", "Jet pt", 1000, 0., 1000., 100, 0., 3.);
      h_RebCorrection_vsReco->Sumw2();
      h_RebCorrection_vsReco_b = fs->make<TH2F> ("RebCorrection_vsReco_b", "Jet pt", 1000, 0., 1000., 100, 0., 3.);
      h_RebCorrection_vsReco_b->Sumw2();

      h_RecJetMatched_Pt = fs->make<TH1F> ("RecJetMatched_Pt", "RecJetMatched_Pt", 1000, 0., 1000.);
      h_RecJetMatched_Pt->Sumw2();
      h_RecJetNotMatched_Pt = fs->make<TH1F> ("RecJetNotMatched_Pt", "RecJetNotMatched_Pt", 1000, 0., 1000.);
      h_RecJetNotMatched_Pt->Sumw2();
      h_RecJetMatched_JetBin1_Pt = fs->make<TH1F> ("RecJetMatched_JetBin1_Pt", "RecJetMatched_JetBin1_Pt", 1000, 0., 1000.);
      h_RecJetMatched_JetBin1_Pt->Sumw2();
      h_RecJetNotMatched_JetBin1_Pt = fs->make<TH1F> ("RecJetNotMatched_JetBin1_Pt", "RecJetNotMatched_JetBin1_Pt", 1000, 0., 1000.);
      h_RecJetNotMatched_JetBin1_Pt->Sumw2();
      h_RecJetMatched_JetBin2_Pt = fs->make<TH1F> ("RecJetMatched_JetBin2_Pt", "RecJetMatched_JetBin2_Pt", 1000, 0., 1000.);
      h_RecJetMatched_JetBin2_Pt->Sumw2();
      h_RecJetNotMatched_JetBin2_Pt = fs->make<TH1F> ("RecJetNotMatched_JetBin2_Pt", "RecJetNotMatched_JetBin2_Pt", 1000, 0., 1000.);
      h_RecJetNotMatched_JetBin2_Pt->Sumw2();
      h_RecJetMatched_JetBin3_Pt = fs->make<TH1F> ("RecJetMatched_JetBin3_Pt", "RecJetMatched_JetBin3_Pt", 1000, 0., 1000.);
      h_RecJetMatched_JetBin3_Pt->Sumw2();
      h_RecJetNotMatched_JetBin3_Pt = fs->make<TH1F> ("RecJetNotMatched_JetBin3_Pt", "RecJetNotMatched_JetBin3_Pt", 1000, 0., 1000.);
      h_RecJetNotMatched_JetBin3_Pt->Sumw2();
      h_RecJetMatched_JetBin4_Pt = fs->make<TH1F> ("RecJetMatched_JetBin4_Pt", "RecJetMatched_JetBin4_Pt", 1000, 0., 1000.);
      h_RecJetMatched_JetBin4_Pt->Sumw2();
      h_RecJetNotMatched_JetBin4_Pt = fs->make<TH1F> ("RecJetNotMatched_JetBin4_Pt", "RecJetNotMatched_JetBin4_Pt", 1000, 0., 1000.);
      h_RecJetNotMatched_JetBin4_Pt->Sumw2();

      h_RecJetRes_Pt = fs->make<TH2F> ("RecJetRes_Pt", "RecJetRes_Pt", 100, 0., 1000., 100, 0., 3.);
      h_RecJetRes_Pt->Sumw2();
      h_RecJetRes_Eta = fs->make<TH2F> ("RecJetRes_Eta", "RecJetRes_Eta", 100, -5., 5., 100, 0., 3.);
      h_RecJetRes_Eta->Sumw2();
      h_RebJetRes_Pt = fs->make<TH2F> ("RebJetRes_Pt", "RebJetRes_Pt", 100, 0., 1000., 100, 0., 3.);
      h_RebJetRes_Pt->Sumw2();
      h_RebJetRes_Eta = fs->make<TH2F> ("RebJetRes_Eta", "RebJetRes_Eta", 100, -5., 5., 100, 0., 3.);
      h_RebJetRes_Eta->Sumw2();

      h_SmearedJetRes_Pt = fs->make<TH2F> ("SmearedJetRes_Pt", "SmearedJetRes_Pt", 100, 0., 1000., 100, 0., 3.);
      h_SmearedJetRes_Pt->Sumw2();
      h_SmearedJetRes_Eta = fs->make<TH2F> ("SmearedJetRes_Eta", "SmearedJetRes_Eta", 100, -5., 5., 100, 0., 3.);
      h_SmearedJetRes_Eta->Sumw2();

      h_HT_gen = fs->make<TH1F> ("HT_gen", "HT_gen", NbinsHT_, HTmin_, HTmax_);
      h_HT_gen->Sumw2();
      h_HT_rec = fs->make<TH1F> ("HT_rec", "HT_rec", NbinsHT_, HTmin_, HTmax_);
      h_HT_rec->Sumw2();
      h_HT_reb = fs->make<TH1F> ("HT_reb", "HT_reb", NbinsHT_, HTmin_, HTmax_);
      h_HT_reb->Sumw2();
      h_HT_smeared = fs->make<TH1F> ("HT_smeared", "HT_smeared", NbinsHT_, HTmin_, HTmax_);
      h_HT_smeared->Sumw2();

      h_HT_JetBin1_gen = fs->make<TH1F> ("HT_JetBin1_gen", "HT_JetBin1_gen", NbinsHT_, HTmin_, HTmax_);
      h_HT_JetBin1_gen->Sumw2();
      h_HT_JetBin1_rec = fs->make<TH1F> ("HT_JetBin1_rec", "HT_JetBin1_rec", NbinsHT_, HTmin_, HTmax_);
      h_HT_JetBin1_rec->Sumw2();
      h_HT_JetBin1_reb = fs->make<TH1F> ("HT_JetBin1_reb", "HT_JetBin1_reb", NbinsHT_, HTmin_, HTmax_);
      h_HT_JetBin1_reb->Sumw2();
      h_HT_JetBin1_smeared = fs->make<TH1F> ("HT_JetBin1_smeared", "HT_JetBin1_smeared", NbinsHT_, HTmin_, HTmax_);
      h_HT_JetBin1_smeared->Sumw2();

      h_HT_JetBin2_gen = fs->make<TH1F> ("HT_JetBin2_gen", "HT_JetBin2_gen", NbinsHT_, HTmin_, HTmax_);
      h_HT_JetBin2_gen->Sumw2();
      h_HT_JetBin2_rec = fs->make<TH1F> ("HT_JetBin2_rec", "HT_JetBin2_rec", NbinsHT_, HTmin_, HTmax_);
      h_HT_JetBin2_rec->Sumw2();
      h_HT_JetBin2_reb = fs->make<TH1F> ("HT_JetBin2_reb", "HT_JetBin2_reb", NbinsHT_, HTmin_, HTmax_);
      h_HT_JetBin2_reb->Sumw2();
      h_HT_JetBin2_smeared = fs->make<TH1F> ("HT_JetBin2_smeared", "HT_JetBin2_smeared", NbinsHT_, HTmin_, HTmax_);
      h_HT_JetBin2_smeared->Sumw2();

      h_HT_JetBin3_gen = fs->make<TH1F> ("HT_JetBin3_gen", "HT_JetBin3_gen", NbinsHT_, HTmin_, HTmax_);
      h_HT_JetBin3_gen->Sumw2();
      h_HT_JetBin3_rec = fs->make<TH1F> ("HT_JetBin3_rec", "HT_JetBin3_rec", NbinsHT_, HTmin_, HTmax_);
      h_HT_JetBin3_rec->Sumw2();
      h_HT_JetBin3_reb = fs->make<TH1F> ("HT_JetBin3_reb", "HT_JetBin3_reb", NbinsHT_, HTmin_, HTmax_);
      h_HT_JetBin3_reb->Sumw2();
      h_HT_JetBin3_smeared = fs->make<TH1F> ("HT_JetBin3_smeared", "HT_JetBin3_smeared", NbinsHT_, HTmin_, HTmax_);
      h_HT_JetBin3_smeared->Sumw2();

      h_HT_JetBin4_gen = fs->make<TH1F> ("HT_JetBin4_gen", "HT_JetBin4_gen", NbinsHT_, HTmin_, HTmax_);
      h_HT_JetBin4_gen->Sumw2();
      h_HT_JetBin4_rec = fs->make<TH1F> ("HT_JetBin4_rec", "HT_JetBin4_rec", NbinsHT_, HTmin_, HTmax_);
      h_HT_JetBin4_rec->Sumw2();
      h_HT_JetBin4_reb = fs->make<TH1F> ("HT_JetBin4_reb", "HT_JetBin4_reb", NbinsHT_, HTmin_, HTmax_);
      h_HT_JetBin4_reb->Sumw2();
      h_HT_JetBin4_smeared = fs->make<TH1F> ("HT_JetBin4_smeared", "HT_JetBin4_smeared", NbinsHT_, HTmin_, HTmax_);
      h_HT_JetBin4_smeared->Sumw2();

      h_HTall_gen = fs->make<TH1F> ("HTall_gen", "HTall_gen", NbinsHT_, HTmin_, HTmax_);
      h_HTall_gen->Sumw2();
      h_HTall_rec = fs->make<TH1F> ("HTall_rec", "HTall_rec", NbinsHT_, HTmin_, HTmax_);
      h_HTall_rec->Sumw2();
      h_HTall_reb = fs->make<TH1F> ("HTall_reb", "HTall_reb", NbinsHT_, HTmin_, HTmax_);
      h_HTall_reb->Sumw2();
      h_HTall_smeared = fs->make<TH1F> ("HTall_smeared", "HTall_smeared", NbinsHT_, HTmin_, HTmax_);
      h_HTall_smeared->Sumw2();

      h_HThigh_gen = fs->make<TH1F> ("HThigh_gen", "HThigh_gen", NbinsHT_, HTmin_, HTmax_);
      h_HThigh_gen->Sumw2();
      h_HThigh_rec = fs->make<TH1F> ("HThigh_rec", "HThigh_rec", NbinsHT_, HTmin_, HTmax_);
      h_HThigh_rec->Sumw2();
      h_HThigh_reb = fs->make<TH1F> ("HThigh_reb", "HThigh_reb", NbinsHT_, HTmin_, HTmax_);
      h_HThigh_reb->Sumw2();
      h_HThigh_smeared = fs->make<TH1F> ("HThigh_smeared", "HThigh_smeared", NbinsHT_, HTmin_, HTmax_);
      h_HThigh_smeared->Sumw2();

      h_MHTall_gen = fs->make<TH1F> ("MHTall_gen", "MHTall_gen", NbinsMHT_, MHTmin_, MHTmax_);
      h_MHTall_gen->Sumw2();
      h_MHTall_rec = fs->make<TH1F> ("MHTall_rec", "MHTall_rec", NbinsMHT_, MHTmin_, MHTmax_);
      h_MHTall_rec->Sumw2();
      h_MHTall_reb = fs->make<TH1F> ("MHTall_reb", "MHTall_reb", NbinsMHT_, MHTmin_, MHTmax_);
      h_MHTall_reb->Sumw2();
      h_MHTall_smeared = fs->make<TH1F> ("MHTall_smeared", "MHTall_smeared", NbinsMHT_, MHTmin_, MHTmax_);
      h_MHTall_smeared->Sumw2();

      h_MHThigh_gen = fs->make<TH1F> ("MHThigh_gen", "MHThigh_gen", NbinsMHT_, MHTmin_, MHTmax_);
      h_MHThigh_gen->Sumw2();
      h_MHThigh_rec = fs->make<TH1F> ("MHThigh_rec", "MHThigh_rec", NbinsMHT_, MHTmin_, MHTmax_);
      h_MHThigh_rec->Sumw2();
      h_MHThigh_reb = fs->make<TH1F> ("MHThigh_reb", "MHThigh_reb", NbinsMHT_, MHTmin_, MHTmax_);
      h_MHThigh_reb->Sumw2();
      h_MHThigh_smeared = fs->make<TH1F> ("MHThigh_smeared", "MHThigh_smeared", NbinsMHT_, MHTmin_, MHTmax_);
      h_MHThigh_smeared->Sumw2();

      h_deltaPhiJet1Jet2_gen = fs->make<TH1F> ("deltaPhiJet1Jet2_gen", "deltaPhiJet1Jet2_gen", 50, 0, TMath::Pi());
      h_deltaPhiJet1Jet2_gen->Sumw2();
      h_deltaPhiJet1Jet2_rec = fs->make<TH1F> ("deltaPhiJet1Jet2_rec", "deltaPhiJet1Jet2_rec", 50, 0, TMath::Pi());
      h_deltaPhiJet1Jet2_rec->Sumw2();
      h_deltaPhiJet1Jet2_reb = fs->make<TH1F> ("deltaPhiJet1Jet2_reb", "deltaPhiJet1Jet2_reb", 50, 0, TMath::Pi());
      h_deltaPhiJet1Jet2_reb->Sumw2();
      h_deltaPhiJet1Jet2_smeared = fs->make<TH1F> ("deltaPhiJet1Jet2_smeared", "deltaPhiJet1Jet2_smeared", 50, 0, TMath::Pi());
      h_deltaPhiJet1Jet2_smeared->Sumw2();

      h_fitProb = fs->make<TH1F> ("h_fitProb", "h_fitProb", 100, 0., 1.);
      h_fitProb->Sumw2();
      h_weight = fs->make<TH1F> ("h_weight", "h_weight", 70, -1., 6.);
      h_weight->Sumw2();
      h_weightedWeight = fs->make<TH1F> ("h_weightedWeight", "h_weightedWeight", 70, -1., 6.);
      h_weightedWeight->Sumw2();
   }

   //// get rebalance correction histo
   if( useRebalanceCorrectionFactors_ ) {
      TFile *f_rebCorr = new TFile(RebalanceCorrectionFile_.c_str(), "READ", "", 0);
      h_RebCorrectionFactor = (TH1F*) f_rebCorr->FindObjectAny("RebCorrection_vsReco_px");
      h_RebCorrectionFactor_b = (TH1F*) f_rebCorr->FindObjectAny("RebCorrection_vsReco_b_px");
      h_2DRebCorrectionFactor = (TH2F*) f_rebCorr->FindObjectAny("RebCorrection_vsReco");
      h_2DRebCorrectionFactor_b = (TH2F*) f_rebCorr->FindObjectAny("RebCorrection_vsReco_b");
      //// Do projections for each x-bin
      for (int ii = 1; ii <= h_2DRebCorrectionFactor_b->GetXaxis()->GetNbins(); ++ii){
         TH1D* tmp_py = new TH1D(*h_2DRebCorrectionFactor_b->ProjectionY("py", ii, ii));
         h_2DRebCorrectionFactor_b_py.push_back(tmp_py);
      }
      for (int ii = 1; ii <= h_2DRebCorrectionFactor->GetXaxis()->GetNbins(); ++ii){
         TH1D* tmp_py = new TH1D(*h_2DRebCorrectionFactor->ProjectionY("py", ii, ii));
         h_2DRebCorrectionFactor_py.push_back(tmp_py);
      }
   }

   //------------------------------------------------------
   //---------------  BTag Corrections  -------------------
   //------------------------------------------------------
   //
   cout << "will read b-tag information from file" << BTagEfficiencyFile_.c_str() << endl;
   if( useBTagEfficiencyFactors_ ){
       char hname[100];
       TFile *f_btags = new TFile(BTagEfficiencyFile_.c_str(), "READ", "", 0);
       // Get histogram for each eta bin:
       for( int e_eta = 0; e_eta < 12; ++e_eta ){

            // get btag efficiency correction histo
            cout << "Reading efficiency factors from file..." << endl;
            sprintf(hname, "BTrue_BTag_vs_RecoPt_Eta%i;1", e_eta);
            if(f_btags->FindObjectAny(hname)){
                BTagEfficiencyFactors.push_back((TH1F*) f_btags->FindObjectAny(hname));
            }

            // get bmistag efficiency correction histo
            cout << "Reading efficiency factors from file..." << endl;
            sprintf(hname, "no_BTrue_BTag_vs_RecoPt_Eta%i;1", e_eta);
            if(f_btags->FindObjectAny(hname)){
                BMisTagEfficiencyFactors.push_back((TH1F*) f_btags->FindObjectAny(hname));
            }

            // get nbtrue efficiency correction histo
            cout << "Reading NBTrue from file..." << endl;
            sprintf(hname, "BTrue_vs_RecoPt_Eta%i;1", e_eta);
            if(f_btags->FindObjectAny(hname)){
                NBTrues.push_back((TH1F*) f_btags->FindObjectAny(hname));
            }

            // get nbtag efficiency correction histo
            cout << "Reading NBTag from file..." << endl;
            sprintf(hname, "BTag_vs_RecoPt_Eta%i;1", e_eta);
            if(f_btags->FindObjectAny(hname)){
                NBTags.push_back((TH1F*) f_btags->FindObjectAny(hname));
            }

            // get n no btag efficiency correction histo
            cout << "Reading NnoBTag from file..." << endl;
            sprintf(hname, "no_BTag_vs_RecoPt_Eta%i;1", e_eta);
            if(f_btags->FindObjectAny(hname)){
                NnoBTags.push_back((TH1F*) f_btags->FindObjectAny(hname));
            }

       }
       cout << "Read efficiency factors from file." << endl;

       //h_trueb_SmearedPt.resize(EtaBinEdges.size() - 1);
       //h_no_trueb_SmearedPt.resize(EtaBinEdges.size() - 1);
       //h_trueb_btag_SmearedPt.resize(EtaBinEdges.size() - 1);
       //h_no_trueb_btag_SmearedPt.resize(EtaBinEdges.size() - 1);
       //h_btag_SmearedPt.resize(EtaBinEdges.size() - 1);
       //h_btag_old_SmearedPt.resize(EtaBinEdges.size() - 1);
   }

   // define output tree
   PredictionTree = fs->make<TTree> ("QCDPrediction", "QCDPrediction", 0);
   PredictionTree->SetAutoSave(10000000000);
   PredictionTree->SetAutoFlush(1000000);
   cout << "Defined output tree." << endl;

   // set branches for output tree
   PredictionTree->Branch("NVtx", &vtxN);
   PredictionTree->Branch("Ntries",&Ntries_pred);
   PredictionTree->Branch("NJets",&Njets_pred);
   PredictionTree->Branch("BTags",&BTags_pred);
   PredictionTree->Branch("Weight",&weight);
   PredictionTree->Branch("HT_seed", &HT_seed);
   PredictionTree->Branch("HT", &HT_pred);
   PredictionTree->Branch("MHT", &MHT_pred);
   PredictionTree->Branch("Jet1Pt", &Jet1Pt_pred);
   PredictionTree->Branch("Jet2Pt", &Jet2Pt_pred);
   PredictionTree->Branch("Jet3Pt", &Jet3Pt_pred);
   PredictionTree->Branch("Jet4Pt", &Jet4Pt_pred);
   PredictionTree->Branch("Jet1Eta", &Jet1Eta_pred);
   PredictionTree->Branch("Jet2Eta", &Jet2Eta_pred);
   PredictionTree->Branch("Jet3Eta", &Jet3Eta_pred);
   PredictionTree->Branch("Jet4Eta", &Jet4Eta_pred);
   PredictionTree->Branch("DeltaPhi1", &DeltaPhi1_pred);
   PredictionTree->Branch("DeltaPhi2", &DeltaPhi2_pred);
   PredictionTree->Branch("DeltaPhi3", &DeltaPhi3_pred);
   PredictionTree->Branch("minDeltaPhiN", &DeltaPhiMinN_pred);
   cout << "Set branches for output tree." << endl;

}
//--------------------------------------------------------------------------



// ------------ method called once each job just after ending the event loop  ------------
void QCDBkgRS::endJob()
{

}

//define this as a plug-in
DEFINE_FWK_MODULE(QCDBkgRS);
