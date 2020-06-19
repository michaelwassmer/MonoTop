// -*- C++ -*-
//
// Package:    MonoTop/BoostedAnalyzer
// Class:      MonoTopSkim
//
/**\class MonoTopSkim MonoTopSkim.cc
 MonoTop/BoostedProducer/plugins/MonoTopSkim.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michael Wassmer
//
//

// system include files
#include <iostream>
#include <memory>

// user include files
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "TVector2.h"

//
// class declaration
//

class MonoTopSkim : public edm::EDFilter {
   public:
    explicit MonoTopSkim(const edm::ParameterSet &);
    ~MonoTopSkim();

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

   private:
    virtual void beginJob() override;
    virtual bool filter(edm::Event &, const edm::EventSetup &) override;
    virtual void endJob() override;

    // ----------member data ---------------------------

    // data access tokens
    edm::EDGetTokenT< pat::ElectronCollection > EDMElectronsToken;  // electrons
    edm::EDGetTokenT< pat::MuonCollection >     EDMMuonsToken;      // muons
    edm::EDGetTokenT< pat::PhotonCollection >   EDMPhotonsToken;    // photons
    edm::EDGetTokenT< pat::JetCollection >      EDMAK4JetsToken;    // AK4 jets
    edm::EDGetTokenT< pat::JetCollection >      EDMAK8JetsToken;    // AK8 jets
    edm::EDGetTokenT< pat::JetCollection >      EDMAK15JetsToken_;  // AK15 jets
    edm::EDGetTokenT< std::vector< pat::MET > > EDMMETToken;        // PF MET
    edm::EDGetTokenT< std::vector< pat::MET > > EDMPuppiMETToken;   // PUPPI MET
    edm::EDGetTokenT< reco::VertexCollection >  EDMVertexToken;     // vertex
    // edm::EDGetTokenT< double >                  EDMRhoToken;     //  pileup density

    const int    minJetsAK4_;
    const int    minJetsAK8_;
    const int    minJetsAK15_;
    const int    maxJetsAK4_;
    const int    maxJetsAK8_;
    const int    maxJetsAK15_;
    const double AK4jetPtMin_;
    const double AK4jetEtaMax_;
    const double AK8jetPtMin_;
    const double AK8jetEtaMax_;
    const double muonPtMin_;
    const double muonEtaMax_;
    const double electronPtMin_;
    const double electronEtaMax_;
    const double photonPtMin_;
    const double photonEtaMax_;
    const double metPtMin_;

    const bool        isData;
    const std::string era;
};

//
// constructors and destructor
//
MonoTopSkim::MonoTopSkim(const edm::ParameterSet &iConfig) :
    // now do what ever initialization is needed
    EDMElectronsToken{consumes< pat::ElectronCollection >(iConfig.getParameter< edm::InputTag >("electrons"))},
    EDMMuonsToken{consumes< pat::MuonCollection >(iConfig.getParameter< edm::InputTag >("muons"))},
    EDMPhotonsToken{consumes< std::vector< pat::Photon > >(iConfig.getParameter< edm::InputTag >("photons"))},
    EDMAK4JetsToken{consumes< pat::JetCollection >(iConfig.getParameter< edm::InputTag >("AK4jets"))},
    EDMAK8JetsToken{consumes< pat::JetCollection >(iConfig.getParameter< edm::InputTag >("AK8jets"))},
    EDMAK15JetsToken_{consumes< pat::JetCollection >(iConfig.getParameter< edm::InputTag >("AK15jets"))},
    EDMMETToken{consumes< std::vector< pat::MET > >(iConfig.getParameter< edm::InputTag >("met"))},
    EDMPuppiMETToken{consumes< std::vector< pat::MET > >(iConfig.getParameter< edm::InputTag >("met_puppi"))},
    EDMVertexToken{consumes< reco::VertexCollection >(iConfig.getParameter< edm::InputTag >("vertices"))},
    // EDMRhoToken{consumes< double >(iConfig.getParameter< edm::InputTag >("rho"))},
    minJetsAK4_{iConfig.getParameter< int >("minJetsAK4")},
    minJetsAK8_{iConfig.getParameter< int >("minJetsAK8")},
    minJetsAK15_{iConfig.getParameter< int >("minJetsAK15")},
    maxJetsAK4_{iConfig.getParameter< int >("maxJetsAK4")},
    maxJetsAK8_{iConfig.getParameter< int >("maxJetsAK8")},
    maxJetsAK15_{iConfig.getParameter< int >("maxJetsAK15")},
    AK4jetPtMin_{iConfig.getParameter< double >("AK4jetPtMin")},
    AK4jetEtaMax_{iConfig.getParameter< double >("AK4jetEtaMax")},
    AK8jetPtMin_{iConfig.getParameter< double >("AK8jetPtMin")},
    AK8jetEtaMax_{iConfig.getParameter< double >("AK8jetEtaMax")},
    muonPtMin_{iConfig.getParameter< double >("muonPtMin")},
    muonEtaMax_{iConfig.getParameter< double >("muonEtaMax")},
    electronPtMin_{iConfig.getParameter< double >("electronPtMin")},
    electronEtaMax_{iConfig.getParameter< double >("electronEtaMax")},
    photonPtMin_{iConfig.getParameter< double >("photonPtMin")},
    photonEtaMax_{iConfig.getParameter< double >("photonEtaMax")},
    metPtMin_{iConfig.getParameter< double >("metPtMin")},
    isData{iConfig.getParameter< bool >("isData")},
    era{iConfig.getParameter< std::string >("era")}
{
}

MonoTopSkim::~MonoTopSkim()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called on each new Event  ------------
bool MonoTopSkim::filter(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
    // get slimmedMETs (PFMET)
    edm::Handle< std::vector< pat::MET > > hMETs;
    iEvent.getByToken(EDMMETToken, hMETs);

    // get slimmedMETs Puppi (PUPPIMET)
    edm::Handle< std::vector< pat::MET > > hPuppiMETs;
    iEvent.getByToken(EDMPuppiMETToken, hPuppiMETs);

    // get PF and Puppi MET 4-vector also considering combined JES and JER variations
    auto met          = hMETs->at(0).corP4(pat::MET::Type1);
    auto met_jes_up   = hMETs->at(0).shiftedP4(pat::MET::JetEnUp, pat::MET::Type1);
    auto met_jes_down = hMETs->at(0).shiftedP4(pat::MET::JetEnDown, pat::MET::Type1);
    auto met_jer_up   = hMETs->at(0).shiftedP4(pat::MET::JetResUp, pat::MET::Type1);
    auto met_jer_down = hMETs->at(0).shiftedP4(pat::MET::JetResDown, pat::MET::Type1);

    auto met_puppi          = hPuppiMETs->at(0).corP4(pat::MET::Type1);
    auto met_puppi_jes_up   = hPuppiMETs->at(0).shiftedP4(pat::MET::JetEnUp, pat::MET::Type1);
    auto met_puppi_jes_down = hPuppiMETs->at(0).shiftedP4(pat::MET::JetEnDown, pat::MET::Type1);
    auto met_puppi_jer_up   = hPuppiMETs->at(0).shiftedP4(pat::MET::JetResUp, pat::MET::Type1);
    auto met_puppi_jer_down = hPuppiMETs->at(0).shiftedP4(pat::MET::JetResDown, pat::MET::Type1);

    // determine the maximum met within the nominal value and the variations
    auto met_max       = std::max({met.pt(), met_jes_up.pt(), met_jes_down.pt(), met_jer_up.pt(), met_jer_down.pt()});
    auto met_puppi_max = std::max({met_puppi.pt(), met_puppi_jes_up.pt(), met_puppi_jes_down.pt(), met_puppi_jer_up.pt(), met_puppi_jer_down.pt()});

    // std::cout << "MET: " << met.pt() << std::endl;
    // std::cout << "MET up: " << met_up.pt() << std::endl;
    // std::cout << "MET down: " << met_down.pt() << std::endl;
    // std::cout << "Puppi MET: " << met_puppi.pt() << std::endl;
    // std::cout << "Puppi MET jes up: " << met_puppi_jes_up.pt() << std::endl;
    // std::cout << "Puppi MET jes down: " << met_puppi_jes_down.pt() << std::endl;
    // std::cout << "Puppi MET jer up: " << met_puppi_jer_up.pt() << std::endl;
    // std::cout << "Puppi MET jer down: " << met_puppi_jer_down.pt() << std::endl;

    // check if we have sizeable MET (200 GeV) in the event and if so, keep the event (hadronic analysis)
    bool met_criterium = (met_max >= metPtMin_) || (met_puppi_max >= metPtMin_);

    // get AK4 jets
    edm::Handle< pat::JetCollection > ak4Jets;
    iEvent.getByToken(EDMAK4JetsToken, ak4Jets);

    // get AK8 jets
    edm::Handle< pat::JetCollection > ak8Jets;
    iEvent.getByToken(EDMAK8JetsToken, ak8Jets);

    // get AK15 jets
    edm::Handle< pat::JetCollection > ak15Jets;
    iEvent.getByToken(EDMAK15JetsToken_, ak15Jets);

    int n_ak4jets  = ak4Jets->size();
    int n_ak8jets  = ak8Jets->size();
    int n_ak15jets = ak15Jets->size();

    // std::cout << "Number of AK4 jets: " << n_ak4jets << std::endl;
    // std::cout << "Number of AK8 jets: " << n_ak8jets << std::endl;
    // std::cout << "Number of AK15 jets: " << n_ak15jets << std::endl;

    auto leading_jet_pt = n_ak4jets > 0 ? ak4Jets->at(0).pt() : 0.;

    // want at least one fat jet for hadronic monotop regions (hadronic analysis)
    bool jet_criterium = ((n_ak8jets >= minJetsAK8_) || (n_ak15jets >= minJetsAK15_) || (leading_jet_pt >= 100.));

    // if met criterium and fat jet criterium is fulfilled, keep the event (hadronic analysis)
    if (met_criterium && jet_criterium) return true;

    // get slimmedElectrons
    edm::Handle< pat::ElectronCollection > hElectrons;
    iEvent.getByToken(EDMElectronsToken, hElectrons);

    // select those electrons satifsying pt and eta cuts and loose cut-based
    // electron ID
    std::vector< pat::Electron > selectedElectrons = *hElectrons;
    selectedElectrons.erase(
        std::remove_if(selectedElectrons.begin(), selectedElectrons.end(),
                       [&](pat::Electron ele) {
                           return (ele.pt() < electronPtMin_ || fabs(ele.eta()) > electronEtaMax_ || !ele.electronID("cutBasedElectronID-Fall17-94X-V2-veto"));
                       }),
        selectedElectrons.end());

    // get slimmedMuons
    edm::Handle< pat::MuonCollection > hMuons;
    iEvent.getByToken(EDMMuonsToken, hMuons);

    // select those muons satisfying pt and eta cuts and loose cut-base muon ID and loose PF isolation
    std::vector< pat::Muon > selectedMuons = *hMuons;
    selectedMuons.erase(
        std::remove_if(selectedMuons.begin(), selectedMuons.end(),
                       [&](pat::Muon mu) {
                           return (mu.pt() < muonPtMin_ || fabs(mu.eta()) > muonEtaMax_ || !muon::isLooseMuon(mu) || !mu.passed(pat::Muon::PFIsoLoose));
                       }),
        selectedMuons.end());

    // get slimmedPhotons
    edm::Handle< pat::PhotonCollection > hPhotons;
    iEvent.getByToken(EDMPhotonsToken, hPhotons);

    // select those photons satisfying pt and eta cuts and loose cut-based photon
    // ID
    std::vector< pat::Photon > selectedPhotons = *hPhotons;
    selectedPhotons.erase(
        std::remove_if(
            selectedPhotons.begin(), selectedPhotons.end(),
            [&](pat::Photon ph) { return (ph.pt() < photonPtMin_ || fabs(ph.eta()) > photonEtaMax_ || !ph.photonID("cutBasedPhotonID-Fall17-94X-V2-loose")); }),
        selectedPhotons.end());

    // std::cout << "Number of selected electrons: " << selectedElectrons.size() << std::endl;
    // std::cout << "Number of selected muons: " << selectedMuons.size() << std::endl;
    // std::cout << "Number of selected photons: " << selectedPhotons.size() << std::endl;

    // calculate very conservative hadronic recoil
    auto hadr_recoil          = met.pt();
    auto hadr_recoil_jes_up   = met_jes_up.pt();
    auto hadr_recoil_jes_down = met_jes_down.pt();
    auto hadr_recoil_jer_up   = met_jer_up.pt();
    auto hadr_recoil_jer_down = met_jer_down.pt();

    auto hadr_recoil_puppi          = met_puppi.pt();
    auto hadr_recoil_puppi_jes_up   = met_puppi_jes_up.pt();
    auto hadr_recoil_puppi_jes_down = met_puppi_jes_down.pt();
    auto hadr_recoil_puppi_jer_up   = met_puppi_jer_up.pt();
    auto hadr_recoil_puppi_jer_down = met_puppi_jer_down.pt();

    for (const auto &ele : selectedElectrons) {
        hadr_recoil += ele.pt();
        hadr_recoil_jes_up += ele.pt();
        hadr_recoil_jes_down += ele.pt();
        hadr_recoil_jer_up += ele.pt();
        hadr_recoil_jer_down += ele.pt();

        hadr_recoil_puppi += ele.pt();
        hadr_recoil_puppi_jes_up += ele.pt();
        hadr_recoil_puppi_jes_down += ele.pt();
        hadr_recoil_puppi_jer_up += ele.pt();
        hadr_recoil_puppi_jer_down += ele.pt();
    }
    for (const auto &mu : selectedMuons) {
        hadr_recoil += mu.pt();
        hadr_recoil_jes_up += mu.pt();
        hadr_recoil_jes_down += mu.pt();
        hadr_recoil_jer_up += mu.pt();
        hadr_recoil_jer_down += mu.pt();

        hadr_recoil_puppi += mu.pt();
        hadr_recoil_puppi_jes_up += mu.pt();
        hadr_recoil_puppi_jes_down += mu.pt();
        hadr_recoil_puppi_jer_up += mu.pt();
        hadr_recoil_puppi_jer_down += mu.pt();
    }
    for (const auto &ph : selectedPhotons) {
        hadr_recoil += ph.pt();
        hadr_recoil_jes_up += ph.pt();
        hadr_recoil_jes_down += ph.pt();
        hadr_recoil_jer_up += ph.pt();
        hadr_recoil_jer_down += ph.pt();

        hadr_recoil_puppi += ph.pt();
        hadr_recoil_puppi_jes_up += ph.pt();
        hadr_recoil_puppi_jes_down += ph.pt();
        hadr_recoil_puppi_jer_up += ph.pt();
        hadr_recoil_puppi_jer_down += ph.pt();
    }

    // std::cout << "Hadronic recoil: " << hadr_recoil.pt() << std::endl;
    // std::cout << "Puppi Hadronic recoil: " << hadr_recoil_puppi.pt() << std::endl;

    // determine the maximum hadronic recoil within the nominal value and the JES/JER variations
    auto hadr_recoil_max = std::max({hadr_recoil, hadr_recoil_jes_up, hadr_recoil_jes_down, hadr_recoil_jer_up, hadr_recoil_jer_down});
    auto hadr_recoil_puppi_max =
        std::max({hadr_recoil_puppi, hadr_recoil_puppi_jes_up, hadr_recoil_puppi_jes_down, hadr_recoil_puppi_jer_up, hadr_recoil_puppi_jer_down});

    // check if we have sizeable hadronic recoil (200 GeV) in the event (hadronic analysis)
    bool recoil_criterium = (hadr_recoil_max >= metPtMin_) || (hadr_recoil_puppi_max >= metPtMin_);

    // keep the event if recoil and fatjet criteria are fulfilled (hadronic analysis)
    if (recoil_criterium && jet_criterium) return true;

    // -----------------------------------------------------------------------------------------------------------
    // Here begins the skimming part for the leptonic analysis

    // get slimmedVertices
    edm::Handle< reco::VertexCollection > hVertices;
    iEvent.getByToken(EDMVertexToken, hVertices);
    auto vertex = hVertices->empty() ? reco::Vertex() : hVertices->at(0);

    // reset lepton collections
    selectedElectrons = *hElectrons;
    selectedMuons     = *hMuons;

    // for leptonic monotop events, find tight electrons
    selectedElectrons.erase(std::remove_if(selectedElectrons.begin(), selectedElectrons.end(),
                                           [&](pat::Electron ele) {
                                               return (ele.pt() < (electronPtMin_ + 10.) || fabs(ele.eta()) > electronEtaMax_ ||
                                                       !ele.electronID("cutBasedElectronID-Fall17-94X-V2-tight"));
                                           }),
                            selectedElectrons.end());

    // for leptonic monotop events, find tight muons
    selectedMuons.erase(std::remove_if(selectedMuons.begin(), selectedMuons.end(),
                                       [&](pat::Muon mu) {
                                           return (mu.pt() < (muonPtMin_ + 10.) || fabs(mu.eta()) > muonEtaMax_ || !muon::isTightMuon(mu, vertex) ||
                                                   !mu.passed(pat::Muon::PFIsoTight));
                                       }),
                        selectedMuons.end());

    // number of tight leptons (electrons and muons)
    int n_electrons = selectedElectrons.size();
    int n_muons     = selectedMuons.size();
    int n_leptons   = n_electrons + n_muons;

    // int n_harder_jets = std::count_if(ak4Jets->begin(), ak4Jets->end(), [&](pat::Jet jet) { return (jet.pt() >= 40.); });

    std::vector< float > m_transverse;
    for (const auto &ele : selectedElectrons) {
        auto  cos_dphi_met_lep = TMath::Cos(fabs(TVector2::Phi_mpi_pi(met.phi() - ele.phi())));
        float m_transv         = TMath::Sqrt(2 * ele.pt() * met.pt() * (1 - cos_dphi_met_lep));
        m_transverse.push_back(m_transv);
    }
    for (const auto &mu : selectedMuons) {
        auto  cos_dphi_met_lep = TMath::Cos(fabs(TVector2::Phi_mpi_pi(met.phi() - mu.phi())));
        float m_transv         = TMath::Sqrt(2 * mu.pt() * met.pt() * (1 - cos_dphi_met_lep));
        m_transverse.push_back(m_transv);
    }
    for (const auto &ele : selectedElectrons) {
        auto  cos_dphi_met_lep = TMath::Cos(fabs(TVector2::Phi_mpi_pi(met_puppi.phi() - ele.phi())));
        float m_transv         = TMath::Sqrt(2 * ele.pt() * met_puppi.pt() * (1 - cos_dphi_met_lep));
        m_transverse.push_back(m_transv);
    }
    for (const auto &mu : selectedMuons) {
        auto  cos_dphi_met_lep = TMath::Cos(fabs(TVector2::Phi_mpi_pi(met_puppi.phi() - mu.phi())));
        float m_transv         = TMath::Sqrt(2 * mu.pt() * met_puppi.pt() * (1 - cos_dphi_met_lep));
        m_transverse.push_back(m_transv);
    }
    // check if at least one of the transverse masses is above 20 GeV
    bool m_t_criterium = (std::count_if(m_transverse.begin(), m_transverse.end(), [&](float m_t) { return (m_t >= 20.); })) > 0;

    // leading lepton pts
    // auto leading_ele    = n_electrons > 0 ? selectedElectrons.at(0).p4() : math::XYZTLorentzVector(0., 0., 0., 0.);
    // auto leading_muon   = n_muons > 0 ? selectedMuons.at(0).p4() : math::XYZTLorentzVector(0., 0., 0., 0.);
    // auto leading_lepton = leading_ele.pt() > leading_muon.pt() ? leading_ele : leading_muon;

    // auto cos_dphi_met_lep       = TMath::Cos(fabs(TVector2::Phi_mpi_pi(met.phi() - leading_lepton.phi())));
    // auto m_W_transv             = TMath::Sqrt(2 * leading_lepton.pt() * met.pt() * (1 - cos_dphi_met_lep));
    // auto cos_dphi_met_lep_puppi = TMath::Cos(fabs(TVector2::Phi_mpi_pi(met_puppi.phi() - leading_lepton.phi())));
    // auto m_W_transv_puppi       = TMath::Sqrt(2 * leading_lepton.pt() * met_puppi.pt() * (1 - cos_dphi_met_lep_puppi));

    // require at least one tight lepton, at least one jet with 50 GeV and at least around 100 GeV MET
    bool lepton_jet_met_criterium =
        (n_leptons >= 1) && (n_ak4jets >= minJetsAK4_) && (leading_jet_pt >= 50.) && (met_max >= 90. || met_puppi_max >= 90.) && m_t_criterium;

    if (lepton_jet_met_criterium) return true;

    return false;
}

// ------------ method called once each job just before starting event loop
// ------------
void MonoTopSkim::beginJob() {}

// ------------ method called once each job just after ending the event loop
// ------------
void MonoTopSkim::endJob() {}

// ------------ method called when starting to processes a run  ------------
/*
void
MonoTopSkim::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
MonoTopSkim::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block
// ------------
/*
void
MonoTopSkim::beginLuminosityBlock(edm::LuminosityBlock const&,
edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block
// ------------
/*
void
MonoTopSkim::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup
const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------
void MonoTopSkim::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
    // The following says we do not know what parameters are allowed so do no
    // validation
    // Please change this to state exactly what you do use, even if it is no
    // parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}
// define this as a plug-in
DEFINE_FWK_MODULE(MonoTopSkim);
