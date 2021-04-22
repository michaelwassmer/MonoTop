// -*- C++ -*-
//
// Package:    MonoTop/EventCounter
// Class:      EventCounter
//
/**\class EventCounter EventCounter.cc MonoTop/EventCounter/plugins/EventCounter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michael Wassmer
//         Created:  Fri, 25 Sep 2020 13:41:53 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "TFile.h"
#include "TH1D.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class EventCounter : public edm::one::EDAnalyzer< edm::one::SharedResources > {
   public:
    explicit EventCounter(const edm::ParameterSet&);
    ~EventCounter();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------

    TH1D* entry_histo;
    TH1D* lhe_weight_histo;
    TH1D* gen_weight_histo;

    /** LHE data access token **/
    edm::EDGetTokenT< LHEEventProduct > lheInfoToken;
    /** gen info data access token **/
    edm::EDGetTokenT< GenEventInfoProduct > genInfoToken;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
EventCounter::EventCounter(const edm::ParameterSet& iConfig) :
    lheInfoToken{consumes< LHEEventProduct >(iConfig.getParameter< edm::InputTag >("lheInfo"))},
    genInfoToken{consumes< GenEventInfoProduct >(iConfig.getParameter< edm::InputTag >("genInfo"))}

{
    // now do what ever initialization is needed
    entry_histo      = new TH1D("entry_histo", "entry_histo", 1, 0, 2);
    lhe_weight_histo = new TH1D("lhe_weight_histo", "lhe_weight_histo", 1, 0, 2);
    gen_weight_histo = new TH1D("gen_weight_histo", "gen_weight_histo", 1, 0, 2);
}

EventCounter::~EventCounter()
{
    TFile* output_file = new TFile("EventCounterHisto.root", "RECREATE");
    output_file->WriteTObject(entry_histo);
    output_file->WriteTObject(lhe_weight_histo);
    output_file->WriteTObject(gen_weight_histo);
    output_file->Close();
}

//
// member functions
//

// ------------ method called for each event  ------------
void EventCounter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    edm::Handle< GenEventInfoProduct > h_genInfo;
    edm::Handle< LHEEventProduct >     h_lheInfo;
    iEvent.getByToken(genInfoToken, h_genInfo);
    iEvent.getByToken(lheInfoToken, h_lheInfo);
    auto         lhe_product = h_lheInfo.product();
    auto         gen_product = h_genInfo.product();
    const double lhe_weight  = lhe_product->originalXWGTUP();
    const double gen_weight  = gen_product->weight();
    entry_histo->Fill(1.0);
    lhe_weight_histo->Fill(1.0, lhe_weight);
    gen_weight_histo->Fill(1.0, gen_weight);
}

// ------------ method called once each job just before starting event loop  ------------
void EventCounter::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void EventCounter::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EventCounter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    // The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(EventCounter);
