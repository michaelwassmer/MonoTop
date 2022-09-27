// -*- C++ -*-
//
// Package:    MonoTop/FlagFilter
// Class:      FlagFilter
//
/**\class FlagFilter FlagFilter.cc MonoTop/FlagFilter/plugins/FlagFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michael Wassmer
//         Created:  Sat, 27 Jun 2020 08:05:09 GMT
//
//

// system include files
#include <memory>
#include <iostream>
#include <string>
#include <regex>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

//
// class declaration
//
class FlagFilter : public edm::stream::EDFilter<> {
   public:
    explicit FlagFilter(const edm::ParameterSet&);
    ~FlagFilter();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
    virtual void beginStream(edm::StreamID) override;
    virtual bool filter(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    // token for the filter bits
    const edm::EDGetTokenT< edm::TriggerResults > filterBitsToken;
    // vector for the filter names
    const std::vector< std::string > filterNames;
    // flag for tagging mode
    const bool taggingMode;
    // AND or OR mode
    const bool OR_mode;
    // single filter decisions
    std::vector< std::string > filters;
    std::vector< int > indices;

    bool run_start;

    virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    // virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    // virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
    // virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

    // ----------member data ---------------------------
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
FlagFilter::FlagFilter(const edm::ParameterSet& iConfig) :
    filterBitsToken{consumes< edm::TriggerResults >(iConfig.getParameter< edm::InputTag >("filterData"))},  // filter bits
    filterNames{iConfig.getParameter< std::vector< std::string > >("filterNames")},                          // names of filters
    taggingMode{iConfig.getParameter< bool >("taggingMode")},
    OR_mode{iConfig.getParameter< bool >("OR_mode")}
{
    // now do what ever initialization is needed
    if (taggingMode) {
        produces< std::vector< bool > >("decisions");
        produces< std::vector< std::string > >("filters");
    }
}

FlagFilter::~FlagFilter()
{
    // do anything here that needs to be done at destruction time
    // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called on each new Event  ------------
bool FlagFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // get filter information
    edm::Handle< edm::TriggerResults > filterData;
    iEvent.getByToken(filterBitsToken, filterData);

    // get filter names and indices in first event of currently processed run
    if(run_start) {
        const edm::TriggerNames& names = iEvent.triggerNames(*filterData);
        for (size_t i = 0; i < filterData->size(); i++) {
            // get filter name
            std::string name = names.triggerName(i);
            // check if the filter at hand is in the list of desired filters
            bool filterfound = false;
            bool regexfound = false;
            for(auto filterName : filterNames) {
                filterfound = filterfound or (filterName == name);
                if (filterfound) break;
                regexfound = regexfound or (std::regex_match(name, std::regex(filterName)));
                if (regexfound) break;
            }
            if (not (filterfound or regexfound)) continue;
            // if the filter is in the list, check if the event passes the desired filters
            filters.push_back(name);
            indices.push_back(i);
        }
        run_start = false;
    }
    // single filter decisions
    std::vector< bool > decisions;
    // final decisions
    bool decision_AND = true;
    bool decision_OR = false;
    // loop over desired filters
    for (size_t i = 0; i < filters.size(); i++) {
        // get filter index
        int& index = indices.at(i);
        // check if the event passes the desired filter
        bool filter_decision = filterData->accept(index);
        if ((not OR_mode) and (not filter_decision)) return false;
        decisions.push_back(filter_decision);
        // update final decisions
        decision_AND = decision_AND and filter_decision;
        decision_OR = decision_OR or filter_decision;
    }
    // determine final event decision
    if (OR_mode) {
        if (not decision_OR) return false;
    }
    else {
        if (not decision_AND) return false;
    }
    if (taggingMode) {
        std::unique_ptr< std::vector< bool > > decisions_ptr = std::make_unique< std::vector< bool > >(decisions);
        std::unique_ptr< std::vector< std::string > > filters_ptr = std::make_unique< std::vector< std::string > >(filters);
        iEvent.put(std::move(decisions_ptr), "decisions");
        iEvent.put(std::move(filters_ptr), "filters");
    }
    // if the event passes the desired filters, keep the event
    return true;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void FlagFilter::beginStream(edm::StreamID) {}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void FlagFilter::endStream() {}

// ------------ method called when starting to processes a run  ------------

void
FlagFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{
   run_start = true;
   filters.clear();
   indices.clear();
}


// ------------ method called when ending the processing of a run  ------------
/*
void
FlagFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
FlagFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
FlagFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void FlagFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    // The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}
// define this as a plug-in
DEFINE_FWK_MODULE(FlagFilter);
