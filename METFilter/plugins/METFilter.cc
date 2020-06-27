// -*- C++ -*-
//
// Package:    MonoTop/METFilter
// Class:      METFilter
// 
/**\class METFilter METFilter.cc MonoTop/METFilter/plugins/METFilter.cc

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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include <string>
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

//
// class declaration
//

class METFilter : public edm::stream::EDFilter<> {
   public:
      explicit METFilter(const edm::ParameterSet&);
      ~METFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
      
      edm::EDGetTokenT< edm::TriggerResults > filterBitsToken;
      std::vector<std::string> filterNames;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

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
METFilter::METFilter(const edm::ParameterSet& iConfig) :
  filterBitsToken{consumes< edm::TriggerResults >(iConfig.getParameter< edm::InputTag >("filterData"))},
  filterNames{iConfig.getParameter< std::vector< std::string > >("filterNames")}
{
   //now do what ever initialization is needed

}


METFilter::~METFilter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
METFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle< edm::TriggerResults > filterData;
    iEvent.getByToken(filterBitsToken, filterData);
    const edm::TriggerNames& names = iEvent.triggerNames(*filterData);
    for (size_t i = 0; i < filterData->size(); i++) {
        std::string name   = names.triggerName(i);
        if(std::find(filterNames.begin(), filterNames.end(), name) == filterNames.end()) continue;
        bool filter_decision = filterData->accept(i);
        if(not filter_decision) return false;
    }
    return true;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
METFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
METFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
METFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
METFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
METFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
METFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
METFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(METFilter);
