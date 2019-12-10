#ifndef CSVHelperSkim_h__
#define CSVHelperSkim_h__

#include <string>
#include <vector>

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "TString.h"

class CSVHelperSkim {
   public:
    // standard constructor
    CSVHelperSkim();
    // destructor
    ~CSVHelperSkim();

    enum class CSVwp { Tight, Medium, Loose, None };

    static float GetJetCSV(const pat::Jet &jet, const std::string taggername);
    static bool  PassesCSV(const pat::Jet &iJet, std::string taggername, const CSVwp iCSVworkingPoint, std::string dataEra);
    static float GetWP(std::string dataEra, const CSVwp iCSVworkingPoint, std::string taggername);

   private:
};

#endif
