
#include "MonoTop/MonoTopSkim/interface/CSVHelperSkim.h"

CSVHelperSkim::CSVHelperSkim() {}

float CSVHelperSkim::GetJetCSV(const pat::Jet &jet, const std::string taggername)
{
    float defaultFailure = -.1;
    float bTagVal        = 0;
    if (taggername == "DeepCSV") { bTagVal = jet.bDiscriminator("pfDeepCSVJetTags:probb") + jet.bDiscriminator("pfDeepCSVJetTags:probbb"); }
    else if (taggername == "DeepJet") {
        bTagVal = jet.bDiscriminator("pfDeepFlavourJetTags:probb") + jet.bDiscriminator("pfDeepFlavourJetTags:probbb") +
                  jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
    }
    else if (taggername == "CSVv2" or TString(taggername).BeginsWith("pfDeep")) {
        bTagVal = jet.bDiscriminator(taggername);
    }
    else {
        throw cms::Exception("CSVHelper: Invalid taggername ")
            << "Taggername '" << taggername << "' not recognized, only DeepCSV/DeepJet/CSVv2 possible" << std::endl;
        bTagVal = defaultFailure;
    }

    if (isnan(bTagVal)) return defaultFailure;

    if (bTagVal > 1.) return 1.;
    if (bTagVal < 0.) return defaultFailure;

    return bTagVal;
}

bool CSVHelperSkim::PassesCSV(const pat::Jet &iJet, std::string taggername, const CSVHelperSkim::CSVwp iCSVworkingPoint, std::string dataEra)
{
    float csvValue = CSVHelperSkim::GetJetCSV(iJet, taggername);
    // CSV b-tagging requirement
    if (csvValue > CSVHelperSkim::GetWP(dataEra, iCSVworkingPoint, taggername)) { return true; }
    else {
        return false;
    }
}

float CSVHelperSkim::GetWP(std::string dataEra, const CSVHelperSkim::CSVwp iCSVworkingPoint, std::string taggername)
{
    if (TString(dataEra).Contains("2016")) {
        if (taggername == "DeepCSV") {
            switch (iCSVworkingPoint) {
                case CSVwp::Loose: {
                    return 0.2217;
                } break;
                case CSVwp::Medium: {
                    return 0.6321;
                } break;
                case CSVwp::Tight: {
                    return 0.8953;
                } break;
                case CSVwp::None: return 0;
            }
        }
        else if (taggername == "DeepJet") {
            switch (iCSVworkingPoint) {
                case CSVwp::Loose: {
                    return 0.0614;
                } break;
                case CSVwp::Medium: {
                    return 0.3093;
                } break;
                case CSVwp::Tight: {
                    return 0.7221;
                } break;
                case CSVwp::None: return 0;
            }
        }
        else if (taggername == "CSVv2") {
            switch (iCSVworkingPoint) {  // CSVv2 not supported for 2016 Legacy->WP are
                                         // assumed to be the same ones as for 2017 by me
                case CSVwp::Loose: {
                    return 0.5803;
                } break;
                case CSVwp::Medium: {
                    return 0.8838;
                } break;
                case CSVwp::Tight: {
                    return 0.9693;
                } break;
                case CSVwp::None: return 0;
            }
        }
        else {
            throw cms::Exception("CSVHelper: Invalid taggername ")
                << "Taggername '" << taggername << "' not recognized, only DeepCSV/DeepJet/CSVv2 possible" << std::endl;
            return 0;
        }
        return 0;
    }
    else if (TString(dataEra).Contains("2017")) {
        if (taggername == "DeepCSV") {
            switch (iCSVworkingPoint) {
                case CSVwp::Loose: {
                    return 0.1522;
                } break;
                case CSVwp::Medium: {
                    return 0.4941;
                } break;
                case CSVwp::Tight: {
                    return 0.8001;
                } break;
                case CSVwp::None: return 0;
            }
        }
        else if (taggername == "DeepJet") {
            switch (iCSVworkingPoint) {
                case CSVwp::Loose: {
                    return 0.0521;
                } break;
                case CSVwp::Medium: {
                    return 0.3033;
                } break;
                case CSVwp::Tight: {
                    return 0.7489;
                } break;
                case CSVwp::None: return 0;
            }
        }
        else if (taggername == "CSVv2") {  // CAREFUL: no WP avaiable !!!
            switch (iCSVworkingPoint) {
                case CSVwp::Loose: {
                    return 0.5803;
                } break;
                case CSVwp::Medium: {
                    return 0.8838;
                } break;
                case CSVwp::Tight: {
                    return 0.9693;
                } break;
                case CSVwp::None: return 0;
            }
        }
        else {
            throw cms::Exception("CSVHelper: Invalid taggername ")
                << "Taggername '" << taggername << "' not recognized, only DeepCSV/DeepJet/CSVv2 possible" << std::endl;
            return 0;
        }
    }
    else if (TString(dataEra).Contains("2018")) {
        if (taggername == "DeepCSV") {
            switch (iCSVworkingPoint) {
                case CSVwp::Loose: {
                    return 0.1241;
                } break;
                case CSVwp::Medium: {
                    return 0.4184;
                } break;
                case CSVwp::Tight: {
                    return 0.7527;
                } break;
                case CSVwp::None: return 0;
            }
        }
        else if (taggername == "DeepJet") {
            switch (iCSVworkingPoint) {
                case CSVwp::Loose: {
                    return 0.0494;
                } break;
                case CSVwp::Medium: {
                    return 0.2770;
                } break;
                case CSVwp::Tight: {
                    return 0.7264;
                } break;
                case CSVwp::None: return 0;
            }
        }
        else if (taggername == "CSVv2") {  // CAREFUL: no WP avaiable !!!
            switch (iCSVworkingPoint) {
                case CSVwp::Loose: {
                    return 0.5803;
                } break;
                case CSVwp::Medium: {
                    return 0.8838;
                } break;
                case CSVwp::Tight: {
                    return 0.9693;
                } break;
                case CSVwp::None: return 0;
            }
        }
        else {
            throw cms::Exception("CSVHelper: Invalid taggername ")
                << "Taggername '" << taggername << "' not recognized, only DeepCSV/DeepJet/CSVv2 possible" << std::endl;
            return 0;
        }
    }
    else {
        throw cms::Exception("CSVHelper: Invalid dataEra ") << "dataEra '" << dataEra << "' not recognized, only 2016/2017/2018 data possible" << std::endl;
        return 0;
    }
    return 0;
}
