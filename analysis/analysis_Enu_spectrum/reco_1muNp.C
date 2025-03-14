
#ifndef recoEnu_efficiency_C
#define recoEnu_efficiency_C

#include "helper.h"
#include "selection.h"

template <class T> 
struct variables {
    std::string name;
    std::vector<T> vars;
};

struct tree {
    ana::Tree *tree;
    std::string name;
};


void reco_1muNp() {
    
    // REMARK: not all events are in both processing, some are only in the 
    // not cheated (17700 evt), nd not in the cheated sample (17500)...
    // These are accounted for in the hard_code namespace (hard_code::is_bad_event...)
    ana::SpectrumLoader loader_non_cheated("msotgia_v09_89_01_01p03_BNB_production_non_cheated_reco_ana_stage1tocaf_flatcafs");
    ana::SpectrumLoader loader_cheated("msotgia_v09_89_01_01p03_BNB_production_cheated_reco_ana_stage1tocaf_flatcafs");


    std::unique_ptr<ana::Tree> cheated(new ana::Tree(
        "cheated", 
        {"reco_E"}, 
        loader_cheated, 
        {vars::reco::slice_neutrino_energy_1muNp}, 
        cheating::cut_bad_events, 
        cuts::truth::slice_numuCC && cuts::reco::slice_at_least_mu
    )); 

    std::unique_ptr<ana::Tree> non_cheated(new ana::Tree(
        "non_cheated", 
        {"reco_E"}, 
        loader_non_cheated, 
        {vars::reco::slice_neutrino_energy_1muNp}, 
        cheating::cut_bad_events, 
        cuts::truth::slice_numuCC && cuts::reco::slice_at_least_mu
    )); 

    loader_cheated.Go();
    loader_non_cheated.Go();

    std::unique_ptr<TFile> file(new TFile("reco1muNp_alldata.root", "RECREATE"));
    cheated->SaveTo(file->mkdir("cheated"));
    non_cheated->SaveTo(file->mkdir("non_cheated"));
};

#endif