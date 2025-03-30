
#ifndef MARTERO_UTIL_H
#define MARTERO_UTIL_H

#include "martero_selection.h"
// #include "selection.h"

// return (reco_E - slice.truth.E)/slice.truth.E [for sanity check]
const SpillVar count_slices_maria ([](const caf::SRSpillProxy *sr) -> double {
    bool debug = false;
    double recoE = -9999.9;
    double trueE = -9999.9;
    double slicetrueE = -9999.9;

   
    int count = 0;
    for (auto const& islc : sr->slc){
        // FIX ME check the classification_type
        if( islc.truth.index >= 0 && classification_type(sr,islc) == 2 ) {
            int ipfp_mu = -1;
            
            // if(automatic_selection_1muNp(sr, islc, 10, 100 )){
                ipfp_mu = find_muon(islc,10);
                // if(ipfp_mu != -1){
                    // recoE = Neutrino_energy_reco_Np(islc,ipfp_mu,10);
                    // slicetrueE = islc.truth.E;
                    if( debug )
                        std::cout << "[k1mu1p_RecoMinusSlicetruthOverSlicetruth] Sanity check: true energy (mc.nu) is " << trueE << "\n"
                                  << "True energy from truth-matching is " << slicetrueE << "\n"
                                  << "This is reco energy " << recoE
                                  << std::endl;
	    
                    count++;
                // } // found muon
            // } // automatic selection ok
	    } // truth info says this is a neutrino slice
    } // loop over slices
    
    return count; 
});

const SpillVar k1mu1p_RecoMinusSlicetruthOverSlicetruth([](const caf::SRSpillProxy *sr) -> double {
    bool debug = false;
    double recoE = -9999;
    double trueE = -9999;
    double slicetrueE = -9999;

    // if truth level is correct, find the match in reconstruction
    int count = 0;
    for (auto const& islc : sr->slc){
        // FIX ME check the classification_type
        if( islc.truth.index >= 0 && (classification_type(sr,islc) == 2 || classification_type(sr,islc) == 5) ) {
            int ipfp_mu = -1;
            
            if(automatic_selection_1muNp(sr, islc, 10, 100 )){
                ipfp_mu = find_muon(islc,10); 
                if(ipfp_mu != -1){
                    recoE = Neutrino_energy_reco_Np(islc,ipfp_mu,10);
                    slicetrueE = islc.truth.E;
	    
                    count++;
                } // found muon
            } // automatic selection ok
	    } // truth info says this is a neutrino slice
    } // loop over slices
    
    if( recoE < 0)
        return -9999;
    if( debug )
        std::cout<< "[k1mu1p_RecoMinusSlicetruthOverSlicetruth] This is recoE final " << recoE << endl;
    if( count > 1 )
        std::cout << "[k1mu1p_RecoMinusSlicetruthOverSlicetruth] Potential problems!!!!! " << endl;
    double ratio = -9999;
  
    ratio = (recoE - slicetrueE)/slicetrueE;
    if( debug )
        std::cout<< "[k1mu1p_RecoMinusSlicetruthOverSlicetruth] This is ratio " << ratio << std::endl;
    return ratio; 
});
#endif