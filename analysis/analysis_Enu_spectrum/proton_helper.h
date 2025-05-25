#ifndef PROTON_HELPER_H
#define PROTON_HELPER_H

#include "selection.h"

namespace vars {
    namespace truth {

        using namespace particle_data;
        using namespace var_utils;

        const ana::SpillVar spill_Np (
            double __def_ret, 
            const ana::Cut &reco_cut = ana::kNoCut,  
            cut_type_t what_to_cut_on = cut_type_t::RECO, 
            const ana::Cut &truth_cut = ana::kNoCut
        ) {
            return ana::SpillVar([=](const caf::SRSpillProxy *spill) -> double {
                int selected_slices = 0, num_protons_above50 = 0;
                double slice_value = __def_ret;
                bool debug = false;
        
                for (auto const& slice: spill->slc) {
        
                    if (what_to_cut_on == cut_type_t::RECO && !(reco_cut(&slice) && truth_cut(&slice))) {
                        if (debug) std::cout << "The mode was RECO but the cuts rejected this event" << std::endl;
                        continue;
                    } 
                    if (what_to_cut_on == cut_type_t::TRUE_1muN1p && !(
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1muNp && truth_cut(&slice)
                    )) {
                        if (debug) std::cout << "The mode was TRUE_1µ1p but the cuts rejected this event" << std::endl;
                        continue;
                    }
                    if (what_to_cut_on == cut_type_t::TRUE_1mu1p && !(
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1mu1p && truth_cut(&slice)
                    )) {
                        if (debug) std::cout << "The mode was TRUE_1µ1p but the cuts rejected this event" << std::endl;
                        continue;
                    }
                    if (what_to_cut_on == cut_type_t::TRUE_1mu2p && !(
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1mu2p && truth_cut(&slice)
                    )) {
                        if (debug) std::cout << "The mode was TRUE_1µ2p but the cuts rejected this event" << std::endl;
                        continue;
                    }
                    if (what_to_cut_on == cut_type_t::TRUE_1mu3p && !(
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1mu3p && truth_cut(&slice)
                    )) {
                        if (debug) std::cout << "The mode was TRUE_1µ3p but the cuts rejected this event" << std::endl;
                        continue;
                    }
                    if (what_to_cut_on == cut_type_t::TRUE_1muNp && !((
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1mu1p || 
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1muNp) && 
                        truth_cut(&slice)
                    )) {
                        if (debug) std::cout << "The mode was TRUE_1µN1p but the cuts rejected this event" << std::endl;
                        continue;
                    }
                    if (what_to_cut_on == cut_type_t::BOTH_1muNp && !(
                        reco_cut(&slice) && truth_cut(&slice) && 
                        (classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1muNp || 
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1mu1p)
                    )) {
                        if (debug) std::cout << "The mode was BOTH_1µNp (N>=1) but the cuts rejected this event" << std::endl;
                        continue;
                    }
                    if (what_to_cut_on == cut_type_t::BOTH_1muN1p && !(
                        reco_cut(&slice) && truth_cut(&slice) || 
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1muNp
                    )) {
                        if (debug) std::cout << "The mode was BOTH_1µN1p (N>1) but the cuts rejected this event" << std::endl;
                        continue;
                    }
                    if (what_to_cut_on == cut_type_t::BOTH_1mu1p && !(
                        reco_cut(&slice) && truth_cut(&slice) && 
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1mu1p
                    )) {
                        if (debug) std::cout << "The mode was BOTH_1mu1p but the cuts rejected this event" << std::endl;
                        continue;
                    }
                    if (what_to_cut_on == cut_type_t::BOTH_1mu2p && !(
                        reco_cut(&slice) && truth_cut(&slice) && 
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1mu2p
                    )) {
                        if (debug) std::cout << "The mode was BOTH_1mu2p but the cuts rejected this event" << std::endl;
                        continue;
                    }
                    if (what_to_cut_on == cut_type_t::BOTH_1mu3p && !(
                        reco_cut(&slice) && truth_cut(&slice) && 
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1mu3p
                    )) {
                        if (debug) std::cout << "The mode was BOTH_1mu3p but the cuts rejected this event" << std::endl;
                        continue;
                    }
                    if (what_to_cut_on == cut_type_t::BOTH_1muN3p && !(
                        reco_cut(&slice) && truth_cut(&slice) && 
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1muNp &&
                        classification_type(spill, &slice) != particle_data::int_type_t::true_visible_1mu3p &&
                        classification_type(spill, &slice) != particle_data::int_type_t::true_visible_1mu2p &&
                        classification_type(spill, &slice) != particle_data::int_type_t::true_visible_1mu1p            
                    )) {
                        if (debug) std::cout << "The mode was BOTH_1mu3p but the cuts rejected this event" << std::endl;
                        continue;
                    }
    
                    // /*
                    //  * true_visible_1mu1pi     | BOTH_1mu1pi
                    if (what_to_cut_on == cut_type_t::BOTH_1mu1pi && !(
                        reco_cut(&slice) && truth_cut(&slice) && 
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1mu1pi
                    )) {
                        if (debug) std::cout << "The mode was BOTH_1mu1pi but the cuts rejected this event" << std::endl;
                        continue;
                    }
                    //  * true_visible_1mu2pi     | BOTH_1mu2pi
                    if (what_to_cut_on == cut_type_t::BOTH_1mu2pi && !(
                        reco_cut(&slice) && truth_cut(&slice) && 
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1mu2pi
                    )) {
                        if (debug) std::cout << "The mode was BOTH_1mu2pi but the cuts rejected this event" << std::endl;
                        continue;
                    }
                    //  * true_visible_1mu3pi     | BOTH_1mu3pi
                    if (what_to_cut_on == cut_type_t::BOTH_1mu3pi && !(
                        reco_cut(&slice) && truth_cut(&slice) && 
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1mu3pi
                    )) {
                        if (debug) std::cout << "The mode was BOTH_1mu3pi but the cuts rejected this event" << std::endl;
                        continue;
                    }
                    //  * true_visible_1muNpi     | BOTH_1muNpi
                    if (what_to_cut_on == cut_type_t::BOTH_1muNpi && !(
                        reco_cut(&slice) && truth_cut(&slice) && 
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1muNpi
                    )) {
                        if (debug) std::cout << "The mode was BOTH_1muNpi but the cuts rejected this event" << std::endl;
                        continue;
                    }
                    //  * true_visible_2mu        | BOTH_2mu
                    if (what_to_cut_on == cut_type_t::BOTH_2mu && !(
                        reco_cut(&slice) && truth_cut(&slice) && 
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_2mu
                    )) {
                        if (debug) std::cout << "The mode was BOTH_2mu but the cuts rejected this event" << std::endl;
                        continue;
                    }
                    //  * true_visible_2p         | BOTH_2p
                    if (what_to_cut_on == cut_type_t::BOTH_2p && !(
                        reco_cut(&slice) && truth_cut(&slice) && 
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_2p
                    )) {
                        if (debug) std::cout << "The mode was BOTH_2p but the cuts rejected this event" << std::endl;
                        continue;
                    }
                    //  * true_visible_1muShortNp | BOTH_1muShortNp
                    if (what_to_cut_on == cut_type_t::BOTH_1muShortNp && !(
                        reco_cut(&slice) && truth_cut(&slice) && 
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1muShortNp
                    )) {
                        if (debug) std::cout << "The mode was BOTH_1muShortNp but the cuts rejected this event" << std::endl;
                        continue;
                    }
                    //  * true_visible_1pi1p      | BOTH_1pi1p
                    if (what_to_cut_on == cut_type_t::BOTH_1pi1p && !(
                        reco_cut(&slice) && truth_cut(&slice) && 
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1pi1p
                    )) {
                        if (debug) std::cout << "The mode was BOTH_1pi1p but the cuts rejected this event" << std::endl;
                        continue;
                    }                
                    //  * true_visible_1mu1p1pi0  | BOTH_1mu1p1pi0
                    if (what_to_cut_on == cut_type_t::BOTH_1mu1p1pi0 && !(
                        reco_cut(&slice) && truth_cut(&slice) && 
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1mu1p1pi0
                    )) {
                        if (debug) std::cout << "The mode was BOTH_1mu1p1pi0 but the cuts rejected this event" << std::endl;
                        continue;
                    }
                    //  * true_visible_1mu1pNpi0  | BOTH_1mu1pNpi0
                    if (what_to_cut_on == cut_type_t::BOTH_1mu1pNpi0 && !(
                        reco_cut(&slice) && truth_cut(&slice) && 
                        classification_type(spill, &slice) == particle_data::int_type_t::true_visible_1mu1pNpi0
                    )) {
                        if (debug) std::cout << "The mode was BOTH_1mu1pNpi0 but the cuts rejected this event" << std::endl;
                        continue;
                    }
                    // */
                    
                    if (!cuts::reco::slice_at_least_mu(&slice))
                        continue;
    
                    
                    for (auto const& prim: slice.truth.prim) {
                        double dep_E = 0; 
                        int use_plane = 2;
                        if (prim.G4ID < 0)
                            continue;
                    
                        if (prim.cryostat < 0)
                            continue;

                        if (std::abs(prim.pdg) == 2212) {
                            if (prim.daughters.size() > 0) {
                                for(auto const& true_particle: spill->true_particles) {
                                    if (true_particle.parent == prim.G4ID) {
                                        dep_E += true_particle.plane[prim.cryostat][use_plane].visE * particle_data::GeV;
                                    }
                                } // loop trough true_particles
                            } // proton has daughters
                            dep_E += prim.plane[prim.cryostat][use_plane].visE * particle_data::GeV;
                        } // found proton
            
                        if (std::abs(prim.pdg) == 2212 && dep_E > 50.0)
                            num_protons_above50 += 1;
                    }

                    selected_slices ++ ;
                } // loop spill->slc
            
                if (selected_slices > 1) 
                    logger::log(level_t::error) << "Something wrong with run:event = " 
                                                << run(spill) << ":" << event(spill) 
                                                << " => found " << selected_slices
                                                << " slice(s) 1µNp"
                                                << std::endl;
                return num_protons_above50;
            });
        }

    } // namespace truth
} // namespace vars

#endif