#ifndef SELECTION_H
#define SELECTION_H

#include "helper_selection.h"

namespace selection
{
    ana::SpillMultiVar nu_energy_reco_1mu1p([](const caf::Proxy<caf::StandardRecord> *std_record) -> std::vector<double> {

        std::vector<double> slices_1muNp_Ereco_neutrino;

        for (auto const& slice: std_record->slc) {
            int ipfp_mu = -1;
            int ipfp_proton = -1;

            if (!automatic_selection_1muNp(std_record, slice, 10, 100)) // automatic selecting 1µNp0π slice
                continue;
            
            ipfp_mu = automatic_selection_mu_index(std_record, slice, 10, 100);
            double max_lenght = -1; // Looking at the longest proton in the slice...

            for (std::size_t ipfp=0; ipfp<slice.reco.npfp; ++ipfp) {
                if ((int)ipfp == ipfp_mu) 
                    continue;

                if (id_pfp(slice, ipfp, 10) == 1) {
                    if (slice.reco.pfp.at(ipfp).trk.len > max_lenght) {
                        max_lenght = slice.reco.pfp.at(ipfp).trk.len;
                        ipfp_proton = ipfp;
                    }
                } // selecting longhest proton...
            } 

            if (ipfp_mu > 0 && ipfp_proton > 0) 
                slices_1muNp_Ereco_neutrino.push_back(Transverse_mom_reco_Np(slice, ipfp_mu)); 
            
        } // looping over slices

        return slices_1muNp_Ereco_neutrino;
    });
} // namespace selection

#endif