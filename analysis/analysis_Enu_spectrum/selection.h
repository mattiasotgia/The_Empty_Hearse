#ifndef SELECTION_H
#define SELECTION_H

#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TFile.h"
#include "TVector3.h"
#include "TProfile.h"

#include "helper.h"

#include <iostream>
#include <string>
#include <vector>

namespace particle_data {
    namespace masses {
        double pion = 139.57039;        // PDG value 2024 [MeV]
        double proton = 938.27208816;   // PDG value 2024 [MeV]
        double muon = 105.6583755;      // PDG value 2024 [MeV]
    }

    enum int_type {
        true_visible_1muNp,
        true_visible_1mu1p,
        unclassified
    };

    double minimum_gamma_MeV = 25.;
} // namespace particle_data

namespace var_utils {

    using level_t = logger::level;

    struct values_minmax {
        double min;
        double max;
    }; 
    
    struct chi2
    {
        double muon;
        double proton;
        double kaon;
        double pi;
    }; 

    double dist_cut = 10.;
    values_minmax barycenterFM_deltaZ_Trigger = {0., 100.};

    std::string dEdx_temp = 
        "/exp/icarus/app/users/msotgia/analysis/sbnana_v09_93_01_thesis_analysis/analysis/dEdxrestemplates.root";
    TFile* file = TFile::Open(dEdx_temp.c_str());

    auto dedx_range_pro = (TProfile *)file->Get("dedx_range_pro");
    auto dedx_range_ka  = (TProfile *)file->Get("dedx_range_ka");
    auto dedx_range_pi  = (TProfile *)file->Get("dedx_range_pi");
    auto dedx_range_mu  = (TProfile *)file->Get("dedx_range_mu");

    enum cut_type {
        RECO,
        TRUE_1muN1p,
        TRUE_1muNp,
        TRUE_1mu1p,
        BOTH_1muN1p,
        BOTH_1muNp,
        BOTH_1mu1p,
        MC_1muNp,
        MC_1mu1p
    };

    particle_data::int_type classification_type (const caf::SRSpillProxy*, const caf::SRSliceProxy*);
    const ana::SpillVar make_spill_from_slice (
        const ana::Var &slice_var, 
        const ana::Cut &reco_cut = ana::kNoCut,  
        cut_type what_to_cut_on = cut_type::RECO, 
        const ana::Cut &truth_cut = ana::kNoCut,
        bool info = false, double treshold = 0.25,
        std::function<particle_data::int_type(
            const caf::SRSpillProxy*, const caf::SRSliceProxy*
        )> classification = classification_type
    ) {
        return ana::SpillVar([=](const caf::SRSpillProxy *spill) -> double {
            int selected_slices = 0;
            double slice_value = -9999;
    
            // if (what_to_cut_on == cut_type::MC_1muNp || what_to_cut_on == cut_type::MC_1mu1p) {
            //     for (auto const& nu: spill->mc.nu) {
            //         if (what_to_cut_on == cut_type::MC_1muNp && 
            //             classification_type_MC(spill, nu) != particle_data::int_type::true_visible_1muNp) continue;
            //         if (what_to_cut_on == cut_type::MC_1mu1p && 
            //             classification_type_MC(spill, nu) != particle_data::int_type::true_visible_1mu1p) continue;
                    
                    
            //     } // loop spill->nu
            // } // MC classification :)
            for (auto const& slice: spill->slc) {
    
                if (what_to_cut_on == cut_type::RECO && !reco_cut(&slice)) continue; 
                if (what_to_cut_on == cut_type::TRUE_1muNp && (
                    classification(spill, &slice) != particle_data::int_type::true_visible_1muNp || !truth_cut(&slice)
                ))
                    continue;
                if (what_to_cut_on == cut_type::TRUE_1mu1p && (
                    classification(spill, &slice) != particle_data::int_type::true_visible_1mu1p || !truth_cut(&slice)
                ))
                    continue;
                if (what_to_cut_on == cut_type::BOTH_1muNp && (
                    !reco_cut(&slice) || !truth_cut(&slice) || 
                    classification(spill, &slice) != particle_data::int_type::true_visible_1muNp || classification(spill, &slice) != particle_data::int_type::true_visible_1mu1p
                ))
                    continue;
                if (what_to_cut_on == cut_type::BOTH_1muN1p && (
                    !reco_cut(&slice) || !truth_cut(&slice) || 
                    classification(spill, &slice) != particle_data::int_type::true_visible_1muNp
                ))
                    continue;
                if (what_to_cut_on == cut_type::BOTH_1mu1p && (
                    !reco_cut(&slice) || !truth_cut(&slice) || 
                    classification(spill, &slice) != particle_data::int_type::true_visible_1mu1p
                ))
                    continue;
                
                slice_value = slice_var(&slice);
                if (slice_value > treshold && info) {
                    std::cout << "Slice value = " << slice_value << " > " << treshold << " for run.event = " << run(spill) << "." << event(spill) << std::endl;
                }
                selected_slices ++ ;
            } // loop spill->slc
        
            if (selected_slices > 1) 
                logger::log(level_t::error) << "Something wrong with run:event = " 
                                            << run(spill) << ":" << event(spill) 
                                            << " => found " << selected_slices
                                            << " slice(s) 1µNp"
                                            << std::endl;
            return slice_value;
        });
    }

    chi2 chi2_ALG (std::vector<double> &dEdx, std::vector<double> &RR, double rr_min, double rr_max) {
        /* The output is chi2s
         * Those are used to classify final state particles 
        */

        double threshold = 0.5;
        double max_rr = rr_max; // Max value for the residual range (RR)
        double min_rr = rr_min; // Min value for the residual range (RR)

        std::vector<float> trkdedx;
        std::vector<float> trkres;
        std::vector<double> vpida;

        for (std::size_t i(0); i < dEdx.size(); ++i) {
            if (i == 0 || i == dEdx.size() - 1)
                continue;
            if (RR[i] < max_rr && RR[i] > rr_min) {
                trkdedx.push_back(dEdx[i]);
                trkres.push_back(RR[i]);
            }
        }

        int npt = 0;
        double chi2pro = 0;
        double chi2ka = 0;
        double chi2pi = 0;
        double chi2mu = 0;
        double avgdedx = 0;
        double PIDA = 0;

        int used_trkres = 0;
        for (unsigned i = 0; i < trkdedx.size(); ++i) { 
            // hits
            // ignore the first and the last point
            // if (i == 0 || i == trkdedx.size() - 1) continue;
            avgdedx += trkdedx[i];
            if (trkres[i] < 26) {
                PIDA += trkdedx[i] * std::pow(trkres[i], 0.42);
                vpida.push_back(trkdedx[i] * std::pow(trkres[i], 0.42));
                used_trkres++;
            }
            if (trkdedx[i] > 100 || trkdedx[i] < threshold)
                continue; // protect against large pulse height
            
            int bin = dedx_range_pro->FindBin(trkres[i]);
            if (bin >= 1 && bin <= dedx_range_pro->GetNbinsX()) {
                
                double bincpro = dedx_range_pro->GetBinContent(bin);
                if (bincpro < 1e-6)  // for 0 bin content, using neighboring bins
                    bincpro = 
                        (dedx_range_pro->GetBinContent(bin - 1) + dedx_range_pro->GetBinContent(bin + 1)) / 2;
                
                double bincka = dedx_range_ka->GetBinContent(bin);
                if (bincka < 1e-6)
                    bincka =
                        (dedx_range_ka->GetBinContent(bin - 1) + dedx_range_ka->GetBinContent(bin + 1)) / 2;

                double bincpi = dedx_range_pi->GetBinContent(bin);
                if (bincpi < 1e-6)
                    bincpi =
                        (dedx_range_pi->GetBinContent(bin - 1) + dedx_range_pi->GetBinContent(bin + 1)) / 2;
                
                double bincmu = dedx_range_mu->GetBinContent(bin);
                if (bincmu < 1e-6)
                    bincmu =
                        (dedx_range_mu->GetBinContent(bin - 1) + dedx_range_mu->GetBinContent(bin + 1)) / 2;
                
                double binepro = dedx_range_pro->GetBinError(bin);
                if (binepro < 1e-6)
                    binepro =
                        (dedx_range_pro->GetBinError(bin - 1) + dedx_range_pro->GetBinError(bin + 1)) / 2;

                double bineka = dedx_range_ka->GetBinError(bin);
                if (bineka < 1e-6)
                    bineka = (dedx_range_ka->GetBinError(bin - 1) + dedx_range_ka->GetBinError(bin + 1)) / 2;

                double binepi = dedx_range_pi->GetBinError(bin);
                if (binepi < 1e-6)
                    binepi = (dedx_range_pi->GetBinError(bin - 1) + dedx_range_pi->GetBinError(bin + 1)) / 2;

                double binemu = dedx_range_mu->GetBinError(bin);
                if (binemu < 1e-6)
                    binemu = (dedx_range_mu->GetBinError(bin - 1) + dedx_range_mu->GetBinError(bin + 1)) / 2;

                // double errke = 0.05*trkdedx[i];   //5% KE resolution

                double errdedx = 0.04231 + 0.0001783 * trkdedx[i] * trkdedx[i]; // resolution on dE/dx
                
                errdedx *= trkdedx[i];
                
                chi2pro += std::pow((trkdedx[i] - bincpro) / std::sqrt(std::pow(binepro, 2) + std::pow(errdedx, 2)), 2);
                chi2ka  += std::pow((trkdedx[i] - bincka)  / std::sqrt(std::pow(bineka, 2)  + std::pow(errdedx, 2)), 2);
                chi2pi  += std::pow((trkdedx[i] - bincpi)  / std::sqrt(std::pow(binepi, 2)  + std::pow(errdedx, 2)), 2);
                chi2mu  += std::pow((trkdedx[i] - bincmu)  / std::sqrt(std::pow(binemu, 2)  + std::pow(errdedx, 2)), 2);
                // std::cout<<i<<" "<<trkdedx[i]<<" "<<trkres[i]<<" "<<bincpro<<std::endl;
                ++npt;
            }
        } // hits

        return {chi2mu / npt, chi2pro / npt, chi2ka / npt, chi2pi / npt};
    } // chi2 chi2_ALG

    int find_muon (const caf::SRSliceProxy &slice, int dist_mucut) {

        /* Find the pfp number of the muon 
        */

        // Select muon as longest track
        double max_length = -1.0;
        int ipfp_mu = -1;
        TVector3 reco_vtx;
        TVector3 reco_start;
        
        reco_vtx.SetXYZ(slice.vertex.x, slice.vertex.y, slice.vertex.z);
        
        for (std::size_t ipfp = 0; ipfp < slice.reco.npfp; ++ipfp) {
            if (std::isnan(slice.reco.pfp[ipfp].trk.start.x) || std::isnan(slice.reco.pfp[ipfp].trk.len))
                continue;
        
            // if(!isInDetector(slice.reco.pfp[ipfp].trk.start.x,slice.reco.pfp[ipfp].trk.start.y,slice.reco.pfp[ipfp].trk.start.z))   continue;
            // if(!isInDetector(slice.reco.pfp[ipfp].trk.end.x,slice.reco.pfp[ipfp].trk.end.y,slice.reco.pfp[ipfp].trk.end.z))         continue;

            reco_start.SetXYZ(slice.reco.pfp[ipfp].trk.start.x, slice.reco.pfp[ipfp].trk.start.y, slice.reco.pfp[ipfp].trk.start.z);
            
            /* For µ the selection requires a trackScore > 0.5 */
            if (slice.reco.pfp[ipfp].trackScore < 0.5)
                continue;

            // int use_plane = slice.reco.pfp[ipfp].trk.calo[2].nhit>slice.reco.pfp[ipfp].trk.calo[1].nhit ? 2:1;
            int use_plane = 2; // Always using collection?
            
            // compute new chi2
            std::vector<double> dedx;
            std::vector<double> rr;
            for (std::size_t ihit = 0; ihit < slice.reco.pfp[ipfp].trk.calo[use_plane].points.size(); ++ihit) {
                dedx.push_back(slice.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].dedx);
                rr.push_back(slice.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].rr);
            } // calo points

            // input to chi2_ALG vector dedx, vector rr, rr_min, rr_max
            // output chi2s {chi2mu/npt,chi2pro/npt,chi2ka/npt,chi2pi/npt}
            chi2 chi2_values = chi2_ALG(dedx, rr, 0.0, 25.0);

            if (
                slice.reco.pfp[ipfp].trk.len > max_length   && 
                // slice.reco.pfp[ipfp].trk.len > 50           && // The selection of the longest muon being 50+ cm should not be here? // separate Cut below
                (reco_vtx - reco_start).Mag() < dist_mucut  && 
                chi2_values.muon < 30                       && // chi2 cuts on proton and muon 
                chi2_values.proton > 60                     && // chi2 cuts on proton and muon 
                // in_contained (slice.reco.pfp[ipfp].trk.end.x, slice.reco.pfp[ipfp].trk.end.y, slice.reco.pfp[ipfp].trk.end.z, 5.) && // separate Cut below
                // slice.reco.pfp[ipfp].trk.end.x * slice.vertex.x > 0 &&                                                               // separate Cut below
                slice.reco.pfp[ipfp].parent_is_primary
            ) {
                max_length = slice.reco.pfp[ipfp].trk.len;
                ipfp_mu = ipfp;
            }
        } // loop of pfp to find muon
        return ipfp_mu;
    } // int find_muon

    int id_pfp (const caf::SRSliceProxy &slice, int ipfp, int dist_cut) {
        /* This utility select the correct particle of the event slice
         * return 1 PROTONS
         * return 2 PIONS
         * return 3 SHOWER
         * return 9 other -> nan, not primary, too far, below energy threshold...
        */

        TVector3 rec_vtx;
        rec_vtx.SetXYZ(slice.vertex.x, slice.vertex.y, slice.vertex.z);

        if (!(slice.reco.pfp[ipfp].parent_is_primary))
            return 9;

        // Sanity checks :)
        if (
            std::isnan(slice.reco.pfp[ipfp].trk.start.x)    || 
            std::isnan(slice.reco.pfp[ipfp].trk.end.x)      || 
            std::isnan(slice.reco.pfp[ipfp].trk.len)
        ) return 9;
        
        // There is always a muon, for a 1mu1p we need 2 tracks - 1 muon = 1 only proton with threshold
        // if (int(ipfp)==ipfp_mu) continue;
        /* Consider only primary tracks which are 20cm 
         * close to the vertex, either vtx-start or vtx-end
        */

        TVector3 start (slice.reco.pfp[ipfp].trk.start.x, slice.reco.pfp[ipfp].trk.start.y, slice.reco.pfp[ipfp].trk.start.z);
        TVector3 end   (slice.reco.pfp[ipfp].trk.end.x,   slice.reco.pfp[ipfp].trk.end.y,   slice.reco.pfp[ipfp].trk.end.z);
        
        /* Some reminder of ternary operator
         * condition ? result_if_true : result_if_false
        */

        double min_dist = ((start - rec_vtx).Mag() < (end - rec_vtx).Mag()) ? (start - rec_vtx).Mag() : (end - rec_vtx).Mag();
        
        // if (min_dist > 50.0) continue;
        // if (min_dist > 10.0) return 9;
        if (min_dist > 50.0)
            return 9;

        // int use_plane = slice.reco.pfp[ipfp].trk.calo[2].nhit>slice.reco.pfp[ipfp].trk.calo[1].nhit ? 2:1;
        int use_plane = 2;
        // compute new chi2
        std::vector<double> dedx;
        std::vector<double> rr;
        for (std::size_t ihit = 0; ihit < slice.reco.pfp[ipfp].trk.calo[use_plane].points.size(); ++ihit) {
            dedx.push_back(slice.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].dedx);
            rr.push_back(slice.reco.pfp[ipfp].trk.calo[use_plane].points[ihit].rr);
        } // calo points

        // input to chi2_ALG vector dedx, vector rr, rr_min, rr_max
        // output chi2s {chi2mu/npt,chi2pro/npt,chi2ka/npt,chi2pi/npt}
        chi2 chi2_values = chi2_ALG(dedx, rr, 0.0, 25.0);

        if (slice.reco.pfp[ipfp].trackScore >= 0.5) {
            if (
                std::isnan(slice.reco.pfp[ipfp].trk.start.x)    || 
                std::isnan(slice.reco.pfp[ipfp].trk.end.x)      || 
                std::isnan(slice.reco.pfp[ipfp].trk.len)
            ) return 9;

            if (
                std::isnan(slice.reco.pfp[ipfp].trk.start.y) || 
                std::isnan(slice.reco.pfp[ipfp].trk.start.z) || 
                std::isnan(slice.reco.pfp[ipfp].trk.end.y)   || 
                std::isnan(slice.reco.pfp[ipfp].trk.end.z)
            ) return 9;

            // Skip low energy tagged pions
            TVector3 start_mom_V3;
        
            if (chi2_values.proton >= 100) 
                start_mom_V3.SetXYZ(
                    slice.reco.pfp[ipfp].trk.rangeP.p_pion * slice.reco.pfp[ipfp].trk.dir.x, 
                    slice.reco.pfp[ipfp].trk.rangeP.p_pion * slice.reco.pfp[ipfp].trk.dir.y, 
                    slice.reco.pfp[ipfp].trk.rangeP.p_pion * slice.reco.pfp[ipfp].trk.dir.z
                );
            
            if (
                chi2_values.proton >= 100           && 
                (rec_vtx - start).Mag() < dist_cut  && 
                std::sqrt ( std::pow(particle_data::masses::pion, 2) + std::pow(start_mom_V3.Mag() * 1000, 2) ) - particle_data::masses::pion >= 25.0 && 
                slice.reco.pfp[ipfp].parent_is_primary
            )  return 2;

            // Skip low energy protons
            if (chi2_values.proton < 100)
                start_mom_V3.SetXYZ(
                    slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.x, 
                    slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.y, 
                    slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.z
                );

            if (
                chi2_values.proton < 100            && 
                (rec_vtx - start).Mag() < dist_cut  && 
                std::sqrt ( std::pow(particle_data::masses::proton, 2) + std::pow(start_mom_V3.Mag() * 1000, 2) ) - particle_data::masses::proton >= 50.0 && 
                slice.reco.pfp[ipfp].parent_is_primary
            ) return 1;

        }

        if (slice.reco.pfp[ipfp].trackScore < 0.5) {
            if (slice.reco.pfp[ipfp].trackScore >= 0.4 && chi2_values.proton < 100) {
                TVector3 start_mom_V3_2;
                start_mom_V3_2.SetXYZ(
                    slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.x, 
                    slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.y, 
                    slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.z
                );

                if (
                    std::sqrt( std::pow(particle_data::masses::proton, 2) + std::pow(start_mom_V3_2.Mag() * 1000, 2)) - particle_data::masses::proton >= 50.0 && 
                    (rec_vtx - start).Mag() < dist_cut && 
                    slice.reco.pfp[ipfp].parent_is_primary
                ) return 1;
            }

            if (!(slice.reco.pfp[ipfp].trackScore >= 0.4 && chi2_values.proton < 100)) {

                // int use_plane2 = slice.reco.pfp[ipfp].trk.calo[2].nhit>slice.reco.pfp[ipfp].trk.calo[1].nhit ? 2:1;
                int use_plane2 = 2;

                if (std::isnan(slice.reco.pfp[ipfp].shw.plane[use_plane2].energy))
                    return 9;
                if (slice.reco.pfp[ipfp].shw.plane[use_plane2].energy * 1000 < 25.0)
                    return 9;
                if (slice.reco.pfp[ipfp].shw.plane[use_plane2].energy * 1000 > 25.0 && slice.reco.pfp[ipfp].parent_is_primary)
                    return 3;
            }
        }
        return 9;
    } // int id_pfp
} // namespace var_utils

namespace cuts {

    bool in_FV (double x, double y, double z) {
        if (std::isnan(x) || std::isnan(y) || std::isnan(z))
            return false;
        // need to add a check to avoid having the vtx slice in one cryo and the end track in the other cryo

        // fiducial volume for the dangling cable
        if (x > 210.0 && y > 60.0 && z > 290.0 && z < 390.0)
            return false;
        // Cathode padding
        // if(std::abs(x)>208.79 && std::abs(x)<211.79) return false;

        return (((x < -61.94 - 25 && x > -358.49 + 25) || // cheack on x direction cryo E/0
                 (x > 61.94 + 25 && x < 358.49 - 25))  && // cheack on x direction cryo W/1
                ((y > -181.86 + 25 && y < 134.96 - 25) && // cheack on y direction
                 (z > -894.95 + 30 && z < 894.95 - 50))); // cheack on z direction
    } // bool in_FV

    bool in_contained (double x, double y, double z, double dist = 5.) {
        if (std::isnan(x) || std::isnan(y) || std::isnan(z))
            return false;

        // 5 cm containment for Y in both cryo
        return (((x < -61.94 - dist && x > -358.49 + dist) || // cheack on x direction cryo E/0
                 (x > 61.94 + dist && x < 358.49 - dist))  && // cheack on x direction cryo W/1
                ((y > -181.86 + dist && y < 134.96 - dist) && // cheack on y direction
                (z > -894.95 + dist && z < 894.95 - dist)));  // cheack on z direction
    } // bool is_contained  

    bool in_active (double x, double y, double z) {
        if (std::isnan(x) || std::isnan(y) || std::isnan(z))
            return false;

        return (((x < -61.94 && x > -358.49) || // cheack on x direction cryo E/0
                 (x > 61.94 && x < 358.49))  && // cheack on x direction cryo W/1
                ((y > -181.86 && y < 134.96) && // cheack on y direction
                 (z > -894.95 && z < 894.95))); // cheack on z direction
    } // bool in_active


    bool in_detector (double x, double y, double z) {
        if (std::isnan(x) || std::isnan(y) || std::isnan(z))
            return false;
        return (((x < -61.94 + 5 && x > -358.49 - 5) || // cheack on x direction cryo E/0
                 (x > 61.94 - 5 && x < 358.49 + 5))  && // cheack on x direction cryo W/1
                ((y > -181.86 - 5 && y < 134.96 + 5) && // cheack on y direction
                 (z > -894.95 - 5 && z < 894.95 + 5))); // cheack on z direction
    } // bool in_detector
    
    bool all_trk_contained_truth (const caf::SRSpillProxy *spill, const caf::SRSliceProxy *slice) {

        for (auto const& prim: slice->truth.prim) {
            if (prim.G4ID < 0)     continue;
            if (prim.cryostat < 0) continue;
            if (
                std::abs(prim.pdg) != 13   &&   // mu
                std::abs(prim.pdg) != 2212 &&   // p
                std::abs(prim.pdg) != 211  &&   // pi
                std::abs(prim.pdg) != 11        // e 
            ) continue; // why are we selecting only charged primaries?
                    
            if (!in_contained(prim.end.x, prim.end.y, prim.end.z)) return false;

            if (prim.daughters.size() < 0) continue;
            int parent_g4id;

            for (auto const& true_particle: spill->true_particles) {
                parent_g4id = true_particle.parent;

                if (parent_g4id != prim.G4ID) continue;
                if (
                    std::abs(true_particle.pdg) != 13   &&  // mu
                    std::abs(true_particle.pdg) != 2212 &&  // p
                    std::abs(true_particle.pdg) != 211  &&  // pi
                    std::abs(true_particle.pdg) != 11       // e 
                ) continue; // why are we selecting only charged primaries?
                    
                /* L’idea e’ mettere soglia su quello che deposita energia…. 
                 * Quindi tutto quello che sia carico e che di solito abbiamo
                 * Per dire in qualche punto ci deve essere anche una soglia 
                 * sui fotoni singoli, per cui se hanno troppa energia depositata 
                 * possono essere identificati (non avevamo messo soglia sui pi0 
                 * pero si sui due fotoni del pi0)
                */

                if (!in_contained(true_particle.end.x, true_particle.end.y, true_particle.end.z)) return false;
            } // loop spill->true_particles
        } // loop slc.truth.prim
        return true;
    } // bool all_trk_contained_truth

    bool all_trk_contained_MC (const caf::SRSpillProxy *spill, const caf::Proxy<caf::SRTrueInteraction> &nu) {

        for (auto const& prim: nu.prim) {
            if (prim.G4ID < 0)     continue;
            if (prim.cryostat < 0) continue;
            if (
                std::abs(prim.pdg) != 13   &&   // mu
                std::abs(prim.pdg) != 2212 &&   // p
                std::abs(prim.pdg) != 211  &&   // pi
                std::abs(prim.pdg) != 11        // e 
            ) continue; // why are we selecting only charged primaries?
                    
            if (!in_contained(prim.end.x, prim.end.y, prim.end.z)) return false;

            if (prim.daughters.size() < 0) continue;
            int parent_g4id;

            for (auto const& true_particle: spill->true_particles) {
                parent_g4id = true_particle.parent;

                if (parent_g4id != prim.G4ID) continue;
                if (
                    std::abs(true_particle.pdg) != 13   &&  // mu
                    std::abs(true_particle.pdg) != 2212 &&  // p
                    std::abs(true_particle.pdg) != 211  &&  // pi
                    std::abs(true_particle.pdg) != 11       // e 
                ) continue; // why are we selecting only charged primaries?
                    
                /* L’idea e’ mettere soglia su quello che deposita energia…. 
                 * Quindi tutto quello che sia carico e che di solito abbiamo
                 * Per dire in qualche punto ci deve essere anche una soglia 
                 * sui fotoni singoli, per cui se hanno troppa energia depositata 
                 * possono essere identificati (non avevamo messo soglia sui pi0 
                 * pero si sui due fotoni del pi0)
                */

                if (!in_contained(true_particle.end.x, true_particle.end.y, true_particle.end.z)) return false;
            } // loop spill->true_particles
        } // loop slc.truth.prim
        return true;
    } // bool all_trk_contained_MC

    namespace truth {
        const ana::Cut slice_numuCC(ana::kIsNumuCC);

        const ana::Cut slice_vtx_in_FV ([](const caf::SRSliceProxy *slice) -> bool {
            return in_FV (slice->truth.position.x, slice->truth.position.y, slice->truth.position.z);
        }); // const ana::Cut slice_vtx_in_FV

        const ana::Cut slice_vtx_in_contained ([](const caf::SRSliceProxy *slice) -> bool {
            return in_contained (slice->truth.position.x, slice->truth.position.y, slice->truth.position.z, 5.);
        }); // const ana::Cut slice_vtx_in_contained

        const ana::Cut slice_vtx_in_active ([](const caf::SRSliceProxy *slice) -> bool {
            return in_active (slice->truth.position.x, slice->truth.position.y, slice->truth.position.z);
        }); // const ana::Cut slice_vtx_in_active

        const ana::SpillCut spill_has_numuCC ([](const caf::SRSpillProxy *spill) -> bool {
            for (auto const& nu: spill->mc.nu) {
                if (nu.pdg == std::numeric_limits<int>::min()) continue;
                if (nu.iscc && nu.pdg == 14) return true;
            }
            return false;
        });

        const ana::SpillCut spill_all_trk_contained_MC ([](const caf::SRSpillProxy *spill) -> bool {
            // Cut the spill if the spill->mc.nu event (iscc && pdg == 14)) has anything 
            // or daughters going out of the contained volume
            
            for (auto const& nu: spill->mc.nu) {

                if (!nu.iscc || nu.pdg != 14) continue;

                for (auto const& prim: nu.prim) {
                    if (prim.G4ID < 0)     continue;
                    if (prim.cryostat < 0) continue;
                    if (
                        std::abs(prim.pdg) != 13   && 
                        std::abs(prim.pdg) != 2212 && 
                        std::abs(prim.pdg) != 211  && 
                        std::abs(prim.pdg) != 11
                    ) continue; // why are we selecting only charged primaries?
                        
                    if (!in_contained(prim.end.x, prim.end.y, prim.end.z)) return false;

                    if (prim.daughters.size() < 0) continue;
                    int parent_g4id;

                    for (auto const& true_particle: spill->true_particles) {
                        parent_g4id = true_particle.parent;

                        if (parent_g4id != prim.G4ID) continue;
                        if (
                            std::abs(true_particle.pdg) != 13   && 
                            std::abs(true_particle.pdg) != 2212 && 
                            std::abs(true_particle.pdg) != 211  && 
                            std::abs(true_particle.pdg) != 11
                        ) continue; // why are we selecting only charged primaries?

                        /* L’idea e’ mettere soglia su quello che deposita energia…. 
                         * Quindi tutto quello che sia carico e che di solito abbiamo
                         * Per dire in qualche punto ci deve essere anche una soglia 
                         * sui fotoni singoli, per cui se hanno troppa energia depositata 
                         * possono essere identificati (non avevamo messo soglia sui pi0 
                         * pero si sui due fotoni del pi0)
                        */

                        if (!in_contained(true_particle.end.x, true_particle.end.y, true_particle.end.z)) return false;
                    } // loop spill->true_particles
                } // loop nu.prim
            } // loop spill->mc.nu
            return true;
        }); // const ana::SpillCut spill_all_trk_contained_MC

        const ana::SpillCut spill_all_trk_contained_truth ([](const caf::SRSpillProxy *spill) -> bool {
            // Cut the spill if the spill->slc slice (truth.iscc && truth.pdg == 14)) has anything 
            // or daughters going out of the contained volume

            for (auto const& slc: spill->slc) {
                if (!slc.truth.iscc || slc.truth.pdg != 14) continue;

                for (auto const& prim: slc.truth.prim) {
                    if (prim.G4ID < 0)     continue;
                    if (prim.cryostat < 0) continue;
                    if (
                        std::abs(prim.pdg) != 13   &&   // mu
                        std::abs(prim.pdg) != 2212 &&   // p
                        std::abs(prim.pdg) != 211  &&   // pi
                        std::abs(prim.pdg) != 11        // e 
                    ) continue; // why are we selecting only charged primaries?
                        
                    if (!in_contained(prim.end.x, prim.end.y, prim.end.z)) return false;

                    if (prim.daughters.size() < 0) continue;
                    int parent_g4id;

                    for (auto const& true_particle: spill->true_particles) {
                        parent_g4id = true_particle.parent;

                        if (parent_g4id != prim.G4ID) continue;
                        if (
                            std::abs(true_particle.pdg) != 13   &&  // mu
                            std::abs(true_particle.pdg) != 2212 &&  // p
                            std::abs(true_particle.pdg) != 211  &&  // pi
                            std::abs(true_particle.pdg) != 11       // e 
                        ) continue; // why are we selecting only charged primaries?
                        
                        /* L’idea e’ mettere soglia su quello che deposita energia…. 
                         * Quindi tutto quello che sia carico e che di solito abbiamo
                         * Per dire in qualche punto ci deve essere anche una soglia 
                         * sui fotoni singoli, per cui se hanno troppa energia depositata 
                         * possono essere identificati (non avevamo messo soglia sui pi0 
                         * pero si sui due fotoni del pi0)
                        */

                        if (!in_contained(true_particle.end.x, true_particle.end.y, true_particle.end.z)) return false;
                    } // loop spill->true_particles
                } // loop slc.truth.prim
            } // loop spill->slc
            return true;
        }); // const ana::SpillCut spill_all_trk_contained_truth

        const ana::Cut slice_1mu_only ([](const caf::SRSliceProxy *slice) -> bool {
            int muon_n = 0;
            for (auto const& prim: slice->truth.prim) {
                if (std::abs(prim.pdg) == 13) muon_n ++;
            } // loop slice->truth.prim
            return muon_n == 1;
        }); // const ana::Cut slice_1mu_only
    } // namespace truth

    namespace reco { 
        const ana::Cut slice_at_least_mu ([](const caf::SRSliceProxy *slice) -> bool {
            /* We should at least found one muon, otherwise something bad happens
            */

            // Using dist_cut = 10 [cm]
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut); 
            if (ipfp_muon == -1) return false; // no muon found ahaha
            return true;
        }); // const ana::Cut slice_at_least_mu

        const ana::Cut slice_all_trk_contained ([](const caf::SRSliceProxy *slice) -> bool {
            for (auto const& pfp: slice->reco.pfp) {
                if (std::isnan(pfp.trk.start.x) || std::isnan(pfp.trk.end.x) || std::isnan(pfp.trk.len)) 
                    continue;

                // Check meaningful points
                // if (!isInDetector(slice.reco.pfp[ipfp].trk.start.x,slice.reco.pfp[ipfp].trk.start.y,slice.reco.pfp[ipfp].trk.start.z)) continue;
                // if (!isInDetector(slice.reco.pfp[ipfp].trk.end.x,slice.reco.pfp[ipfp].trk.end.y,slice.reco.pfp[ipfp].trk.end.z)) continue;
                // if (!(slice.reco.pfp[ipfp].parent_is_primary )) continue; // Skip secondaries
                // if (slice.reco.pfp[ipfp].trackScore<0.4) continue;        // Want to check only tracks??

                if ((pfp.trk.start.x * slice->vertex.x) < 0) return false; // PFP crossing cryostats :(
                if (!in_contained(pfp.trk.end.x, pfp.trk.end.y, pfp.trk.end.z, 5.)) 
                    return false;
            } // loop pfp
            return true;
        }); // const ana::Cut slice_all_contained
        
        const ana::Cut slice_vtx_in_FV ([](const caf::SRSliceProxy *slice) -> bool {
            return in_FV (slice->vertex.x, slice->vertex.y, slice->vertex.z);
        }); // const ana::Cut slice_vtx_in_FV

        const ana::Cut slice_vtx_in_contained ([](const caf::SRSliceProxy *slice) -> bool {
            return in_contained (slice->vertex.x, slice->vertex.y, slice->vertex.z, 5.);
        }); // const ana::Cut slice_vtx_in_contained

        const ana::Cut slice_vtx_in_active ([](const caf::SRSliceProxy *slice) -> bool {
            return in_active (slice->vertex.x, slice->vertex.y, slice->vertex.z);
        }); // const ana::Cut slice_vtx_in_active

        const ana::Cut slice_mu_50_length ([](const caf::SRSliceProxy *slice) -> bool {
            /* The selected muon must also be longer than 50 cm
            */

            // Using dist_cut = 10 [cm]
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut); 
            if (ipfp_muon == -1) return false; // no muon found ahaha

            return slice->reco.pfp[ipfp_muon].trk.len > 50;
        }); // const ana::Cut slice_mu_50_length

        const ana::Cut slice_mu_in_contained ([](const caf::SRSliceProxy *slice) -> bool {
            /* The selected muon must also be in contained volume
            */

            // Using dist_cut = 10 [cm]
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut); 
            if (ipfp_muon == -1) return false; // no muon found ahaha

            return in_contained (slice->reco.pfp[ipfp_muon].trk.end.x, slice->reco.pfp[ipfp_muon].trk.end.y, slice->reco.pfp[ipfp_muon].trk.end.z, 5.);
        }); // const ana::Cut slice_mu_in_contained

        const ana::Cut slice_mu_not_crossing ([](const caf::SRSliceProxy *slice) -> bool {
            /* The selected muon must also be not crossing ;)
            */

            // Using dist_cut = 10 [cm]
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut); 
            if (ipfp_muon == -1) return false; // no muon found ahaha

            return slice->reco.pfp[ipfp_muon].trk.end.x * slice->vertex.x > 0;
        }); // const ana::Cut slice_mu_not_crossing

        const ana::Cut slice_barycenter ([](const caf::SRSliceProxy *slice) -> bool {
            return (
                slice->barycenterFM.deltaZ_Trigger < var_utils::barycenterFM_deltaZ_Trigger.max && 
                slice->barycenterFM.deltaZ_Trigger > var_utils::barycenterFM_deltaZ_Trigger.min
            );
        }); // const ana::Cut slice_barycenter

        const ana::Cut slice_1muNp ([](const caf::SRSliceProxy *slice) -> bool {
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            if (ipfp_muon == -1) return false; // redundant, btw, but who cares...

            int num_protons = 0;
            int num_pions = 0;
            int num_showers = 0;

            for (std::size_t ipfp = 0; ipfp < slice->reco.npfp; ++ipfp) {
                if (int(ipfp) == ipfp_muon)
                    continue;

                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == 1) num_protons++;
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == 2) num_pions++;
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == 3) num_showers++;
            } // loop pfp

            return num_protons > 1 && num_pions == 0 && num_showers == 0;
        });

        const ana::Cut slice_1mu1p ([](const caf::SRSliceProxy *slice) -> bool {
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            if (ipfp_muon == -1) return false; // redundant, btw, but who cares...

            int num_protons = 0;
            int num_pions = 0;
            int num_showers = 0;

            for (std::size_t ipfp = 0; ipfp < slice->reco.npfp; ++ipfp) {
                if (int(ipfp) == ipfp_muon)
                    continue;

                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == 1) num_protons++;
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == 2) num_pions++;
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == 3) num_showers++;
            } // loop pfp

            return num_protons == 1 && num_pions == 0 && num_showers == 0;
        });

        double flashtime = 0;
        const ana::SpillCut spill_CRTPMTNeutrino ([](const caf::SRSpillProxy *spill) -> bool {

            
            for(const auto& crtpmt_match: spill->crtpmt_matches) {
                //Define the interval depending on Data or MC files
                double min_time = -1, max_time = -1;
                
                if (spill->hdr.ismc) {
                    min_time = 0.0; 
                    max_time = 1.6; 
                } else {
                    min_time = -0.4; 
                    max_time = 1.5;
                }

                if (crtpmt_match.flashGateTime > min_time && crtpmt_match.flashGateTime < max_time && crtpmt_match.flashClassification == 0) {
                    flashtime = crtpmt_match.flashGateTime;
                    return true;
                } 
                // if (match.flashGateTime > -0.3 && match.flashGateTime < 1.3 && match.flashClassification == 0) { 
                //     return true; 
                // } 
                
            } // loop spill->crtpmt_matches
            
            return false; 
        }); // const ana::SpillCut CRTPMTNeutrino
    } // namespace reco
} //namespace cuts

namespace vars {
    namespace reco {
        double neutrino_energy_Np (const caf::Proxy<caf::SRSlice> &slice, int ipfp_mu, int dist_emucut) {
            /* Returns the neutrino energy in MeV 
            */

            float p_mu_x = -1, p_mu_y = -1, p_mu_z = -1;
            // float p_p_x = -1, p_p_y = -1, p_p_z = -1;
            // float p_tot_x = -1, p_tot_y = -1, p_tot_z = -1;
            double E_mu = 0, E_p = 0;

            int ipfp_pro = -1;

            p_mu_x = slice.reco.pfp[ipfp_mu].trk.rangeP.p_muon * slice.reco.pfp[ipfp_mu].trk.dir.x;   // Momenta are in GeV
            p_mu_y = slice.reco.pfp[ipfp_mu].trk.rangeP.p_muon * slice.reco.pfp[ipfp_mu].trk.dir.y;   // Momenta are in GeV
            p_mu_z = slice.reco.pfp[ipfp_mu].trk.rangeP.p_muon * slice.reco.pfp[ipfp_mu].trk.dir.z;   // Momenta are in GeV
            
            double p_mu_tot = std::sqrt(p_mu_x * p_mu_x + p_mu_y * p_mu_y + p_mu_z * p_mu_z); // [GeV]
            
            E_mu = 1000. * std::sqrt(p_mu_tot * p_mu_tot + std::pow(particle_data::masses::muon, 2) / (1000. * 1000.));

            for (std::size_t ipfp = 0; ipfp < slice.reco.npfp; ++ipfp) {
                if (int(ipfp) == ipfp_mu)
                    continue;
                if (var_utils::id_pfp(slice, ipfp, dist_emucut) == 1) {
                    TVector3 start_mom;
                    start_mom.SetXYZ(
                        slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.x, 
                        slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.y, 
                        slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.z
                    );
                    E_p += std::sqrt(
                        std::pow(particle_data::masses::proton, 2) + 
                        std::pow(start_mom.Mag() * 1000, 2)
                    ) - particle_data::masses::proton;
                    // ipfp_pro = ipfp;
                } // this pfp is proton-like
            } // loop pfp

            return (E_mu + E_p) / 1000.;
        } // double neutrino_energy_Np

        double neutrino_pT_Np (const caf::Proxy<caf::SRSlice> &islc, int ipfp_mu, int dist_cut) {

            float p_p_x = 0, p_p_y = 0, p_p_z = 0;
            float p_tot_x = 0, p_tot_y = 0, p_tot_z = 0;
            
            int ipfp_pro = -1;

            float p_mu_x = (islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon) * islc.reco.pfp[ipfp_mu].trk.dir.x;   // Momenta are in GeV
            float p_mu_y = (islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon) * islc.reco.pfp[ipfp_mu].trk.dir.y;   // Momenta are in GeV
            float p_mu_z = (islc.reco.pfp[ipfp_mu].trk.rangeP.p_muon) * islc.reco.pfp[ipfp_mu].trk.dir.z;   // Momenta are in GeV

            for (std::size_t ipfp = 0; ipfp < islc.reco.npfp; ++ipfp) {
                if (int(ipfp) == ipfp_mu)
                    continue;
                if (var_utils::id_pfp(islc, ipfp, dist_cut) == 1) {
                    p_p_x += (islc.reco.pfp[ipfp].trk.rangeP.p_proton) * islc.reco.pfp[ipfp].trk.dir.x;
                    p_p_y += (islc.reco.pfp[ipfp].trk.rangeP.p_proton) * islc.reco.pfp[ipfp].trk.dir.y;
                    p_p_z += (islc.reco.pfp[ipfp].trk.rangeP.p_proton) * islc.reco.pfp[ipfp].trk.dir.z;
                    ipfp_pro = ipfp;
                } // this pfp is proton-like
            } // loop pfp

            // if (ipfp_mu != -1 && ipfp_pro != -1) {
                p_tot_x = p_p_x + p_mu_x;
                p_tot_y = p_p_y + p_mu_y;
                p_tot_z = p_p_z + p_mu_z;
            // }

            return std::sqrt( std::pow(p_tot_x, 2) + std::pow(p_tot_y, 2) );
        } // double Transverse_mom_reco_Np

        const ana::Var slice_neutrino_energy_1muNp ([](const caf::SRSliceProxy *slice) -> double {
            /* This slice is the selected slice (so reco 1µNp or truth 1µNp)
             * Here only the computation should be performed, since the cuts will be applied
             * at the Tree/Spectrum stage
            */
            
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            if (ipfp_muon == -1) {
                std::cout << "Did not find any µ" << std::endl;
                return -1; // negative energy backed up by cut
            }

            if(std::isnan(neutrino_energy_Np (*slice, ipfp_muon, var_utils::dist_cut))) {
                std::cout << "Neutrino energy is nan" << std::endl;
                return -1;
            }

            return neutrino_energy_Np (*slice, ipfp_muon, var_utils::dist_cut);
        }); // const ana::Var slice_neutrino_energy_reco_1muNp

        const ana::Var slice_neutrino_pT_1muNp ([](const caf::SRSliceProxy *slice) -> double {
            /* This slice is the selected slice (so reco 1µNp or truth 1µNp)
             * Here only the computation should be performed, since the cuts will be applied
             * at the Tree/Spectrum stage
            */
            
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            if (ipfp_muon == -1) return -1; // negative energy backed up by cut

            return neutrino_pT_Np (*slice, ipfp_muon, var_utils::dist_cut);
        }); // const ana::Var slice_neutrino_energy_reco_1muNp

        const ana::Var slice_muon_hit_completeness ([](const caf::SRSliceProxy *slice) -> double {
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            if (ipfp_muon == -1) return -1; // negative energy backed up by cut

            return slice->reco.pfp[ipfp_muon].trk.truth.bestmatch.hit_completeness;
        }); // const ana::Var slice_muon_hit_completeness

        const ana::Var slice_muon_hit_purity ([](const caf::SRSliceProxy *slice)    -> double {
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            if (ipfp_muon == -1) return -1; // negative energy backed up by cut

            return slice->reco.pfp[ipfp_muon].trk.truth.bestmatch.hit_purity;
        }); // const ana::Var slice_muon_hit_purity
        
        const ana::Var slice_pid_muon_reco_length ([](const caf::SRSliceProxy *slice) -> double {
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            if (ipfp_muon == -1) return -1; // negative energy backed up by cut

            return slice->reco.pfp[ipfp_muon].trk.len;
        }); // const ana::Var slice_pid_muon_reco_length

        const ana::Var slice_pid_muon_true_length ([](const caf::SRSliceProxy *slice) -> double {
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            if (ipfp_muon == -1) return -1; // negative energy backed up by cut

            return slice->reco.pfp[ipfp_muon].trk.truth.p.length;
        }); // const ana::Var slice_pid_muon_true_length

        const ana::Var slice_pid_muon_L_reco_true_ratio ([](const caf::SRSliceProxy *slice) -> double {
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            if (ipfp_muon == -1) return -1; // negative energy backed up by cut

            return slice->reco.pfp[ipfp_muon].trk.len / slice->reco.pfp[ipfp_muon].trk.truth.p.length;
        }); // const ana::Var slice_pid_muon_L_reco_true_ratio

        const ana::Var slice_pid_proton_reco_length ([](const caf::SRSliceProxy *slice) -> double {
            double length = -1;
            for (std::size_t ipfp=0; ipfp<slice->reco.npfp; ++ipfp) {
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) != 1) continue;
                if (slice->reco.pfp[ipfp].trk.len>length)
                    length = slice->reco.pfp[ipfp].trk.len;
            } // pfp loops
            
            return length;
        }); // const ana::Var slice_pid_muon_reco_length

        const ana::Var slice_pid_proton_true_length ([](const caf::SRSliceProxy *slice) -> double {
            double length = -1;
            for (std::size_t ipfp=0; ipfp<slice->reco.npfp; ++ipfp) {
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) != 1) continue;
                if (slice->reco.pfp[ipfp].trk.truth.p.length>length)
                    length = slice->reco.pfp[ipfp].trk.truth.p.length;
            } // pfp loops

            return length;
        }); // const ana::Var slice_pid_muon_true_length

        const ana::Var slice_pid_proton_L_reco_true_ratio ([](const caf::SRSliceProxy *slice) -> double {
            double length = -1;
            double ratio = -1;
            for (std::size_t ipfp=0; ipfp<slice->reco.npfp; ++ipfp) {
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) != 1) continue;
                if (slice->reco.pfp[ipfp].trk.len>length)
                    length = slice->reco.pfp[ipfp].trk.len;
                    ratio = slice->reco.pfp[ipfp].trk.len / slice->reco.pfp[ipfp].trk.truth.p.length;
            } // pfp loops
            
            return ratio;
        }); // const ana::Var slice_pid_proton_L_reco_true_ratio

        const ana::Var slice_proton_hit_completeness ([](const caf::SRSliceProxy *slice) -> double {
            double length = -1;
            int ipfp_proton = -1;
            for (std::size_t ipfp=0; ipfp<slice->reco.npfp; ++ipfp) {
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) != 1) continue;
                if (slice->reco.pfp[ipfp].trk.len>length) {
                    length = slice->reco.pfp[ipfp].trk.len;
                    ipfp_proton = ipfp;
                }
            } // pfp loops

            if (ipfp_proton == -1) return -5;
            return slice->reco.pfp[ipfp_proton].trk.truth.bestmatch.hit_completeness;
        }); // const ana::Var slice_proton_hit_completeness

        const ana::Var slice_proton_hit_purity ([](const caf::SRSliceProxy *slice) -> double {
            double length = -1;
            int ipfp_proton = -1;
            for (std::size_t ipfp=0; ipfp<slice->reco.npfp; ++ipfp) {
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) != 1) continue;
                if (slice->reco.pfp[ipfp].trk.len>length) {
                    length = slice->reco.pfp[ipfp].trk.len;
                    ipfp_proton = ipfp;
                }
            } // pfp loops

            if (ipfp_proton == -1) return -5;
            return slice->reco.pfp[ipfp_proton].trk.truth.bestmatch.hit_purity;
        }); // const ana::Var slice_proton_hit_purity

        const ana::Var slice_muon_momentum_rangeP ([](const caf::SRSliceProxy *slice) -> double {
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            if (ipfp_muon == -1) return -1; // negative energy backed up by cut

            return slice->reco.pfp[ipfp_muon].trk.rangeP.p_muon;
        }); // const ana::Var slice_muon_momentum_rangeP

        const ana::Var slice_proton_momentum_rangeP ([](const caf::SRSliceProxy *slice) -> double {
            double length = -1;
            int ipfp_proton = -1;
            for (std::size_t ipfp=0; ipfp<slice->reco.npfp; ++ipfp) {
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) != 1) continue;
                if (slice->reco.pfp[ipfp].trk.len>length) {
                    length = slice->reco.pfp[ipfp].trk.len;
                    ipfp_proton = ipfp;
                }
            } // pfp loops

            if (ipfp_proton == -1) return -5;
            return slice->reco.pfp[ipfp_proton].trk.rangeP.p_proton;
        }); // const ana::Var slice_proton_momentum_rangeP

        const ana::Var slice_muon_P_reco_true_ratio ([](const caf::SRSliceProxy *slice) -> double {
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            if (ipfp_muon == -1) return -1; // negative energy backed up by cut

            double startp_mag = std::sqrt(
                std::pow(slice->reco.pfp[ipfp_muon].trk.truth.p.startp.x, 2) + 
                std::pow(slice->reco.pfp[ipfp_muon].trk.truth.p.startp.y, 2) + 
                std::pow(slice->reco.pfp[ipfp_muon].trk.truth.p.startp.z, 2)
            );
            return slice->reco.pfp[ipfp_muon].trk.rangeP.p_muon / startp_mag;
        }); // const ana::Var slice_muon_P_reco_true_ratio

        const ana::Var slice_proton_P_reco_true_ratio ([](const caf::SRSliceProxy *slice) -> double {
            double length = -1;
            int ipfp_proton = -1;
            for (std::size_t ipfp=0; ipfp<slice->reco.npfp; ++ipfp) {
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) != 1) continue;
                if (slice->reco.pfp[ipfp].trk.len>length) {
                    length = slice->reco.pfp[ipfp].trk.len;
                    ipfp_proton = ipfp;
                }
            } // pfp loops

            if (ipfp_proton == -1) return -5;

            double startp_mag = std::sqrt(
                std::pow(slice->reco.pfp[ipfp_proton].trk.truth.p.startp.x, 2) + 
                std::pow(slice->reco.pfp[ipfp_proton].trk.truth.p.startp.y, 2) + 
                std::pow(slice->reco.pfp[ipfp_proton].trk.truth.p.startp.z, 2)
            );

            return slice->reco.pfp[ipfp_proton].trk.rangeP.p_proton / startp_mag;
        }); // const ana::Var slice_proton_P_reco_true_ratio

        const ana::Var slice_vertex_difference_z ([](const caf::SRSliceProxy *slice) -> double {

            if (
                std::isnan(slice->vertex.z) &&
                std::isnan(slice->truth.position.z) 
            ) return -1; 

            return std::sqrt(
                std::pow(slice->vertex.z - slice->truth.position.z, 2)
            );

        }); // const ana::Var slice_vertex_difference_z

        const ana::Var slice_vertex_difference_y ([](const caf::SRSliceProxy *slice) -> double {

            if (
                std::isnan(slice->vertex.y) &&
                std::isnan(slice->truth.position.y) 
            ) return -1; 

            return std::sqrt(
                std::pow(slice->vertex.y - slice->truth.position.y, 2)
            );

        }); // const ana::Var slice_vertex_difference_y
        
        const ana::Var slice_vertex_difference_x ([](const caf::SRSliceProxy *slice) -> double {

            if (
                std::isnan(slice->vertex.x) &&
                std::isnan(slice->truth.position.x) 
            ) return -1; 

            return std::sqrt(
                std::pow(slice->vertex.x - slice->truth.position.x, 2)
            );

        }); // const ana::Var slice_vertex_difference_x

        const ana::Var slice_vertex_difference ([](const caf::SRSliceProxy *slice) -> double {

            if (
                std::isnan(slice->vertex.x) &&
                std::isnan(slice->vertex.y) &&
                std::isnan(slice->vertex.z) &&
                std::isnan(slice->truth.position.x) &&
                std::isnan(slice->truth.position.y) &&
                std::isnan(slice->truth.position.z) 
            ) return -1; 

            return std::sqrt(
                std::pow(slice->vertex.x - slice->truth.position.x, 2) + 
                std::pow(slice->vertex.y - slice->truth.position.y, 2) + 
                std::pow(slice->vertex.z - slice->truth.position.z, 2)
            );

        }); // const ana::Var slice_vertex_difference

        const ana::Var slice_muon_R ([](const caf::SRSliceProxy *slice) -> double {
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            if (ipfp_muon == -1) return -1;

            double E_dep_true = slice->reco.pfp.at(ipfp_muon).trk.truth.p.startE - slice->reco.pfp.at(ipfp_muon).trk.truth.p.endE;
            double E_true_in_hits = slice->reco.pfp.at(ipfp_muon).trk.truth.bestmatch.energy/3. * slice->reco.pfp.at(ipfp_muon).trk.truth.bestmatch.energy_completeness;

            return E_true_in_hits / E_dep_true;

        });

        const ana::Var slice_leading_proton_R ([](const caf::SRSliceProxy *slice) -> double {
            double length = -1;
            int ipfp_proton = -1;
            for (std::size_t ipfp=0; ipfp<slice->reco.npfp; ++ipfp) {
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) != 1) continue;
                if (slice->reco.pfp[ipfp].trk.len>length) {
                    length = slice->reco.pfp[ipfp].trk.len;
                    ipfp_proton = ipfp;
                }
            } // pfp loops

            if (ipfp_proton == -1) return -5;

            double E_dep_true = slice->reco.pfp.at(ipfp_proton).trk.truth.p.startE - slice->reco.pfp.at(ipfp_proton).trk.truth.p.endE;
            double E_true_in_hits = slice->reco.pfp.at(ipfp_proton).trk.truth.bestmatch.energy/3. * slice->reco.pfp.at(ipfp_proton).trk.truth.bestmatch.energy_completeness;

            return E_true_in_hits / E_dep_true;

        });
    } // namespace reco

    namespace truth {
        const ana::Var slice_neutrino_energy ([](const caf::SRSliceProxy *slice) -> double {
            /* This slice is the selected slice (so reco 1µNp or truth 1µNp)
             * Here only the computation should be performed, since the cuts will be applied
             * at the Tree/Spectrum stage
            */
                
            if (std::isnan(slice->truth.E)) return -1;
            return slice->truth.E; // true neutrino energy in GeV
        }); // const ana::Var slice_neutrino_energy

        const ana::Var slice_neutrino_dE ([](const caf::SRSliceProxy *slice) -> double {
            
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            if (ipfp_muon == -1) {
                std::cout << "Did not find any µ" << std::endl;
                return -9999; // negative energy backed up by cut
            }

            double reco_E = vars::reco::neutrino_energy_Np (*slice, ipfp_muon, var_utils::dist_cut);
            if(std::isnan(reco_E)) {
                std::cout << "Neutrino energy is nan" << std::endl;
                return -9999;
            }

            if (std::isnan(slice->truth.E)) return -9999;

            return (reco_E - slice->truth.E)/slice->truth.E;
        }); // const ana::Var slice_neutrino_dE
    } // namespace truth
} // namespace vars

namespace var_utils {
    particle_data::int_type classification_type (const caf::SRSpillProxy *spill, const caf::SRSliceProxy *slice) {
        /* Utility to help classify different true interaction types
         *  - 1µNp are particle_data::int_type::true_visible_1muNp
         *  - 1µ1p are particle_data::int_type::true_visible_1mu1p
         *  - unclassified
        */

        int num_protons_above50 = 0;
        int num_muons = 0;
        double length_muon = 0;
        double dep_E = 0;

        int G4ID_parent;
        int use_plane = 2;

        for (auto const& prim: slice->truth.prim) {

            dep_E = 0;
            
            if (prim.G4ID < 0)
                continue;
            
            if (prim.cryostat < 0)
                continue;

            if (std::abs(prim.pdg) == 211) {
                return particle_data::int_type::unclassified;
            } // found charged pion, returning

            /* Looking at neutral pions: trickier
             * Reject if any of daughter photon is > 25 MeV
            */ 
            if (std::abs(prim.pdg) == 111 && prim.daughters.size() >= 1) {
                for (auto const& true_particle: spill->true_particles) {
                    G4ID_parent = true_particle.parent;

                    if (
                        G4ID_parent == prim.G4ID && std::abs(true_particle.pdg) == 22 &&
                        true_particle.plane[prim.cryostat][use_plane].visE * 1000 > particle_data::minimum_gamma_MeV
                    ){
                        return particle_data::int_type::unclassified;
                    }
                } // loop spill->true_particles
            } // neutral pion found

            if (std::abs(prim.pdg) == 13) {
                num_muons += 1;
                length_muon = prim.length;
            } // found muon
            
            if (std::abs(prim.pdg) == 22) {
                if (prim.daughters.size() >= 1) {
                    for (auto const& true_particle: spill->true_particles) {
                        G4ID_parent = true_particle.parent;
                        if (G4ID_parent == prim.G4ID)
                            dep_E += true_particle.plane[prim.cryostat][use_plane].visE * 1000;
                    } // loop trough true_particles
                } // photon has daughters
                dep_E += prim.plane[prim.cryostat][use_plane].visE * 1000;
            } // found photon

            // if primary is photon with dep energy > 25 MeV skip...
            if (std::abs(prim.pdg) == 22 && dep_E > particle_data::minimum_gamma_MeV)
                return particle_data::int_type::unclassified;

            if (std::abs(prim.pdg) == 2212) {
                if (prim.daughters.size() >= 1) {
                    for(auto const& true_particle: spill->true_particles) {
                        G4ID_parent = true_particle.parent;
                        if (G4ID_parent == prim.G4ID) 
                            dep_E += true_particle.plane[prim.cryostat][use_plane].visE * 1000;
                    } // loop trough true_particles
                } // proton has daughters
                dep_E += prim.plane[prim.cryostat][use_plane].visE * 1000;
            } // found proton

            if (std::abs(prim.pdg) == 2212 && dep_E > 50.0)
                num_protons_above50 += 1;
            
            dep_E = 0;
        } // loop slice->truth.prim

        if (
            num_muons == 1 && num_protons_above50 == 1 && 
            length_muon > 50 && cuts::all_trk_contained_truth(spill, slice)
        ) return particle_data::int_type::true_visible_1mu1p;

        if (
            num_muons == 1 && num_protons_above50 > 1 && 
            length_muon > 50 && cuts::all_trk_contained_truth(spill, slice)
        ) return particle_data::int_type::true_visible_1muNp;

        return particle_data::int_type::unclassified;
    } // int classification_type

    particle_data::int_type classification_type_MC (
        const caf::SRSpillProxy *spill, const caf::Proxy<caf::SRTrueInteraction> &nu
    ) {
        /* Utility to help classify different true interaction types
         *  - 1µNp are particle_data::int_type::true_visible_1muNp
         *  - 1µ1p are particle_data::int_type::true_visible_1mu1p
         *  - all others are particle_data::int_type::unclassified
        */

        int num_protons_above50 = 0;
        int num_muons = 0;
        double length_muon = 0;
        double dep_E = 0;

        int G4ID_parent;
        int use_plane = 2;

        for (auto const& prim: nu.prim) {

            dep_E = 0;
            
            if (prim.G4ID < 0)
                continue;

            if (prim.cryostat < 0)
                continue;

            if (std::abs(prim.pdg) == 211) {
                return particle_data::int_type::unclassified;
            } // found charged pion, returning

            /* Looking at neutral pions: trickier
             * Reject if any of daughter photon is > 25 MeV
            */ 
            if (std::abs(prim.pdg) == 111 && prim.daughters.size() >= 1) {
                for (auto const& true_particle: spill->true_particles) {
                    G4ID_parent = true_particle.parent;

                    if (
                        G4ID_parent == prim.G4ID && std::abs(true_particle.pdg) == 22 &&
                        true_particle.plane[prim.cryostat][use_plane].visE * 1000 > particle_data::minimum_gamma_MeV
                    ){
                        return particle_data::int_type::unclassified;
                    }
                } // loop spill->true_particles
            } // neutral pion found

            if (std::abs(prim.pdg) == 13) {
                num_muons += 1;
                length_muon = prim.length;
            } // found muon
            
            if (std::abs(prim.pdg) == 22) {
                if (prim.daughters.size() >= 1) {
                    for (auto const& true_particle: spill->true_particles) {
                        G4ID_parent = true_particle.parent;
                        if (G4ID_parent == prim.G4ID)
                            dep_E += true_particle.plane[prim.cryostat][use_plane].visE * 1000;
                    } // loop trough true_particles
                } // photon has daughters
                dep_E += prim.plane[prim.cryostat][use_plane].visE * 1000;
            } // found photon

            // if primary is photon with dep energy > 25 MeV skip...
            if (std::abs(prim.pdg) == 22 && dep_E > particle_data::minimum_gamma_MeV)
                return particle_data::int_type::unclassified;

            if (std::abs(prim.pdg) == 2212) {
                if (prim.daughters.size() >= 1) {
                    for(auto const& true_particle: spill->true_particles) {
                        G4ID_parent = true_particle.parent;
                        if (G4ID_parent == prim.G4ID) 
                            dep_E += true_particle.plane[prim.cryostat][use_plane].visE * 1000;
                    } // loop trough true_particles
                } // proton has daughters
                dep_E += prim.plane[prim.cryostat][use_plane].visE * 1000;
            } // found proton

            if (std::abs(prim.pdg) == 2212 && dep_E > 50.0)
                num_protons_above50 += 1;
            
            dep_E = 0;
        } // loop slice->truth.prim

        if (
            num_muons == 1 && num_protons_above50 == 1 && 
            length_muon > 50 && cuts::all_trk_contained_MC(spill, nu)
        ) return particle_data::int_type::true_visible_1mu1p;

        if (
            num_muons == 1 && num_protons_above50 > 1 && 
            length_muon > 50 && cuts::all_trk_contained_MC(spill, nu)
        ) return particle_data::int_type::true_visible_1muNp;

        return particle_data::int_type::unclassified;
    } // int classification_type
} // namespace var_utils

#endif