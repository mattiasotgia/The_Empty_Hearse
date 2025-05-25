#ifndef SELECTION_H
#define SELECTION_H

#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"

#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "TFile.h"
#include "TVector3.h"
#include "TProfile.h"
#include "TMath.h"

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

    enum int_type_t {
        true_visible_1muNp,
        true_visible_1mu1p,
        true_visible_1mu2p,
        true_visible_1mu3p,
        true_visible_1mu1pi,
        true_visible_1mu2pi,
        true_visible_1mu3pi,
        true_visible_1muNpi,
        true_visible_2mu,
        true_visible_2p,
        true_visible_1muShortNp,
        true_visible_1pi1p, 
        true_visible_1mu1p1pi0,
        true_visible_1mu1pNpi0,
        unclassified
    };

    enum particle_t {
        proton,
        pion,
        shower,
        muon,
        electron,
        not_primary,
        too_far,
        low_energy,
        undefined
    };

    double minimum_gamma_MeV = 25.;
    double GeV = 1000.;
} // namespace particle_data

namespace var_utils {
    int find_muon(const caf::SRSliceProxy&, double);
    double dist_cut = 10.;
}

namespace cuts {
    bool in_contained(double, double, double, double);
    namespace reco {
        const ana::Cut slice_at_least_mu ([](const caf::SRSliceProxy *slice) -> bool {
            /* We should at least found one muon, otherwise something bad happens
            */

            // Using dist_cut = 10 [cm]
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut); 
            if (ipfp_muon == -1) return false; // no muon found ahaha
            return true;
        }); // const ana::Cut slice_at_least_mu
    }
}

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

    values_minmax barycenterFM_deltaZ_Trigger = {0., 100.};

    std::string dEdx_temp = 
        "/exp/icarus/app/users/msotgia/analysis/sbnana_v09_93_01_thesis_analysis/analysis/dEdxrestemplates.root";
    TFile* file = TFile::Open(dEdx_temp.c_str());

    auto dedx_range_pro = (TProfile *)file->Get("dedx_range_pro");
    auto dedx_range_ka  = (TProfile *)file->Get("dedx_range_ka");
    auto dedx_range_pi  = (TProfile *)file->Get("dedx_range_pi");
    auto dedx_range_mu  = (TProfile *)file->Get("dedx_range_mu");

    enum cut_type_t {
        RECO,
        TRUE_1muN1p,
        TRUE_1muNp,
        TRUE_1mu1p,
        TRUE_1mu2p,
        TRUE_1mu3p,
        BOTH_1muN1p,
        BOTH_1muNp,
        BOTH_1mu1p,
        BOTH_1mu2p,
        BOTH_1mu3p,
        BOTH_1muN3p,
        BOTH_1mu1pi,
        BOTH_1mu2pi,
        BOTH_1mu3pi,
        BOTH_1muNpi,
        BOTH_2mu,
        BOTH_2p,
        BOTH_1muShortNp,
        BOTH_1pi1p,
        BOTH_1mu1p1pi0,
        BOTH_1mu1pNpi0,
        MC_1muNp,
        MC_1mu1p
    };

    particle_data::int_type_t classification_type (const caf::SRSpillProxy*, const caf::SRSliceProxy*);
    template<class R, class T, class A>
    const R make_spill_from_slice (
        A __def_ret, 
        const T &slice_var, 
        const ana::Cut &reco_cut = ana::kNoCut,  
        cut_type_t what_to_cut_on = cut_type_t::RECO, 
        const ana::Cut &truth_cut = ana::kNoCut
    ) {
        return R([=](const caf::SRSpillProxy *spill) -> A {
            int selected_slices = 0;
            A slice_value = __def_ret;
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
                
                // if (!cuts::reco::slice_at_least_mu(&slice))
                //     continue;

                slice_value = slice_var(&slice);

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
                if (bincpro < 1e-6) { // for 0 bin content, using neighboring bins
                    bincpro = 
                        (dedx_range_pro->GetBinContent(bin - 1) + dedx_range_pro->GetBinContent(bin + 1)) / 2;
                }
                
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
    
    int find_muon (const caf::SRSliceProxy &slice, double dist_mucut) {

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
                slice.reco.pfp[ipfp].trk.len > 50           && // The selection of the longest muon being 50+ cm should not be here? // separate Cut below
                (reco_vtx - reco_start).Mag() < dist_mucut  && 
                chi2_values.muon < 30                       && // chi2 cuts on proton and muon 
                chi2_values.proton > 60                     && // chi2 cuts on proton and muon 
                cuts::in_contained (slice.reco.pfp[ipfp].trk.end.x, slice.reco.pfp[ipfp].trk.end.y, slice.reco.pfp[ipfp].trk.end.z, 5.) && // separate Cut below
                slice.reco.pfp[ipfp].trk.end.x * slice.vertex.x > 0 &&                                                               // separate Cut below
                slice.reco.pfp[ipfp].parent_is_primary
            ) {
                max_length = slice.reco.pfp[ipfp].trk.len;
                ipfp_mu = ipfp;
            }
        } // loop of pfp to find muon
        return ipfp_mu;
    } // int find_muon

    particle_data::particle_t id_pfp (const caf::SRSliceProxy &slice, int ipfp, int dist_cut) {
        /* This utility select the correct particle of the event slice
         * return 1 PROTONS
         * return 2 PIONS
         * return 3 SHOWER
         * return 9 other -> nan, not primary, too far, below energy threshold...
        */

        TVector3 rec_vtx;
        rec_vtx.SetXYZ(slice.vertex.x, slice.vertex.y, slice.vertex.z);

        if (!(slice.reco.pfp[ipfp].parent_is_primary))
            return particle_data::particle_t::not_primary;

        // Sanity checks :)
        if (
            std::isnan(slice.reco.pfp[ipfp].trk.start.x)    || 
            std::isnan(slice.reco.pfp[ipfp].trk.end.x)      || 
            std::isnan(slice.reco.pfp[ipfp].trk.len)
        ) return particle_data::particle_t::undefined;
        
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
            return particle_data::particle_t::too_far;

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
            ) return particle_data::particle_t::undefined;

            if (
                std::isnan(slice.reco.pfp[ipfp].trk.start.y) || 
                std::isnan(slice.reco.pfp[ipfp].trk.start.z) || 
                std::isnan(slice.reco.pfp[ipfp].trk.end.y)   || 
                std::isnan(slice.reco.pfp[ipfp].trk.end.z)
            ) return particle_data::particle_t::undefined;

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
                std::sqrt ( 
                    std::pow(particle_data::masses::pion, 2) + std::pow(start_mom_V3.Mag() * particle_data::GeV, 2) 
                ) - particle_data::masses::pion >= 25.0 && 
                slice.reco.pfp[ipfp].parent_is_primary
            )  return particle_data::particle_t::pion;

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
                std::sqrt ( 
                    std::pow(particle_data::masses::proton, 2) + std::pow(start_mom_V3.Mag() * particle_data::GeV, 2) 
                ) - particle_data::masses::proton >= 50.0 && 
                slice.reco.pfp[ipfp].parent_is_primary
            ) return particle_data::particle_t::proton;

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
                    std::sqrt( 
                        std::pow(particle_data::masses::proton, 2) + std::pow(start_mom_V3_2.Mag() * particle_data::GeV, 2)
                    ) - particle_data::masses::proton >= 50.0 && 
                    (rec_vtx - start).Mag() < dist_cut && 
                    slice.reco.pfp[ipfp].parent_is_primary
                ) return particle_data::particle_t::proton;
            }

            if (!(slice.reco.pfp[ipfp].trackScore >= 0.4 && chi2_values.proton < 100)) {

                // int use_plane2 = slice.reco.pfp[ipfp].trk.calo[2].nhit>slice.reco.pfp[ipfp].trk.calo[1].nhit ? 2:1;
                int use_plane2 = 2;

                if (std::isnan(slice.reco.pfp[ipfp].shw.plane[use_plane2].energy))
                    return particle_data::particle_t::undefined;
                if (slice.reco.pfp[ipfp].shw.plane[use_plane2].energy * particle_data::GeV < 25.0)
                    return particle_data::particle_t::low_energy;
                if (
                    slice.reco.pfp[ipfp].shw.plane[use_plane2].energy * particle_data::GeV > 25.0 && 
                    slice.reco.pfp[ipfp].parent_is_primary)
                    return particle_data::particle_t::shower;
            }
        }
        return particle_data::particle_t::undefined;
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
                std::abs(prim.pdg) == 13   ||   // mu
                std::abs(prim.pdg) == 2212 ||   // p
                std::abs(prim.pdg) == 211  ||   // pi
                std::abs(prim.pdg) == 11        // e 
            ) {
                if (!in_contained(prim.end.x, prim.end.y, prim.end.z)) {
                    return false;
                }
            }
                    

            if (prim.daughters.size() > 0) {
                int parent_g4id;

                for (auto const& true_particle: spill->true_particles) {
                    parent_g4id = true_particle.parent;

                    if (parent_g4id == prim.G4ID) {
                        if (
                            std::abs(true_particle.pdg) == 13   ||  // mu
                            std::abs(true_particle.pdg) == 2212 ||  // p
                            std::abs(true_particle.pdg) == 211  ||  // pi
                            std::abs(true_particle.pdg) == 11       // e 
                        ) {
                            if (!in_contained(true_particle.end.x, true_particle.end.y, true_particle.end.z)) {
                                return false;
                            }
                        }
                    }
                    
                    /* L’idea e’ mettere soglia su quello che deposita energia…. 
                     * Quindi tutto quello che sia carico e che di solito abbiamo
                     * Per dire in qualche punto ci deve essere anche una soglia 
                     * sui fotoni singoli, per cui se hanno troppa energia depositata 
                     * possono essere identificati (non avevamo messo soglia sui pi0 
                     * pero si sui due fotoni del pi0)
                    */

                } // loop spill->true_particles
            } // has daughters
        } // loop slc.truth.prim
        return true;
    } // bool all_trk_contained_truth

    bool all_trk_contained_MC (const caf::SRSpillProxy *spill, const caf::Proxy<caf::SRTrueInteraction> &nu) {

        for (auto const& prim: nu.prim) {
            if (prim.G4ID < 0)     continue;
            if (prim.cryostat < 0) continue;
            if (
                std::abs(prim.pdg) == 13   ||   // mu
                std::abs(prim.pdg) == 2212 ||   // p
                std::abs(prim.pdg) == 211  ||   // pi
                std::abs(prim.pdg) == 11        // e 
            ) {
                if (!in_contained(prim.end.x, prim.end.y, prim.end.z)) {
                    return false;
                }
            }
                    

            if (prim.daughters.size() > 0) {
                int parent_g4id;

                for (auto const& true_particle: spill->true_particles) {
                    parent_g4id = true_particle.parent;

                    if (parent_g4id == prim.G4ID) {
                        if (
                            std::abs(true_particle.pdg) == 13   ||  // mu
                            std::abs(true_particle.pdg) == 2212 ||  // p
                            std::abs(true_particle.pdg) == 211  ||  // pi
                            std::abs(true_particle.pdg) == 11       // e 
                        ) {
                            if (!in_contained(true_particle.end.x, true_particle.end.y, true_particle.end.z)) {
                                return false;
                            }
                        }
                    }
                    
                    /* L’idea e’ mettere soglia su quello che deposita energia…. 
                     * Quindi tutto quello che sia carico e che di solito abbiamo
                     * Per dire in qualche punto ci deve essere anche una soglia 
                     * sui fotoni singoli, per cui se hanno troppa energia depositata 
                     * possono essere identificati (non avevamo messo soglia sui pi0 
                     * pero si sui due fotoni del pi0)
                    */

                } // loop spill->true_particles
            } // has daughters
        } // loop slc.truth.prim
        return true;
    } // bool all_trk_contained_MC

    namespace truth {
        const ana::Cut slice_numuCC (ana::kIsNumuCC);

        const ana::Cut slice_numuNC ([](const caf::SRSliceProxy *slice) -> bool {
            return slice->truth.isnc && slice->truth.pdg == 14;
        });

        const ana::Cut slice_nueCC ([](const caf::SRSliceProxy *slice) -> bool {
            return slice->truth.iscc && slice->truth.pdg == 12;
        });

        // Genie_interaction_mode_
        const ana::Cut slice_QE ([](const caf::SRSliceProxy *slice) -> bool {
            return slice->truth.genie_mode == caf::genie_interaction_mode_::kQE;
        });

        const ana::Cut slice_Res ([](const caf::SRSliceProxy *slice) -> bool {
            return slice->truth.genie_mode == caf::genie_interaction_mode_::kRes;
        });

        const ana::Cut slice_DIS ([](const caf::SRSliceProxy *slice) -> bool {
            return slice->truth.genie_mode == caf::genie_interaction_mode_::kDIS;
        });

        const ana::Cut slice_Coh ([](const caf::SRSliceProxy *slice) -> bool {
            return slice->truth.genie_mode == caf::genie_interaction_mode_::kCoh;
        });

        const ana::Cut slice_CohElastic ([](const caf::SRSliceProxy *slice) -> bool {
            return slice->truth.genie_mode == caf::genie_interaction_mode_::kCohElastic;
        });

        const ana::Cut slice_ElectronScattering ([](const caf::SRSliceProxy *slice) -> bool {
            return slice->truth.genie_mode == caf::genie_interaction_mode_::kElectronScattering;
        });

        const ana::Cut slice_IMDAnnihilation ([](const caf::SRSliceProxy *slice) -> bool {
            return slice->truth.genie_mode == caf::genie_interaction_mode_::kIMDAnnihilation;
        });

        const ana::Cut slice_InverseBetaDecay ([](const caf::SRSliceProxy *slice) -> bool {
            return slice->truth.genie_mode == caf::genie_interaction_mode_::kInverseBetaDecay;
        });

        const ana::Cut slice_GlashowResonance ([](const caf::SRSliceProxy *slice) -> bool {
            return slice->truth.genie_mode == caf::genie_interaction_mode_::kGlashowResonance;
        });
        
        const ana::Cut slice_AMNuGamma ([](const caf::SRSliceProxy *slice) -> bool {
            return slice->truth.genie_mode == caf::genie_interaction_mode_::kAMNuGamma;
        });
        
        const ana::Cut slice_MEC ([](const caf::SRSliceProxy *slice) -> bool {
            return slice->truth.genie_mode == caf::genie_interaction_mode_::kMEC;
        });

        const ana::Cut slice_Diffractive ([](const caf::SRSliceProxy *slice) -> bool {
            return slice->truth.genie_mode == caf::genie_interaction_mode_::kDiffractive;
        });
        
        const ana::Cut slice_EM ([](const caf::SRSliceProxy *slice) -> bool {
            return slice->truth.genie_mode == caf::genie_interaction_mode_::kEM;
        });

        const ana::Cut slice_WeakMix ([](const caf::SRSliceProxy *slice) -> bool {
            return slice->truth.genie_mode == caf::genie_interaction_mode_::kWeakMix;
        });
        // --------------------------

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

        const ana::Cut slice_1mu_only ([](const caf::SRSliceProxy *slice) -> bool {
            int muon_n = 0;
            for (auto const& prim: slice->truth.prim) {
                if (std::abs(prim.pdg) == 13) muon_n ++;
            } // loop slice->truth.prim
            return muon_n == 1;
        }); // const ana::Cut slice_1mu_only
    } // namespace truth

    namespace reco {
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

                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == particle_data::particle_t::proton)  num_protons++;
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == particle_data::particle_t::pion)    num_pions++;
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == particle_data::particle_t::shower)  num_showers++;
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

                    if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == particle_data::particle_t::proton)  num_protons++;
                    if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == particle_data::particle_t::pion)    num_pions++;
                    if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == particle_data::particle_t::shower)  num_showers++;
            } // loop pfp

            return num_protons == 1 && num_pions == 0 && num_showers == 0;
        });

        const ana::Cut slice_1mu2p ([](const caf::SRSliceProxy *slice) -> bool {
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            if (ipfp_muon == -1) return false; // redundant, btw, but who cares...

            int num_protons = 0;
            int num_pions = 0;
            int num_showers = 0;

            for (std::size_t ipfp = 0; ipfp < slice->reco.npfp; ++ipfp) {
                if (int(ipfp) == ipfp_muon)
                    continue;

                    if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == particle_data::particle_t::proton)  num_protons++;
                    if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == particle_data::particle_t::pion)    num_pions++;
                    if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == particle_data::particle_t::shower)  num_showers++;
            } // loop pfp

            return num_protons == 2 && num_pions == 0 && num_showers == 0;
        });

        const ana::Cut slice_1mu3p ([](const caf::SRSliceProxy *slice) -> bool {
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            if (ipfp_muon == -1) return false; // redundant, btw, but who cares...

            int num_protons = 0;
            int num_pions = 0;
            int num_showers = 0;

            for (std::size_t ipfp = 0; ipfp < slice->reco.npfp; ++ipfp) {
                if (int(ipfp) == ipfp_muon)
                    continue;

                    if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == particle_data::particle_t::proton)  num_protons++;
                    if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == particle_data::particle_t::pion)    num_pions++;
                    if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) == particle_data::particle_t::shower)  num_showers++;
            } // loop pfp

            return num_protons == 3 && num_pions == 0 && num_showers == 0;
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
            
            E_mu = particle_data::GeV * std::sqrt(p_mu_tot * p_mu_tot + std::pow(particle_data::masses::muon, 2) / (particle_data::GeV * particle_data::GeV));

            for (std::size_t ipfp = 0; ipfp < slice.reco.npfp; ++ipfp) {
                if (int(ipfp) == ipfp_mu)
                    continue;
                if (var_utils::id_pfp(slice, ipfp, dist_emucut) == particle_data::particle_t::proton) {
                    TVector3 start_mom;
                    start_mom.SetXYZ(
                        slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.x, 
                        slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.y, 
                        slice.reco.pfp[ipfp].trk.rangeP.p_proton * slice.reco.pfp[ipfp].trk.dir.z
                    );
                    E_p += std::sqrt(
                        std::pow(particle_data::masses::proton, 2) + 
                        std::pow(start_mom.Mag() * particle_data::GeV, 2)
                    ) - particle_data::masses::proton;
                    // ipfp_pro = ipfp;
                } // this pfp is proton-like
            } // loop pfp

            return (E_mu + E_p) / particle_data::GeV;
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
                if (var_utils::id_pfp(islc, ipfp, dist_cut) == particle_data::particle_t::proton) {
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
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            for (std::size_t ipfp=0; ipfp<slice->reco.npfp; ++ipfp) {
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) != particle_data::particle_t::proton) continue;
                if (ipfp == ipfp_muon)
                    continue;
                if (slice->reco.pfp[ipfp].trk.len>length)
                    length = slice->reco.pfp[ipfp].trk.len;
            } // pfp loops
            
            return length;
        }); // const ana::Var slice_pid_muon_reco_length

        const ana::Var slice_pid_proton_true_length ([](const caf::SRSliceProxy *slice) -> double {
            double length = -1;
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            for (std::size_t ipfp=0; ipfp<slice->reco.npfp; ++ipfp) {
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) != particle_data::particle_t::proton) continue;
                if (ipfp == ipfp_muon)
                    continue;
                if (slice->reco.pfp[ipfp].trk.truth.p.length>length)
                    length = slice->reco.pfp[ipfp].trk.truth.p.length;
            } // pfp loops

            return length;
        }); // const ana::Var slice_pid_muon_true_length

        const ana::Var slice_pid_proton_L_reco_true_ratio ([](const caf::SRSliceProxy *slice) -> double {
            double length = -1;
            double ratio = -1;
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            for (std::size_t ipfp=0; ipfp<slice->reco.npfp; ++ipfp) {
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) != particle_data::particle_t::proton) continue;
                if (ipfp == ipfp_muon)
                    continue;
                if (slice->reco.pfp[ipfp].trk.len>length)
                    length = slice->reco.pfp[ipfp].trk.len;
                    ratio = slice->reco.pfp[ipfp].trk.len / slice->reco.pfp[ipfp].trk.truth.p.length;
            } // pfp loops
            
            return ratio;
        }); // const ana::Var slice_pid_proton_L_reco_true_ratio

        const ana::Var slice_proton_hit_completeness ([](const caf::SRSliceProxy *slice) -> double {
            double length = -1;
            int ipfp_proton = -1;
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            for (std::size_t ipfp=0; ipfp<slice->reco.npfp; ++ipfp) {
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) != particle_data::particle_t::proton) continue;
                if (ipfp == ipfp_muon)
                    continue;
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
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            for (std::size_t ipfp=0; ipfp<slice->reco.npfp; ++ipfp) {
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) != particle_data::particle_t::proton) continue;
                if (ipfp == ipfp_muon)
                    continue;
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
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            for (std::size_t ipfp=0; ipfp<slice->reco.npfp; ++ipfp) {
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) != particle_data::particle_t::proton) continue;
                if (ipfp == ipfp_muon)
                    continue;
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
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            for (std::size_t ipfp=0; ipfp<slice->reco.npfp; ++ipfp) {
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) != particle_data::particle_t::proton) continue;
                if (ipfp == ipfp_muon)
                    continue;
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
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            for (std::size_t ipfp=0; ipfp<slice->reco.npfp; ++ipfp) {
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) != particle_data::particle_t::proton) continue;
                if (ipfp_proton == ipfp_muon)
                    continue;
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

        const ana::Var slice_CT3D_rangeP_muon_proton ([](const caf::SRSliceProxy *slice) -> double {
            double length = -1;
            int ipfp_proton = -1;
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            
            for (std::size_t ipfp = 0; ipfp < slice->reco.npfp; ++ipfp)
            {
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) != particle_data::particle_t::proton)
                    continue;
                if (ipfp_proton == ipfp_muon)
                    continue;
                if (slice->reco.pfp[ipfp].trk.len > length)
                {
                    length = slice->reco.pfp[ipfp].trk.len;
                    ipfp_proton = ipfp;
                }
            } // pfp loops
            
            if (ipfp_muon == -1 || ipfp_proton == -1)
                return -9999;
            
            TVector3 muon_p, proton_p;
            muon_p.SetXYZ(
                slice->reco.pfp.at(ipfp_muon).trk.rangeP.p_muon * slice->reco.pfp.at(ipfp_muon).trk.dir.x,
                slice->reco.pfp.at(ipfp_muon).trk.rangeP.p_muon * slice->reco.pfp.at(ipfp_muon).trk.dir.y,
                slice->reco.pfp.at(ipfp_muon).trk.rangeP.p_muon * slice->reco.pfp.at(ipfp_muon).trk.dir.z
            );
            proton_p.SetXYZ(
                slice->reco.pfp.at(ipfp_proton).trk.rangeP.p_proton * slice->reco.pfp.at(ipfp_proton).trk.dir.x,
                slice->reco.pfp.at(ipfp_proton).trk.rangeP.p_proton * slice->reco.pfp.at(ipfp_proton).trk.dir.y,
                slice->reco.pfp.at(ipfp_proton).trk.rangeP.p_proton * slice->reco.pfp.at(ipfp_proton).trk.dir.z
            );
            return TMath::Cos(muon_p.Angle(proton_p));
        }); // const ana::Var slice_CT3D_muon_leading_proton

        const ana::Var slice_CT3D_trueP_muon_proton ([](const caf::SRSliceProxy *slice) -> double {
            double length = -1;
            int ipfp_proton = -1;
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            
            for (std::size_t ipfp = 0; ipfp < slice->reco.npfp; ++ipfp)
            {
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) != particle_data::particle_t::proton)
                    continue;
                if (ipfp_proton == ipfp_muon)
                    continue;
                if (slice->reco.pfp[ipfp].trk.len > length)
                {
                    length = slice->reco.pfp[ipfp].trk.len;
                    ipfp_proton = ipfp;
                }
            } // pfp loops
            
            if (ipfp_muon == -1 || ipfp_proton == -1)
                return -9999;
            
            TVector3 muon_p, proton_p;
            muon_p.SetXYZ(
                slice->reco.pfp.at(ipfp_muon).trk.truth.p.startp.x,
                slice->reco.pfp.at(ipfp_muon).trk.truth.p.startp.y,
                slice->reco.pfp.at(ipfp_muon).trk.truth.p.startp.z
            );
            proton_p.SetXYZ(
                slice->reco.pfp.at(ipfp_proton).trk.truth.p.startp.x,
                slice->reco.pfp.at(ipfp_proton).trk.truth.p.startp.y,
                slice->reco.pfp.at(ipfp_proton).trk.truth.p.startp.z
            );
            return TMath::Cos(muon_p.Angle(proton_p));
        }); // const ana::Var slice_CT3D_trueP_muon_leading_proton

        std::vector<int> ordered_protons_by_length (const caf::SRSliceProxy *slice) {
            std::vector<std::pair<int, double>> ipfp_length_pairs_all;
            std::vector<int> ordered_ipfp;

            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            
            for (std::size_t ipfp = 0; ipfp < slice->reco.npfp; ++ipfp)
            {
                if (var_utils::id_pfp(*slice, ipfp, var_utils::dist_cut) != particle_data::particle_t::proton)
                    continue;
                if (ipfp == ipfp_muon)
                    continue;
            
                // Add all the pairs (ipfp, lenght)
                ipfp_length_pairs_all.emplace_back(ipfp, slice->reco.pfp[ipfp].trk.len);
    
            } // loop ipfp<slice->npfp
    
            std::sort(ipfp_length_pairs_all.begin(), ipfp_length_pairs_all.end(), 
                [](const std::pair<int, double> &a, const std::pair<int, double> &b) { return a.second > b.second; }
            );
    
            for (auto const& ipfp_length: ipfp_length_pairs_all) 
                ordered_ipfp.push_back(ipfp_length.first);
    
            return ordered_ipfp;
        } // std::vector<int> ordered_pfps_by_length

        const ana::Var slice_muon_track_score ([](const caf::SRSliceProxy *slice) -> double {
            int ipfp_muon = var_utils::find_muon(*slice, var_utils::dist_cut);
            if (ipfp_muon == -1) return -1;

            return slice->reco.pfp[ipfp_muon].trackScore;
        });

        const ana::Var slice_proton_track_score ([](const caf::SRSliceProxy *slice) -> double {
            std::vector<int> ipfps = ordered_protons_by_length(slice);
            if (ipfps.size() == 0) return -1;
            int ipfp_proton = ipfps.at(0);
            if (ipfp_proton == -1) return -1;

            return slice->reco.pfp[ipfp_proton].trackScore;
        });

        const ana::Var slice_second_proton_track_score ([](const caf::SRSliceProxy *slice) -> double {
            std::vector<int> ipfps = ordered_protons_by_length(slice);
            if (ipfps.size() < 2) return -1;
            int ipfp_proton = ipfps.at(1);
            if (ipfp_proton == -1) return -1;

            return slice->reco.pfp[ipfp_proton].trackScore;
        });

        const ana::Var slice_Np ([](const caf::SRSliceProxy *slice) -> double {
            std::vector<int> ipfps = ordered_protons_by_length(slice);
            // if (ipfps.size() == 0) return -1;
            return ipfps.size();
        }); // const ana::Var slice_Np
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
    particle_data::int_type_t classification_type (const caf::SRSpillProxy *spill, const caf::SRSliceProxy *slice) {
        /* Utility to help classify different true interaction types
         *  - 1µNp are particle_data::int_type_t::true_visible_1muNp
         *  - 1µ1p are particle_data::int_type_t::true_visible_1mu1p
         *  - unclassified
        */

       
       if (std::isnan(slice->vertex.x) || std::isnan(slice->vertex.y) || std::isnan(slice->vertex.z))
       return particle_data::int_type_t::unclassified;
       
       if (std::isnan(slice->truth.position.x) || std::isnan(slice->truth.position.y) || std::isnan(slice->truth.position.z))
       return particle_data::int_type_t::unclassified;
       
       if (slice->truth.index < 0)
       return particle_data::int_type_t::unclassified;
       
       if (!(
           std::abs(slice->truth.pdg) == 14            && 
           slice->truth.iscc                           && 
           !std::isnan(slice->truth.position.x)        && 
           !std::isnan(slice->truth.position.y)        && 
           !std::isnan(slice->truth.position.z)        // &&
        //    cuts::in_FV (slice->truth.position.x, slice->truth.position.y, slice->truth.position.z) &&
        //    cuts::in_active (slice->truth.position.x, slice->truth.position.y, slice->truth.position.z)
        )) return particle_data::int_type_t::unclassified;
        
        int num_protons_above50 = 0;
        int num_muons = 0;
        int num_pions = 0;
        int num_neutral_pions = 0;
        int num_gamma = 0;
        double length_muon = 0;
        double dep_E = 0;
 
        int use_plane = 2;

        for (auto const& prim: slice->truth.prim) {

            dep_E = 0; 
            
            if (prim.G4ID < 0)
                continue;
            
            if (prim.cryostat < 0)
                continue;

            if (std::abs(prim.pdg) == 211) {
                num_pions += 1;
            } // found charged pion, returning

            /* Looking at neutral pions: trickier
             * Reject if any of daughter photon is > 25 MeV
            */ 
            if (std::abs(prim.pdg) == 111 && prim.daughters.size() > 0) {
                for (auto const& true_particle: spill->true_particles) {
                    if (
                        true_particle.parent == prim.G4ID && std::abs(true_particle.pdg) == 22 &&
                        true_particle.plane[prim.cryostat][use_plane].visE * particle_data::GeV > particle_data::minimum_gamma_MeV
                    ){
                        num_neutral_pions += 1;
                    }
                } // loop spill->true_particles
            } // neutral pion found

            if (std::abs(prim.pdg) == 13) {
                num_muons += 1;
                length_muon = prim.length;
            } // found muon
            
            if (std::abs(prim.pdg) == 22) {
                if (prim.daughters.size() > 0) {
                    for (auto const& true_particle: spill->true_particles) {
                        if (true_particle.parent == prim.G4ID)
                            dep_E += true_particle.plane[prim.cryostat][use_plane].visE * particle_data::GeV;
                    } // loop trough true_particles
                } // photon has daughters
                dep_E += prim.plane[prim.cryostat][use_plane].visE * particle_data::GeV;
            } // found photon

            // if primary is photon with dep energy > 25 MeV skip...
            if (std::abs(prim.pdg) == 22 && dep_E > particle_data::minimum_gamma_MeV)
                num_gamma += 1;

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
            
            dep_E = 0;
        } // loop slice->truth.prim

        if (num_gamma == 0 && num_pions == 0 && num_neutral_pions == 0) {
            if (
                num_muons == 1 && num_protons_above50 == 1 && 
                length_muon > 50 && cuts::all_trk_contained_truth(spill, slice)
            ) return particle_data::int_type_t::true_visible_1mu1p;
    
            if (
                num_muons == 1 && num_protons_above50 == 2 && 
                length_muon > 50 && cuts::all_trk_contained_truth(spill, slice)
            ) return particle_data::int_type_t::true_visible_1mu2p;
    
            if (
                num_muons == 1 && num_protons_above50 == 3 && 
                length_muon > 50 && cuts::all_trk_contained_truth(spill, slice)
            ) return particle_data::int_type_t::true_visible_1mu3p;
    
            if (
                num_muons == 1 && num_protons_above50 > 1 && 
                length_muon > 50 && cuts::all_trk_contained_truth(spill, slice)
            ) return particle_data::int_type_t::true_visible_1muNp;
        } else {
            // here fallback for additional tests
            // interesting topologies to search:
            //  - 1mu1pi
            if (
                num_muons == 1 && num_pions == 1 && 
                length_muon > 50 && cuts::all_trk_contained_truth(spill, slice)
            ) return particle_data::int_type_t::true_visible_1mu1pi;
            //  - 1mu2pi
            if (
                num_muons == 1 && num_pions == 2 && 
                length_muon > 50 && cuts::all_trk_contained_truth(spill, slice)
            ) return particle_data::int_type_t::true_visible_1mu2pi;
            //  - 1mu3pi
            if (
                num_muons == 1 && num_pions == 3 && 
                length_muon > 50 && cuts::all_trk_contained_truth(spill, slice)
            ) return particle_data::int_type_t::true_visible_1mu3pi;
            //  - 1muNpi
            if (
                num_muons == 1 && num_pions > 0 && 
                length_muon > 50 && cuts::all_trk_contained_truth(spill, slice)
            ) return particle_data::int_type_t::true_visible_1muNpi;
            //  - 2mu
            if (
                num_muons == 2 &&  
                length_muon > 50 && cuts::all_trk_contained_truth(spill, slice)
            ) return particle_data::int_type_t::true_visible_2mu;
            //  - 2p
            if (
                num_protons_above50 == 2 &&
                cuts::all_trk_contained_truth(spill, slice)
            ) return particle_data::int_type_t::true_visible_2p;
            //  - 1muShortNp (no muon 50cm lenght)
            if (
                num_muons == 1 &&
                cuts::all_trk_contained_truth(spill, slice)
            ) return particle_data::int_type_t::true_visible_1muShortNp;
            //  - 1pi1p
            if (
                num_protons_above50 == 1 && num_pions == 1 &&
                cuts::all_trk_contained_truth(spill, slice)
            ) return particle_data::int_type_t::true_visible_1pi1p;
            //  - 1mu1p1pi0
            if (
                num_muons == 1 && num_pions == 1 && num_neutral_pions == 1 &&
                length_muon > 50 && cuts::all_trk_contained_truth(spill, slice)
            ) return particle_data::int_type_t::true_visible_1mu1p1pi0;
            //  - 1mu1pNpi0
            if (
                num_muons == 1 && num_pions == 1 && num_neutral_pions > 1 &&
                length_muon > 50 && cuts::all_trk_contained_truth(spill, slice)
            ) return particle_data::int_type_t::true_visible_1mu1pNpi0;

        }

        return particle_data::int_type_t::unclassified;
    } // int classification_type

    particle_data::int_type_t classification_type_MC (
        const caf::SRSpillProxy *spill, const caf::Proxy<caf::SRTrueInteraction> &nu
    ) {
        /* Utility to help classify different true interaction types
         *  - 1µNp are particle_data::int_type_t::true_visible_1muNp
         *  - 1µ1p are particle_data::int_type_t::true_visible_1mu1p
         *  - all others are particle_data::int_type_t::unclassified
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
                return particle_data::int_type_t::unclassified;
            } // found charged pion, returning

            /* Looking at neutral pions: trickier
             * Reject if any of daughter photon is > 25 MeV
            */ 
            if (std::abs(prim.pdg) == 111 && prim.daughters.size() >= 1) {
                for (auto const& true_particle: spill->true_particles) {
                    G4ID_parent = true_particle.parent;

                    if (
                        G4ID_parent == prim.G4ID && std::abs(true_particle.pdg) == 22 &&
                        true_particle.plane[prim.cryostat][use_plane].visE * particle_data::GeV > particle_data::minimum_gamma_MeV
                    ){
                        return particle_data::int_type_t::unclassified;
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
                            dep_E += true_particle.plane[prim.cryostat][use_plane].visE * particle_data::GeV;
                    } // loop trough true_particles
                } // photon has daughters
                dep_E += prim.plane[prim.cryostat][use_plane].visE * particle_data::GeV;
            } // found photon

            // if primary is photon with dep energy > 25 MeV skip...
            if (std::abs(prim.pdg) == 22 && dep_E > particle_data::minimum_gamma_MeV)
                return particle_data::int_type_t::unclassified;

            if (std::abs(prim.pdg) == 2212) {
                if (prim.daughters.size() >= 1) {
                    for(auto const& true_particle: spill->true_particles) {
                        G4ID_parent = true_particle.parent;
                        if (G4ID_parent == prim.G4ID) 
                            dep_E += true_particle.plane[prim.cryostat][use_plane].visE * particle_data::GeV;
                    } // loop trough true_particles
                } // proton has daughters
                dep_E += prim.plane[prim.cryostat][use_plane].visE * particle_data::GeV;
            } // found proton

            if (std::abs(prim.pdg) == 2212 && dep_E > 50.0)
                num_protons_above50 += 1;
            
            dep_E = 0;
        } // loop slice->truth.prim

        if (
            num_muons == 1 && num_protons_above50 == 1 && 
            length_muon > 50 && cuts::all_trk_contained_MC(spill, nu)
        ) return particle_data::int_type_t::true_visible_1mu1p;

        if (
            num_muons == 1 && num_protons_above50 > 1 && 
            length_muon > 50 && cuts::all_trk_contained_MC(spill, nu)
        ) return particle_data::int_type_t::true_visible_1muNp;

        return particle_data::int_type_t::unclassified;
    } // int classification_type
} // namespace var_utils

#endif
