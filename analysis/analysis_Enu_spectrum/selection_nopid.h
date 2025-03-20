#ifndef SELECTION_NOPID_H
#define SELECTION_NOPID_H

#include "helper.h"
#include "selection.h"

#include <functional>
#include <algorithm>

namespace cuts {
    namespace truth {
        using int_type_t = particle_data::int_type;
        const ana::SpillCut spill_1muNp_MC ([](const caf::SRSpillProxy *spill) -> bool {
            for (auto const& nu: spill->mc.nu) {
                if (var_utils::classification_type_MC(spill, nu) == int_type_t::true_visible_1muNp) 
                    return true;
            } // loop spill->mc.nu
            return false;
        }); // const ana::SpillCut spill_1muNp_MC

        const ana::SpillCut spill_1mu1p_MC ([](const caf::SRSpillProxy *spill) -> bool {
            for (auto const& nu: spill->mc.nu) {
                if (var_utils::classification_type_MC(spill, nu) == int_type_t::true_visible_1mu1p) 
                    return true;
            } // loop spill->mc.nu
            return false;
        }); // const ana::SpillCut spill_1mu1p_MC

        const ana::SpillCut spill_sliceless ([](const caf::SRSpillProxy *spill) -> bool {
            return spill->nslc != 0;
        });
    } // namespace truth
} // namespace cuts

namespace no_pid {

    using int_type_t = particle_data::int_type;
    template<class T, class A>
    const T make_spill_from_slice (
        std::function<A(const caf::SRSliceProxy*, const int&)> nopid_slice_var, 
        const ana::Cut &truth_cut, 
        int_type_t int_type = int_type_t::true_visible_1muNp
    ) {
        return T([=](const caf::SRSpillProxy *spill) -> A {
            
            // std::cout << "This is run.event = " << run(spill) << "." << event(spill) << std::endl;
            int N_prim_true = -1;
            for (auto const& nu: spill->mc.nu) {
                if (var_utils::classification_type_MC(spill, nu) == int_type) {
                    N_prim_true = nu.prim.size();
                    break;
                } // found 1µNp/1µ1p
            } // loop spill->mc.nu

            double slice_value = -9999, longest_track = -1;
            int slice_with_longest_track = -1;
            for (std::size_t islice=0; islice<spill->nslc; ++islice) {
                if (!truth_cut(&(spill->slc.at(islice)))) continue;
                double tmp_length = -1;
                for (auto const& pfp: spill->slc.at(islice).reco.pfp) {
                    if (pfp.trk.len > tmp_length) 
                        tmp_length = pfp.trk.len;
                }
                if (tmp_length > longest_track) {
                    longest_track = tmp_length;
                    slice_with_longest_track = islice;
                }
            } // loop spill->nslc

            if (slice_with_longest_track < 0 || N_prim_true < 0) 
                return nopid_slice_var(nullptr, N_prim_true);

            return nopid_slice_var(&(spill->slc[slice_with_longest_track]), N_prim_true);
        });
    } // const ana::SpillVar make_spill_from_slice


    int longest_pfp (const caf::SRSliceProxy *slice) {
        double longest_track = -1;
        int longest_ipfp = -1;
        for (std::size_t ipfp=0; ipfp<slice->reco.npfp; ++ipfp) {
            if (slice->reco.pfp[ipfp].trk.len > longest_track) {
                longest_track = slice->reco.pfp[ipfp].trk.len;
                longest_ipfp = ipfp;
            }
        } // loop ipfp<slice->npfp
        return longest_ipfp;
    } // int longest_pfp

    std::vector<int> ordered_pfps_by_length (const caf::SRSliceProxy *slice, const int& N) {
        std::vector<std::pair<int, double>> ipfp_length_pairs_all;
        std::vector<int> ordered_ipfp;

        for (std::size_t ipfp=0; ipfp<slice->reco.npfp; ++ipfp) {
            if (!std::isnan(slice->reco.pfp[ipfp].trk.len))
                ipfp_length_pairs_all.emplace_back(ipfp, slice->reco.pfp[ipfp].trk.len);
                // safety check
        } // loop ipfp<slice->npfp

        std::sort(ipfp_length_pairs_all.begin(), ipfp_length_pairs_all.end(), 
            [](const std::pair<int, double> &a, const std::pair<int, double> &b) { return a.second > b.second; }
        );

        for (auto const& ipfp_length: ipfp_length_pairs_all) 
            ordered_ipfp.push_back(ipfp_length.first);

        if (slice->reco.npfp < N || ordered_ipfp.size() < N) 
            return ordered_ipfp;
        
        std::vector<int> tmp_ipfp;
        for (std::size_t i=0; i<ordered_ipfp.size(); ++i)
            tmp_ipfp.push_back(ordered_ipfp[i]);

        return tmp_ipfp;
    } // std::vector<int> ordered_pfps_by_length

    namespace vars {
        double slice_muon_hit_completeness (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::numeric_limits<double>::min();
            int ipfp_muon = longest_pfp(slice);
            if (ipfp_muon < 0) return std::numeric_limits<double>::min();

            return slice->reco.pfp.at(ipfp_muon).trk.truth.bestmatch.hit_completeness;
        } // double slice_muon_hit_completeness

        double slice_leading_proton_hit_completeness (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::numeric_limits<double>::min();
            std::vector<int> ipfps = ordered_pfps_by_length(slice, N);
            int ipfp_proton = -1;
            if (ipfps.size() > 1) 
                ipfp_proton = ipfps.at(1);
            else 
                return std::numeric_limits<double>::min();
            
            return slice->reco.pfp.at(ipfp_proton).trk.truth.bestmatch.hit_completeness;
        } // double slice_leading_proton_hit_completeness
        
        std::vector<double> slice_protons_hit_completeness (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::vector<double>();
            std::vector<int> ipfps = ordered_pfps_by_length(slice, N);
            if (ipfps.size() < 1) 
                return std::vector<double>();

            std::vector<double> values;
            for (std::size_t i=1; i<ipfps.size(); i++) 
                values.push_back(slice->reco.pfp.at(ipfps.at(i)).trk.truth.bestmatch.hit_completeness);

            return values;
        } // std::vector<double> slice_protons_hit_completeness
        
        double slice_muon_hit_purity (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::numeric_limits<double>::min();
            int ipfp_muon = longest_pfp(slice);
            if (ipfp_muon < 0) return std::numeric_limits<double>::min();

            return slice->reco.pfp.at(ipfp_muon).trk.truth.bestmatch.hit_purity;
        } // double slice_muon_hit_purity

        double slice_leading_proton_hit_purity (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::numeric_limits<double>::min();
            std::vector<int> ipfps = ordered_pfps_by_length(slice, N);
            int ipfp_proton = -1;
            if (ipfps.size() > 1) 
                ipfp_proton = ipfps.at(1);
            else 
                return std::numeric_limits<double>::min();

            return slice->reco.pfp.at(ipfp_proton).trk.truth.bestmatch.hit_purity;
        } // double slice_leading_proton_hit_purity

        std::vector<double> slice_protons_hit_purity (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::vector<double>();
            std::vector<int> ipfps = ordered_pfps_by_length(slice, N);
            if (ipfps.size() < 1) 
                return std::vector<double>();

            std::vector<double> values;
            for (std::size_t i=1; i<ipfps.size(); i++) 
                values.push_back(slice->reco.pfp.at(ipfps.at(i)).trk.truth.bestmatch.hit_purity);

            return values;
        } // std::vector<double> slice_protons_hit_purity

        double slice_muon_L_reco_true_ratio (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::numeric_limits<double>::min();
            int ipfp_muon = longest_pfp(slice);
            if (ipfp_muon < 0) return std::numeric_limits<double>::min();

            return slice->reco.pfp.at(ipfp_muon).trk.len / slice->reco.pfp.at(ipfp_muon).trk.truth.p.length;
        } // double slice_muon_L_reco_true_ratio

        double slice_leading_proton_L_reco_true_ratio (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::numeric_limits<double>::min();
            std::vector<int> ipfps = ordered_pfps_by_length(slice, N);
            int ipfp_proton = -1;
            if (ipfps.size() > 1) 
                ipfp_proton = ipfps.at(1);
            else 
                return std::numeric_limits<double>::min();

            return slice->reco.pfp.at(ipfp_proton).trk.len / slice->reco.pfp.at(ipfp_proton).trk.truth.p.length;
        } // double slice_leading_proton_L_reco_true_ratio

        std::vector<double> slice_protons_L_reco_true_ratio (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::vector<double>();
            std::vector<int> ipfps = ordered_pfps_by_length(slice, N);
            if (ipfps.size() < 1) 
                return std::vector<double>();

            std::vector<double> values;
            for (std::size_t i=1; i<ipfps.size(); i++) 
                values.push_back(
                    slice->reco.pfp.at(ipfps.at(i)).trk.len / slice->reco.pfp.at(ipfps.at(i)).trk.truth.p.length
                );

            return values;
        } // std::vector<double> slice_protons_L_reco_true_ratio

        double slice_vertex_delta_3D (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::numeric_limits<double>::min();
            if (
                std::isnan(slice->vertex.x) ||
                std::isnan(slice->vertex.y) ||
                std::isnan(slice->vertex.z) ||
                std::isnan(slice->truth.position.x) ||
                std::isnan(slice->truth.position.y) ||
                std::isnan(slice->truth.position.z) 
            ) return std::numeric_limits<double>::min(); 

            return std::sqrt(
                std::pow(slice->vertex.x - slice->truth.position.x, 2) + 
                std::pow(slice->vertex.y - slice->truth.position.y, 2) + 
                std::pow(slice->vertex.z - slice->truth.position.z, 2)
            );
        } // double slice_vertex_delta_3D

        double slice_vertex_delta_x (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::numeric_limits<double>::min();
            if (
                std::isnan(slice->vertex.x) ||
                std::isnan(slice->truth.position.x) 
            ) return std::numeric_limits<double>::min(); 

            return std::sqrt(
                std::pow(slice->vertex.x - slice->truth.position.x, 2)
            );
        } // double slice_vertex_delta_x

        double slice_vertex_delta_y (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::numeric_limits<double>::min();
            if (
                std::isnan(slice->vertex.y) ||
                std::isnan(slice->truth.position.y) 
            ) return std::numeric_limits<double>::min(); 

            return std::sqrt(
                std::pow(slice->vertex.y - slice->truth.position.y, 2)
            );
        } // double slice_vertex_delta_y

        double slice_vertex_delta_z (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::numeric_limits<double>::min();
            if (
                std::isnan(slice->vertex.z) ||
                std::isnan(slice->truth.position.z) 
            ) return std::numeric_limits<double>::min(); 

            return std::sqrt(
                std::pow(slice->vertex.z - slice->truth.position.z, 2)
            );
        } // double slice_vertex_delta_z

        double slice_muon_P_reco_true_ratio (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::numeric_limits<double>::min();
            int ipfp_muon = longest_pfp(slice);
            if (ipfp_muon < 0) return std::numeric_limits<double>::min();

            double startp_mag = std::sqrt(
                std::pow(slice->reco.pfp.at(ipfp_muon).trk.truth.p.startp.x, 2) + 
                std::pow(slice->reco.pfp.at(ipfp_muon).trk.truth.p.startp.y, 2) + 
                std::pow(slice->reco.pfp.at(ipfp_muon).trk.truth.p.startp.z, 2)
            );

            return slice->reco.pfp.at(ipfp_muon).trk.rangeP.p_muon / startp_mag;
        } // const ana::Var slice_muon_L_reco_true_ratio

        double slice_leading_proton_P_reco_true_ratio (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::numeric_limits<double>::min();
            std::vector<int> ipfps = ordered_pfps_by_length(slice, N);
            int ipfp_proton = -1;
            if (ipfps.size() > 1) 
                ipfp_proton = ipfps.at(1);
            else 
                return std::numeric_limits<double>::min();

            double startp_mag = std::sqrt(
                std::pow(slice->reco.pfp.at(ipfp_proton).trk.truth.p.startp.x, 2) + 
                std::pow(slice->reco.pfp.at(ipfp_proton).trk.truth.p.startp.y, 2) + 
                std::pow(slice->reco.pfp.at(ipfp_proton).trk.truth.p.startp.z, 2)
            );

            return slice->reco.pfp.at(ipfp_proton).trk.rangeP.p_proton / startp_mag;
        } // const ana::Var slice_leading_proton_L_reco_true_ratio

        std::vector<double> slice_protons_P_reco_true_ratio (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::vector<double>();
            std::vector<int> ipfps = ordered_pfps_by_length(slice, N);
            if (ipfps.size() < 1) 
                return std::vector<double>();

            std::vector<double> values;
            double startp_mag = 0;
            for (std::size_t i=1; i<ipfps.size(); i++) {
                
                startp_mag = std::sqrt(
                    std::pow(slice->reco.pfp.at(ipfps.at(i)).trk.truth.p.startp.x, 2) + 
                    std::pow(slice->reco.pfp.at(ipfps.at(i)).trk.truth.p.startp.y, 2) + 
                    std::pow(slice->reco.pfp.at(ipfps.at(i)).trk.truth.p.startp.z, 2)
                );

                values.push_back(
                    slice->reco.pfp.at(ipfps.at(i)).trk.rangeP.p_proton / startp_mag
                );
            } // loop

            return values;
        } // const ana::Var slice_protons_L_reco_true_ratio

        double slice_neutrino_reco_E (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::numeric_limits<double>::min();
            int ipfp_muon = longest_pfp(slice);
            if (ipfp_muon < 0) return std::numeric_limits<double>::min();

            double p_mu_x = slice->reco.pfp[ipfp_muon].trk.rangeP.p_muon * slice->reco.pfp[ipfp_muon].trk.dir.x;   // Momenta are in GeV
            double p_mu_y = slice->reco.pfp[ipfp_muon].trk.rangeP.p_muon * slice->reco.pfp[ipfp_muon].trk.dir.y;   // Momenta are in GeV
            double p_mu_z = slice->reco.pfp[ipfp_muon].trk.rangeP.p_muon * slice->reco.pfp[ipfp_muon].trk.dir.z;   // Momenta are in GeV

            double p_mu_tot = std::sqrt(p_mu_x * p_mu_x + p_mu_y * p_mu_y + p_mu_z * p_mu_z); // [GeV]

            double E_mu = 1000 * std::sqrt(p_mu_tot * p_mu_tot + std::pow(particle_data::masses::muon, 2) / (1000 * 1000));
            double E_p = 0;

            std::vector<int> Np_ipfp = ordered_pfps_by_length(slice, N);

            for (auto const& ipfp: Np_ipfp) {
                if (int(ipfp) == ipfp_muon)
                    continue;
                
                TVector3 start_mom;
                start_mom.SetXYZ(
                    slice->reco.pfp[ipfp].trk.rangeP.p_proton * slice->reco.pfp[ipfp].trk.dir.x, 
                    slice->reco.pfp[ipfp].trk.rangeP.p_proton * slice->reco.pfp[ipfp].trk.dir.y, 
                    slice->reco.pfp[ipfp].trk.rangeP.p_proton * slice->reco.pfp[ipfp].trk.dir.z
                );
                E_p += std::sqrt(
                    std::pow(particle_data::masses::proton, 2) + 
                    std::pow(start_mom.Mag() * 1000, 2)
                ) - particle_data::masses::proton;
            } // loop Np_ipfp

            return (E_mu + E_p) / 1000;
        } // double slice_neutrino_reco_E

        double slice_neutrino_true_E (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::numeric_limits<double>::min();
            if (std::isnan(slice->truth.E)) return std::numeric_limits<double>::min();

            return  slice->truth.E;
        } // double slice_neutrino_true_E

        double slice_neutrino_reco_dE (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::numeric_limits<double>::min();
            int ipfp_muon = longest_pfp(slice);
            if (ipfp_muon < 0) return std::numeric_limits<double>::min();

            double p_mu_x = slice->reco.pfp[ipfp_muon].trk.rangeP.p_muon * slice->reco.pfp[ipfp_muon].trk.dir.x;   // Momenta are in GeV
            double p_mu_y = slice->reco.pfp[ipfp_muon].trk.rangeP.p_muon * slice->reco.pfp[ipfp_muon].trk.dir.y;   // Momenta are in GeV
            double p_mu_z = slice->reco.pfp[ipfp_muon].trk.rangeP.p_muon * slice->reco.pfp[ipfp_muon].trk.dir.z;   // Momenta are in GeV

            double p_mu_tot = std::sqrt(p_mu_x * p_mu_x + p_mu_y * p_mu_y + p_mu_z * p_mu_z); // [GeV]

            double E_mu = 1000 * std::sqrt(p_mu_tot * p_mu_tot + std::pow(particle_data::masses::muon, 2) / (1000 * 1000));
            double E_p = 0;

            std::vector<int> Np_ipfp = ordered_pfps_by_length(slice, N);

            for (auto const& ipfp: Np_ipfp) {
                if (int(ipfp) == ipfp_muon)
                    continue;
                
                TVector3 start_mom;
                start_mom.SetXYZ(
                    slice->reco.pfp[ipfp].trk.rangeP.p_proton * slice->reco.pfp[ipfp].trk.dir.x, 
                    slice->reco.pfp[ipfp].trk.rangeP.p_proton * slice->reco.pfp[ipfp].trk.dir.y, 
                    slice->reco.pfp[ipfp].trk.rangeP.p_proton * slice->reco.pfp[ipfp].trk.dir.z
                );
                E_p += std::sqrt(
                    std::pow(particle_data::masses::proton, 2) + 
                    std::pow(start_mom.Mag() * 1000, 2)
                ) - particle_data::masses::proton;
            } // loop Np_ipfp
            
            double reco_E = (E_mu + E_p) / 1000;

            return slice->truth.E - reco_E/slice->truth.E;
        } // double slice_neutrino_reco_dE

        double slice_neutrino_reco_pT (const caf::SRSliceProxy *slice, const int& N) {
            if (slice == nullptr) return std::numeric_limits<double>::min();
            int ipfp_muon = longest_pfp(slice);
            if (ipfp_muon < 0) return std::numeric_limits<double>::min();

            double p_mu_x = slice->reco.pfp[ipfp_muon].trk.rangeP.p_muon * slice->reco.pfp[ipfp_muon].trk.dir.x;   // Momenta are in GeV
            double p_mu_y = slice->reco.pfp[ipfp_muon].trk.rangeP.p_muon * slice->reco.pfp[ipfp_muon].trk.dir.y;   // Momenta are in GeV
            double p_p_x = 0, p_p_y = 0;

            std::vector<int> Np_ipfp = ordered_pfps_by_length(slice, N);

            for (auto const& ipfp: Np_ipfp) {
                if (int(ipfp) == ipfp_muon)
                    continue;

                    p_p_x += (slice->reco.pfp[ipfp].trk.rangeP.p_proton) * slice->reco.pfp[ipfp].trk.dir.x;
                    p_p_y += (slice->reco.pfp[ipfp].trk.rangeP.p_proton) * slice->reco.pfp[ipfp].trk.dir.y;
            } // loop Np_ipfp
            return std::sqrt(std::pow(p_mu_x + p_p_x, 2) + std::pow(p_mu_y + p_p_y, 2));
        } // double slice_neutrino_reco_pT
    } // namespace vars
} // namespace no_pid



#endif // SELECTION_NOPID_H