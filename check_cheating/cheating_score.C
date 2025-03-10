
#ifndef CHEATING_SCORE_C
#define CHEATING_SCORE_C

#include "helper.h"

#include "TCanvas.h"
#include "TH1D.h"

#include "TStyle.h"
#include "TLegend.h"

template <class T> 
struct plot1D {
    std::string name;
    const ana::Binning bin;
    T var;
};

struct spectra {
    ana::Spectrum *spectra;
    std::string name;
};

void cheating_score() {

    gStyle->SetOptStat(0);

    // REMARK: not all events are in both processing, some are only in the 
    // not cheated (17700 evt), nd not in the cheated sample (17500)...
    // These are accounted for in the hard_code namespace (hard_code::is_bad_event...)
    // 

    // ana::SpectrumLoader loader10evt_noSAM_fullcheat("/exp/icarus/data/users/msotgia/thesis/stage1_runs_single/50evt_tests/test_4_jobsub/stage1_full_cheat_essential_50.flat.caf.root");
    // ana::SpectrumLoader loader10evt_noSAM_nocheat("/exp/icarus/data/users/msotgia/thesis/stage1_runs_single/50evt_tests/test_4_jobsub/stage1_no_cheat_essential_50.flat.caf.root");
    // ana::SpectrumLoader dataloader_cheated("msotgia_v09_89_01_01p03_BNB_production_cheated_stage1tocaf_flatcafs");
    // ana::SpectrumLoader dataloader_non_cheated("msotgia_v09_89_01_01p03_BNB_production_non_cheated_stage1tocaf_flatcafs");
    ana::SpectrumLoader loader2evt_noSAM_fullcheat("/exp/icarus/data/users/msotgia/thesis/stage1_runs_single/single_1muNp_stages1tocafs/stage1_cheat_essential.flat.caf.root");
    ana::SpectrumLoader loader2evt_noSAM_nocheat("/exp/icarus/data/users/msotgia/thesis/stage1_runs_single/single_1muNp_stages1tocafs/stage1_no_cheat_essential.flat.caf.root");

    
    std::vector<plot1D<ana::MultiVar>> plots = {
        {"cheating_checks", bins::simple, cheating::score}, 
    };

    logger::log(logger::level::info) << "Starting RUN CAFAna..." << std::endl;

    std::vector<spectra> cheating_spectra;

        // /// The only \ref MultiVar variant available
        // Spectrum(const std::string& label, const Binning& bins,
        //     SpectrumLoaderBase& loader,
        //     const MultiVar& var,
        //     const SpillCut& spillcut,
        //     const Cut& cut,
        //     const SystShifts& shift = kNoShift,
        //     const Var& wei = kUnweighted);

    for (auto const& plot: plots) {
        cheating_spectra.push_back({
            new ana::Spectrum(plot.name, plot.bin, loader2evt_noSAM_fullcheat, plot.var, cheating::cut_bad_events, cheating::cut_numuCC), plot.name
        });
        cheating_spectra.push_back({
            new ana::Spectrum(plot.name, plot.bin, loader2evt_noSAM_nocheat, plot.var, cheating::cut_bad_events, cheating::cut_numuCC), plot.name
        });
    }

    loader2evt_noSAM_fullcheat.Go();
    loader2evt_noSAM_nocheat.Go();

    auto figure = new TCanvas("", "", 600, 500);

    auto cheated = cheating_spectra.at(0).spectra->ToTH1(cheating_spectra.at(0).spectra->POT());
    cheated->SetLineColor(kRed);
    cheated->SetFillColorAlpha(kRed, 0.5);
    cheated->SetLineWidth(5);
    cheated->SetFillStyle(3004);
    cheated->SetTitle(";#leftarrow Cheating score [arb. u.]; Entries / bin");

    auto non_cheated = cheating_spectra.at(1).spectra->ToTH1(cheating_spectra.at(1).spectra->POT());
    non_cheated->SetLineColor(kBlack);
    non_cheated->SetFillColorAlpha(kBlack, 0.5);
    non_cheated->SetFillStyle(3005);

    cheated->Draw("HIST");
    non_cheated->Draw("HIST SAME");

    auto legend = new TLegend(0.45,0.65,0.85,0.85);
    legend->AddEntry(cheated, "Cheated events", "F");
    legend->AddEntry(non_cheated, "Non cheated events", "F");
    legend->SetLineWidth(0);
    legend->Draw();

    // gPad->SetLogy();
    // gPad->SetLogx();

    figure->Print("plot.pdf");

    return;
}

#endif // CHEATING_SCORE_C
    