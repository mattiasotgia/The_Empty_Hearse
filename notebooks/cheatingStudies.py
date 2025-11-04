import numpy as np
import pandas as pd

import ultraplot as plot

import uproot
import hist
import hist.intervals
from uncertainties import ufloat

import argparse
import tomllib

plot.rc['figure.facecolor'] = 'none'
plot.rc['savefig.facecolor'] = 'none'
plot.rc['legend.frameon'] = False
plot.rc['axes.autolimit_mode'] = 'data'
plot.rc['label.size'] = 11
plot.rc['font.size'] = 11


loader_formatted = {
    'nominal': 'Nominal reconstruction',
    'cheated_2d': 'Injected clustering',
    'cheated_vtx': 'Injected vertex',
    'cheated_vtxSelection': 'Injected vertex selection',
    'cheated_3d': 'Injected 3D matching',
    'cheated_nuH': r'Injected $\nu$-hierarchy',
    'cheated_mva': 'Injected track-score',
    'cheated_2d_vtx': 'Injected clustering + vertex',
    'cheated_2d_vtx_3d': 'Injected clustering + vertex + 3D matching',
    'cheated_2d_vtx_3d_nu': r'Injected clustering + vertex + 3D matching + $\nu$-hierarchy',
    'cheated_2d_vtx_3d_nu_mva': r'Injected clustering + vertex + 3D matching + $\nu$-hierarchy + track-score'
}

class Analysis:
    def __init__(self, data, binning, tree_bases=['reco_true_{}', 'reco_{}', 'true_{}']):
        self.binning = binning
        self.data = data

        self.reco_true_base_string, self.reco_base_string, self.true_base_string = tree_bases

    def efficiency(self, tree_name: str, variable: str = 'true_E'):

        reco_true = self.data[self.reco_true_base_string.format(tree_name)].arrays(library='pd')
        true = self.data[self.true_base_string.format(tree_name)].arrays(library='pd')

        common = pd.merge(reco_true.Evt, true.Evt, on='Evt')
        
        reco_true_H = hist.Hist(self.binning).fill(reco_true[reco_true.Evt.isin(common.Evt)][variable].values)
        true_H = hist.Hist(self.binning).fill(true[true.Evt.isin(common.Evt)][variable].values)

        reco_true_values = reco_true_H.values()
        true_values = true_H.values()
        
        with np.errstate(divide="ignore", invalid="ignore"):
            efficiency = reco_true_values/true_values
            efficiency_uncertainty = hist.intervals.ratio_uncertainty(
                reco_true_values, true_values, 'efficiency'
            )
        
        return efficiency, efficiency_uncertainty
    
    def purity(self, tree_name: str, variable: str = 'true_E'):
        
        reco_true = self.data[self.reco_true_base_string.format(tree_name)].arrays(library='pd')
        reco = self.data[self.reco_base_string.format(tree_name)].arrays(library='pd')

        common = pd.merge(reco_true.Evt, reco.Evt, on='Evt')
        
        reco_true_H = hist.Hist(self.binning).fill(reco_true[reco_true.Evt.isin(common.Evt)][variable].values)
        reco_H = hist.Hist(self.binning).fill(reco[reco.Evt.isin(common.Evt)][variable].values)

        reco_true_values = reco_true_H.values()
        reco_values = reco_H.values()
        
        with np.errstate(divide="ignore", invalid="ignore"):
            purity = reco_true_values/reco_values
            purity_uncertainty = hist.intervals.ratio_uncertainty(
                reco_true_values, reco_values, 'efficiency'
            )
        
        return purity, purity_uncertainty

    def spectra(self, tree_name: str, variable: str = 'true_E'):
        
        reco_true = self.data[self.reco_true_base_string.format(tree_name)].arrays(library='pd')
        reco = self.data[self.reco_base_string.format(tree_name)].arrays(library='pd')
        true = self.data[self.true_base_string.format(tree_name)].arrays(library='pd')

        common = pd.merge(reco_true.Evt, reco.Evt, on='Evt')
        common = pd.merge(common.Evt, true.Evt, on='Evt')

        reco_true_H = hist.Hist(self.binning).fill(reco_true[reco_true.Evt.isin(common.Evt)][variable].values)
        reco_H = hist.Hist(self.binning).fill(reco[reco.Evt.isin(common.Evt)][variable].values)
        true_H = hist.Hist(self.binning).fill(true[true.Evt.isin(common.Evt)][variable].values)

        return reco_true_H, reco_H, true_H



def main():


    cliapp = argparse.ArgumentParser('cheatingStudies.py', description='Study efficiency of reconstruction for multiple algorithms')

    cliapp.add_argument('-c', '--configuration', action='store')

    args = cliapp.parse_args()

    configuration: dict = tomllib.load(open(args.configuration, 'rb'))
    globalConfiguration: dict = configuration.get('global')

    

    CCNp_data = uproot.open(globalConfiguration.get('filePath', 'CCNp_efficiencyNoNu_cheatingSlice_cheatingPid.root'))

    numberOfEnergyBins = len(configuration['energies'].keys())

    fig, ax = plot.subplots(
        width=globalConfiguration.get('figWidth', 7.5), height=globalConfiguration.get('figHeight', 6.5), ncols=numberOfEnergyBins, nrows=1, ylabel=f'Efficiency', sharey='labs',
        titlecolor='k', share=False, grid=False, 
        ltitle=r'ICARUS MC BNB $\nu$-only / true $\nu_\mu$CC QE $1\mu N\mathrm{p}$ selection',
        xticks=[0, 1, 2, 3],
        xticklabels=[
            'Injecting clustering',
            '+ vertex',
            '+ 2D $\\to$ 3D',
            '+ track-score'
        ],
    )

    cheated_dict = {
        'xerr': 0.5, 
        'markersize': 4, 
        'mec': 'w', 
        'markeredgewidth': 0.5, 
        'capsize': 0, 
        'elinewidth': 3
    }

    trees = [
        'cheated_2d',
        'cheated_2d_vt',
        'cheated_2d_vtx_3d',
        'cheated_2d_vtx_3d_mva',
    ]

    for idx, key in enumerate(configuration['energies']):

        energyConfiguration = configuration['energies'][key]

        print(f'Running over energetic space {key}')

        CCNp_analysis = Analysis(CCNp_data, hist.axis.Regular(1, energyConfiguration.get('low', 0.24), energyConfiguration.get('high', 2.4)))
        reconstruction_efficiency = {}
        ## NOMINAL RECO
        efficiency, (low, high) = CCNp_analysis.efficiency('nominal', 'E')
        ax[idx].axhspan(ymin=(efficiency-low)[0], ymax=(efficiency+high)[0], hatch='////', fill=True, color='gray4', ec='gray4', alpha=0.25, linewidth=1, label='Nominal reco.', zorder=-99)
        ax[idx].axhline(efficiency, lw=3, c='gray4')
        print(f'{'nominal'}: {efficiency[0]:.2%}')

        reconstruction_efficiency['nominal'] = (efficiency[0], np.max([low, high]))

        ## Cheated PID
        efficiency, (low, high) = CCNp_analysis.efficiency('nominal_cheatedPid', 'E')
        ax[idx].axhspan(ymin=(efficiency-low)[0], ymax=(efficiency+high)[0], hatch='////', fill=True, color='green4', ec='green4', alpha=0.25, linewidth=1, label='Nominal reco. (modified Pid)', zorder=-99)
        ax[idx].axhline(efficiency, lw=3, c='green4')
        print(f'{'nominal_cheatedPid'}: {efficiency[0]:.2%}')

        reconstruction_efficiency['nominal_cheatedPid'] = (efficiency[0], np.max([low, high]))

        ## LADDER 
        for i, tree in enumerate(trees):
            efficiency, uncertainty = CCNp_analysis.efficiency(tree, 'E')
            ax[idx].errorbar(i, efficiency, yerr=uncertainty, **cheated_dict, color='gray8', fmt='^', label=('Cheated reco.' if i==0 else None), zorder=99)
            print(f'{tree}: {efficiency[0]:.2%}')
            reconstruction_efficiency[tree] = (efficiency[0], np.max(uncertainty))

        ## LADDER cheated pid
        for i, tree in enumerate(trees):
            tree = tree + '_cheatedPid'
            efficiency, uncertainty = CCNp_analysis.efficiency(tree, 'E')
            ax[idx].errorbar(i, efficiency, yerr=uncertainty, **cheated_dict, color='green8', fmt='v', label=('Cheated reco. (modified Pid)' if i==0 else None), zorder=99)
            print(f'{tree}: {efficiency[0]:.2%}')
            reconstruction_efficiency[tree] = (efficiency[0], np.max(uncertainty))

        ax[idx].format(ylim=energyConfiguration.get('ylims', [0.524, 0.935]), ultitle=energyConfiguration.get('label', ''))
        ax[idx].legend(loc='ll', ncols=globalConfiguration.get('legendNCols', 2), order='F')

        # first numerator, second denominator (pid efficiency is inverted
        stageToConfigs = {
            '2d': ('cheated_2d', 'nominal'),
            'vtx': ('cheated_2d_vt', 'cheated_2d'),
            '3d': ('cheated_2d_vtx_3d', 'cheated_2d_vt'),
            'mva': ('cheated_2d_vtx_3d_mva', 'cheated_2d_vtx_3d'),
        }

        for stage in stageToConfigs:
            A, B = stageToConfigs[stage]
            effA, effErrA = reconstruction_efficiency[A]
            effAModifiedPid, effErrAModifiedPid = reconstruction_efficiency[f'{A}_cheatedPid']

            effB, effErrB = reconstruction_efficiency[B]
            effBModifiedPid, effErrBModifiedPid = reconstruction_efficiency[f'{B}_cheatedPid']

            print('\n__________')
            print(f'Computing stage {stage}\n----------')
            print(f'Double checks')
            print(f'To get {stage} it is eff_selection({B})/eff_selection({A}) * eff_pid({A})/eff_pid({B})')
            print(f'To get eff_pid({A}) it is eff_selection({A})/eff_selection({A}_cheatedPid)')
            print(f'To get eff_pid({B}) it is eff_selection({B})/eff_selection({B}_cheatedPid)')
            print(f'   where')
            print(f'      {B} = {ufloat(effB, effErrB)*100:uS}%')
            print(f'      {A} = {ufloat(effA, effErrA)*100:uS}%')
            print(f'   where')
            print(f'      {B}_cheatedPid = {ufloat(effBModifiedPid, effErrBModifiedPid)*100:uS}%')
            print(f'      {A}_cheatedPid = {ufloat(effAModifiedPid, effErrAModifiedPid)*100:uS}%')

            effPidA = effA/effAModifiedPid
            effErrPidA = np.sqrt((effErrA / effA)**2 + (effErrAModifiedPid / effAModifiedPid)**2) * effPidA
            
            effPidB = effB/effBModifiedPid
            effErrPidB = np.sqrt((effErrB / effB)**2 + (effErrBModifiedPid / effBModifiedPid)**2) * effPidB

            print('__________')
            print(f'eff_pid({A}) = {ufloat(effPidA, effErrPidA)*100:uS}%')
            print(f'eff_pid({B}) = {ufloat(effPidB, effErrPidB)*100:uS}%')

            effStage = ( effB / effA ) * ( effPidA / effPidB )
            effErrStage = np.sqrt((effErrA / effA)**2 + (effErrB / effB)**2 + (effErrPidA / effPidA)**2 + (effErrPidB / effPidB)**2) * effStage

            print('__________')
            print(f'eff({stage}) = {ufloat(effStage, effErrStage)*100:uS}%')


    fig.savefig(configuration['global'].get('pdfPlotPath', 'plots_CCNp_test/CCNp_efficiencyNoNu_cheatingPid.pdf'), bbox_inches='tight')


if __name__ == '__main__':
    main()