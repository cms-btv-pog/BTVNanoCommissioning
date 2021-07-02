import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import mplhep as hep
from coffea.util import load
from coffea.hist import plot
import coffea.hist as hist
import re

#output = load('hists_ctag_ctag_AK4_correctID.coffea')
# output = load('hists_ctag_ctag_AK4_allcut.coffea')
#output = load('hists_dilep_dilep_AK4_tightpu.coffea')
output = load('hists_ctag_ctag_AK4_rawpT.coffea')
#output = load('hists_ctag_ctag_AK4_drmet.coffea')
plt.style.use([hep.style.ROOT, {'font.size': 16}])

# plot options for data
data_err_opts = {
    'linestyle': 'none',
    'marker': '.',
    'markersize': 10.,
    'color': 'k',
    'elinewidth': 1,
}


from cycler import cycler
colors=["#F44336","#E91E63","#9C27B0","#673AB7","#3F51B5","#2196F3","#03A9F4","#00BCD4","#009688","#4CAF50","#8BC34A","#CDDC39","#FFEB3B","#FFC107","#FF9800","#FF5722","#795548","#BDBDBD","#9E9E9E","#616161","#90BED4","#607D8B","#455A64"]
print(output.keys())
# h1name = list(output.keys())[len(output.keys())-2]
# h1 = output[h1name]
# plot.plot1d(h1)
# for i in output['DeepCSV_vertexEnergyRatio'],output['DeepCSV_trackPtRel_1']:
for i in range(1,2):
#for i in range(len(output.keys())-8,  len(output.keys())-5):
# for i in range( len(output.keys())-5, len(output.keys())-1):
# regions=['lumimask','trigger', 'musel', 'jetsel','metsel']
#, 'jetdr', 'jetpuid', 'jetjetid','metsel']
# for r in regions:
# for i in range(0,len(output.keys())-1):
        fig, ((ax1),(rax1)) = plt.subplots(2, 1, figsize=(8, 6), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
        fig.subplots_adjust(hspace=.07)
        h1name = list(output.keys())[i]
        # h1name='DeepCSV_trackPtRel_1'
        if any([h1name.startswith('cutflow')]): 
            break
    #     if not any([h1name.startswith('btagCSVV2')]): 
    #         coffe
        h1 = output[h1name]
        dense = False 
        for ax, rax, h in zip([ax1], [rax1], [h1]):
            # ax.set_prop_cycle(cycler('color', colors))

            ###runD
            scales={
#                 'ZZ_TuneCP5_13TeV-pythia8': 3.7054E-02, 
#                 'WZ_TuneCP5_13TeV-pythia8': 5.3457E-02, 
#                 'WW_TuneCP5_13TeV-pythia8':3.5154E-02,
                'QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.3752E+04,
                'QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8':4.5788E+02, 
                'QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8':2.4585E+02, 
                'QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8':7.9470E+01,
                'QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.9696E+01,
                'QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8':5.1390E+00,
                'QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.0724E+01,
                'QCD_Pt-300to470_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.9869E-01,
                'QCD_Pt-470to600_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.4074E-02,
                'QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8':6.2639E-03,
                'QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.1944E-03,
                'QCD_Pt-1000toinf_MuEnrichedPt5_TuneCP5_13TeV_pythia8':5.9589E-04,
                'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8':3.5635E+00,
                'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8':2.6664E-01,        
                'ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8':5.6511E-04, 
                'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8':5.5109E-04, 
                'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8':3.8952E-03, 
                'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8':4.8001E-04,
                'TTToHadronic_TuneCP5_13TeV-powheg-pythia8':3.2535E-02,
                'ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8':9.8205E-02, 
                'ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8':9.5112E-02,
                'ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8':4.4811E-04,
                'ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8':4.3698E-04, 
             }
            ## Small set runD
            # scales={
            #     'QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.35E+04,
            #     'QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.61E+04, 
            #     'QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8':2.14E+03, 
            #     'QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.94E+03,
            #     'QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8':4.66E+02,
            #     'QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.31E+02,
            #     'QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8':3.30E+01,
            #     'QCD_Pt-300to470_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.95E+01,
            #     'QCD_Pt-470to600_MuEnrichedPt5_TuneCP5_13TeV_pythia8':3.36E-01,
            #     'QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8':2.64E-01,
            #     'QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.13E-01,
            #     'QCD_Pt-1000toinf_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.31E-02,
            #     'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8':4.087E+02,
            #     'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8':6.048E+01,
            #     'ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8':8.029E-02, 
            #     'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8':4.904E-02, 
            #     'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8':6.684E-02, 
            #     'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8':1.901E-02,
            #     'TTToHadronic_TuneCP5_13TeV-powheg-pythia8':2.780E+00, 
            #     'ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8':2.511E+00  , 
            #     'ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8':1.984E+00,
            #     'ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8':1.178E-01,
            #     'ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8':5.146E-02 
            # }
            ## Medium, set runD
            # scales={
            #         'QCD_Pt-15to20_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.35E+04,
            #      'QCD_Pt-20to30_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.61E+04, 
            #      'QCD_Pt-30to50_MuEnrichedPt5_TuneCP5_13TeV_pythia8':2.14E+03, 
            #      'QCD_Pt-50to80_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.94E+03,
            #      'QCD_Pt-80to120_MuEnrichedPt5_TuneCP5_13TeV_pythia8':4.66E+02,
            #      'QCD_Pt-120to170_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.31E+02,
            #      'QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8':3.30E+01,
            #      'QCD_Pt-300to470_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.95E+01,
            #      'QCD_Pt-470to600_MuEnrichedPt5_TuneCP5_13TeV_pythia8':3.36E-01,
            #      'QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.04E-02,
            #      'QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.13E-01,
            #      'QCD_Pt-1000toinf_MuEnrichedPt5_TuneCP5_13TeV_pythia8':1.31E-02,
            #      'WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8':3.508E+00,
            #      'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8':6.048E+01,
            #      'ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8':8.029E-02, 
            #      'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8':4.904E-02, 
            #      'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8':3.835E-03, 
            #      'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8':9.273E-04,
            #      'TTToHadronic_TuneCP5_13TeV-powheg-pythia8':4.970E-02,
            #      'ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8':2.511E+00  , 
            #      'ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8':1.984E+00,
            #      'ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8':1.178E-01,
            #      'ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8':5.146E-02 
            #   }
            h.scale(scales,axis='dataset')        
    #         notdata = re.compile('(?!SingleMuon)')
            notdata = re.compile('(?!Data)')
            print(h)

            if hasattr(h, 'dim'):
                h=h.rebin("flav",hist.Bin("flav", "flav", [0,1,4,5,6]))
                #ax=plot.plot1d(h["Data"].sum("flav"),overlay="dataset",ax=ax,  density=dense, error_opts=data_err_opts, clear=False);
                #ax = plot.plot1d(h[notdata].sum("flav"), overlay="dataset", ax=ax,  density=dense, stack=True, clear=False);
                # ax = plot.plot1d(h[notdata].integrate("dataset"), overlay="flav", ax=ax,  density=dense, stack=True, clear=False);
                hlist=[h[notdata].integrate("dataset").integrate("flav",slice(5,6)).values()[()],h[notdata].integrate("dataset").integrate("flav",slice(4,5)).values()[()],h[notdata].integrate("dataset").integrate("flav",slice(1,4)).values()[()],h[notdata].integrate("dataset").integrate("flav",slice(0,1)).values()[()]]
                print(h[notdata].integrate("dataset").integrate("flav",slice(4,5)).values()[()])
                print(h[notdata].integrate("dataset").integrate("flav",slice(5,6)).values()[()])
                print(h[notdata].integrate("dataset").integrate("flav",slice(1,4)).values()[()])
                print(h[notdata].integrate("dataset").integrate("flav",slice(0,1)).values()[()])
                
                hlist=[h[notdata].integrate("dataset").integrate("flav",slice(4,5)).values()[()],h[notdata].integrate("dataset").integrate("flav",slice(1,4)).values()[()],h[notdata].integrate("dataset").integrate("flav",slice(5,6)).values()[()],h[notdata].integrate("dataset").integrate("flav",slice(0,1)).values()[()]]
                print(h1name)
                hep.histplot(hlist,h[notdata].axis(h1name).edges(),stack=True,ax=ax,histtype="fill",label=['c', 'pileup','b','udsg'],color=['tab:green','tab:orange','tab:red','tab:blue'])
                hep.histplot(h["Data"].integrate("dataset").integrate("flav").values()[()],h["Data"].axis(h1name).edges(),ax=ax,histtype="errorbar",label="data",yerr=True,color='black')
                
                # plot.plot1d(h["Data"].integrate("flav"),overlay="dataset",ax=ax,  density=dense, error_opts=data_err_opts, clear=False);
                #plot.plot1d(h["Data"],overlay="dataset",ax=ax,  density=dense, error_opts=data_err_opts, clear=False);
    #             plot.plot1d(h["SingleMuon"].sum("flav"),overlay="dataset",ax=ax,  density=dense, error_opts=data_err_opts, clear=False);
    #             plot.plot1d(h["SingleMuon"].sum("dataset"),overlay="flav",ax=ax,  density=dense, clear=False);
    #             ax.legend(handles=ax.get_legend_handles_labels()[0],labels=['QCD80-120','QCD800-1000','QCD600-800','QCD50-80','QCD470-600','QCD30-50','QCD300-470','QCD20-30','QCD170-300','QCD15-20','QCD120-170','QCD1000-inf','data'])
                # ax.legend(handles=ax.get_legend_handles_labels()[0],labels=['W+Jets','tt_semi','tt_had','tt_lep','tW_top','tW_antitop','ST_t_top','ST_t_antitop','ST_s_top_Lep','ST_s_top_Had','QCD80-120','QCD800-1000','QCD600-800','QCD50-80','QCD470-600','QCD30-50','QCD300-470','QCD20-30','QCD170-300','QCD15-20','QCD120-170','QCD1000-inf','DY','data'],fontsize=10,ncol=3)
                # ax.legend(handles=ax.get_legend_handles_labels()[0],labels=['ZZ','WZ','WW','W+Jets','tt_semi','tt_had','tt_lep','tW_top','tW_antitop','ST_t_top','ST_t_antitop','DY','data'],ncol=2)
                # ax.legend(handles=ax.get_legend_handles_labels()[0],labels=['c','pileup','b','udsg','data'])
                ax.legend(loc="upper right")
    #             now we build the ratio plot
                # print(ax.get_xlim())
                
                plot.plotratio(
                        num=h["Data"].sum("dataset").sum("flav"),
                        denom=h[notdata].sum("dataset").sum("flav"),
                        ax=rax,
                        error_opts=data_err_opts,
                        denom_fill_opts={},
                        guide_opts={},
                        unc='num'
                    )         

            else:
                continue   

        for ax, rax, hname in zip([ax1], [rax1], [h1name]):
    #         if not any([h1name.startswith('btagCSVV2')]): 
    #             continue
            at = AnchoredText(r"$1\mu, 1e$"+"\n"+
                               "+$\geq$2 jets"+"\n"+
                              r"$|\eta| < 2.5$",
                               loc=2, frameon=False)
            # at = AnchoredText(r"$1\mu RelIso\leq0.12$"+"\n"+
            #                 r"$\geq$4 AK4 jets"+"\n",loc=2, frameon=False)
    #         at = AnchoredText(r"2+ jets"+"\n",
    #                            loc=2, frameon=False)
            ax.add_artist(at)
            ax.set_ylim(0.01,1e6)
            #if ("btag" in hname) or hname.startswith("DeepCSV_trackDecayLenVal") or (hname=="pt") or hname.startswith("DeepCSV_trackDeltaR"):
            ax.semilogy()
            rax.set_ylabel('Data/Pred.')
            rax.set_xlabel(hname)
            rax.set_ylim(0.5,1.5)
            # if(hname is "eta"):
            # ax.set_xlim(-2.4,2.4)
            # rax.set_xlim(-2.4,2.4)
            #ax.set_xlim(0,1)
            #rax.set_xlim(0,1)
            ax.set_xlabel(None)
        
        hep.mpl_magic(ax1)
        # hep.mpl_magic(rax1)
        # plt.show()
        # fig.savefig("ctag_med_sample_%s_%s.pdf"%(hname,r))
        fig.savefig("plot/ctag_rawpT_flav_%s.pdf"%(hname))
