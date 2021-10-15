import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from coffea.util import load
from coffea.hist import plot
from coffea import hist
import re
import argparse


from cycler import cycler

# colors = ["#F44336","#E91E63","#9C27B0","#673AB7","#3F51B5","#2196F3","#03A9F4","#00BCD4","#009688","#4CAF50","#8BC34A","#CDDC39","#FFEB3B","#FFC107","#FF9800","#FF5722","#795548","#BDBDBD","#9E9E9E","#616161","#90BED4","#607D8B","#455A64"]
markers = ['o','^','s','+','x','D','*']
parser = argparse.ArgumentParser(description='hist plotter for commissioning')
parser.add_argument('-p','--phase',required=True,choices=['dilep','ctag','Wc','DY','ttsemilep','ttdilep'],dest='phase',help='which phase space')
# parser.add_argument('--ext', type=str, default='data', help='addional name')
parser.add_argument('-o','--output',required=True, type=str, help='files set')
parser.add_argument('-r' ,'--ref',required=True,help='referance dataset')
parser.add_argument('-c','--compared' ,required=True,nargs='*',default=[],action='append', help='compared dataset')
parser.add_argument('--sepflav', default=False, type = bool, help = 'seperate flavour')
parser.add_argument('--log',type=bool,help='log on x axis')
parser.add_argument('-d','--discr_list',nargs='+', default=['btagDeepCvL','btagDeepCvB','btagDeepFlavCvL','btagDeepFlavCvB','btagDeepB_b','btagDeepFlavB'],help='discriminators')
arg = parser.parse_args()
output=load('hists_%s.coffea'%(arg.output))
print(output)
events = output['sumw']
if arg.phase == 'ttdilep' :
    input_txt = 'dilepton ttbar'
    nj=2
elif arg.phase == 'ttsemilep' : 
    input_txt = 'semileptonic ttbar'
    nj=4
else:
    if arg.phase == "Wc"  : input_txt = "W+c"
    elif arg.phase == "DY" : input_txt = "DY+jets"
    elif arg.phase == "ttsemilep" : input_txt  = "semileptonic ttbar"
    elif arg.phase == "ttdilep": input_txt  = "dileptonic ttbar"
    nj=1 
if 'njet' in arg.discr_list or 'nbjet' in arg.discr_list or 'mu' in arg.discr_list: nj=1
for j in range(nj):
    for discr in arg.discr_list:
        hflav = output['%s' %(discr)] 
        if 'btag' in discr or 'CvL' in discr or 'CvB' in discr or 'DeepCSV' in discr:
            hflav=hflav.rebin("flav",hist.Bin("flav", "flav", [0,1,4,5,6]))
        if arg.sepflav:
            fig, (ax,rax,rax2,rax3) = plt.subplots(4,1, figsize=(6, 6), gridspec_kw={"height_ratios": (3, 1,1,1)}, sharex=True)
            fig.subplots_adjust(hspace=.01)
            
        else:
            fig, ((ax),(rax)) = plt.subplots(2, 1, figsize=(6, 6), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
            fig.subplots_adjust(hspace=.01)
            
        if 'mu'  in discr and 'njet' in discr :
            
            ax=plot.plot1d(hflav[arg.ref].sum("dataset"),error_opts=None,ax=ax,  density=True)  
            leg_label = []
            leg_label.append(arg.ref)
            for c in arg.compared:          
                plot.plot1d(hflav[c].sum("dataset"),error_opts=data_err_opts,ax=ax,clear=False,  density=True)
                leg_label.append(c)
            ax.legend(ncol=2,loc="upper right",handles=ax.get_legend_handles_labels()[0],labels=leg_label,fontsize=13)
            ax.set_xlabel(None)
            # ax.set_xticklabels(ax.get_xticklabels(), fontsize=0)
            
                     
           
            for c in len(arg.compared)-1:
                if(c==0):
                    rax = plot.plotratio(
                                                num=hflav[arg.compared[0]].sum("dataset"),
                                                denom=hflav[arg.ref].sum("dataset"),
                                                ax=rax,
                                                error_opts=data_err_opts,
                                                denom_fill_opts={},
                                                guide_opts={},
                                                unc='num',
                                                clear=False) 
                else:
                    plot.plotratio(
                                                num=hflav[arg.compared[c+1]].sum("dataset"),
                                                denom=hflav[arg.ref].sum("dataset"),
                                                ax=rax,
                                                error_opts=data_err_opts,
                                                denom_fill_opts={},
                                                guide_opts={},
                                                unc='num',
                                                clear=False)    
                                              
        else:
            if arg.sepflav:
                
                
                ax = plot.plot1d(hflav[arg.ref].sum("dataset").integrate("flav",slice(0,4)),ax=ax,error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'color':'b'},density=True,clear=False)
                plot.plot1d(hflav[arg.ref].sum("dataset").integrate("flav",slice(4,5)),ax=ax,error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'color':'g'},density=True,clear=False)
                plot.plot1d(hflav[arg.ref].sum("dataset").integrate("flav",slice(5,6)),ax=ax,error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'color':'r'},density=True,clear=False)
                leg_label = []
                leg_label.append("%s-l" %arg.ref)    
                leg_label.append("%s-c" %arg.ref)    
                leg_label.append("%s-b" %arg.ref)    

                for c in range(len(arg.compared[0])):     
                    leg_label.append("%s-l" %arg.compared[0][c])
                    leg_label.append("%s-c" %arg.compared[0][c])
                    leg_label.append("%s-b" %arg.compared[0][c])
                    plot.plot1d(hflav[arg.compared[0][c]].sum("dataset").integrate("flav",slice(0,4)),ax=ax,error_opts={'linestyle': 'none','marker': markers[c], 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'b'},density=True)
                    plot.plot1d(hflav[arg.compared[0][c]].sum("dataset").integrate("flav",slice(4,5)),ax=ax,error_opts={'linestyle': 'none','marker': markers[c], 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'g'},density=True,clear=False)
                    plot.plot1d(hflav[arg.compared[0][c]].sum("dataset").integrate("flav",slice(5,6)),ax=ax,error_opts={'linestyle': 'none','marker': markers[c], 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'r'},density=True,clear=False)
                  
                
                ax.legend(ncol=3,loc="upper right",handles=ax.get_legend_handles_labels()[0],labels=leg_label,fontsize=13)
                
                for c in range(len(arg.compared[0])):
                    if c ==0 :
                        rax = plot.plotratio(
                                                        num=hflav[arg.compared[0]].sum("dataset").integrate("flav",slice(0,4)),
                                                        denom=hflav[arg.ref].sum("dataset").integrate("flav",slice(0,4)),
                                                        ax=rax,
                                                        error_opts={'linestyle': 'none','marker': markers[c], 'markersize': 5.,'color':'b'},
                                                        denom_fill_opts={},
                                                        guide_opts={},
                                                        unc='num',
                                                        clear=False)
                        rax2 = plot.plotratio(
                                                        num=hflav[arg.compared[0]].sum("dataset").integrate("flav",slice(4,5)),
                                                        denom=hflav[arg.ref].sum("dataset").integrate("flav",slice(4,5)),
                                                        ax=rax2,
                                                        error_opts={'linestyle': 'none','marker': markers[c], 'markersize': 5.,'color':'g'},
                                                        denom_fill_opts={},
                                                        guide_opts={},
                                                        unc='num',
                                                        clear=False)
                        rax3 = plot.plotratio(
                                                        num=hflav[arg.compared[0]].sum("dataset").integrate("flav",slice(5,6)),
                                                        denom=hflav[arg.ref].sum("dataset").integrate("flav",slice(5,6)),
                                                        ax=rax3,
                                                        error_opts={'linestyle': 'none','marker': markers[c], 'markersize': 5.,'color':'r'},
                                                        denom_fill_opts={},
                                                        guide_opts={},
                                                        unc='num',
                                                        clear=False)
                    else:
                        plot.plotratio(
                                                        num=hflav[arg.compared[0]].sum("dataset").integrate("flav",slice(0,4)),
                                                        denom=hflav[arg.ref].sum("dataset").integrate("flav",slice(0,4)),
                                                        ax=rax,
                                                        error_opts={'linestyle': 'none','marker': markers[c], 'markersize': 5.,'color':'b'},
                                                        denom_fill_opts={},
                                                        guide_opts={},
                                                        unc='num',
                                                        clear=False)
                        plot.plotratio(
                                                        num=hflav[arg.compared[0]].sum("dataset").integrate("flav",slice(4,5)),
                                                        denom=hflav[arg.ref].sum("dataset").integrate("flav",slice(4,5)),
                                                        ax=rax2,
                                                        error_opts={'linestyle': 'none','marker': markers[c], 'markersize': 5.,'color':'g'},
                                                        denom_fill_opts={},
                                                        guide_opts={},
                                                        unc='num',
                                                        clear=False)
                        plot.plotratio(
                                                        num=hflav[arg.compared[0]].sum("dataset").integrate("flav",slice(5,6)),
                                                        denom=hflav[arg.ref].sum("dataset").integrate("flav",slice(5,6)),
                                                        ax=rax3,
                                                        error_opts={'linestyle': 'none','marker': markers[c], 'markersize': 5.,'color':'r'},
                                                        denom_fill_opts={},
                                                        guide_opts={},
                                                        unc='num',
                                                        clear=False)
                
                discrs=discr
                ax.set_xlabel("A.U.")
                rax3.set_xlabel("%s"%(discrs),fontsize=15)
                rax.set_ylabel('udsg-jets',fontsize=12)
                rax2.set_ylabel('c-jets',fontsize=12)
                rax3.set_ylabel('b-jets',fontsize=12)
                rax.set_ylim(0.5,1.5)
                rax2.set_ylim(0.5,1.5)
                rax3.set_ylim(0.5,1.5)
                # ax.set_ylim(0.1,10)
                # ax.semilogy()
                fig.savefig("plot/%s_lin_inclusive%s_nocut.png" %(arg.phase, discrs))
            else: 
                ax=plot.plot1d(hflav[arg.ref].sum("dataset").sum("flav"),ax=ax,density=True)     
                leg_label = []
                leg_label.append(arg.ref.replace("-02Apr2020_ver2-v1",""))      
                for c in range(len(arg.compared[0])):     
                    leg_label.append(arg.compared[0][c].replace("-02Apr2020-v1",""))
                    plot.plot1d(hflav[arg.compared[0][c]].sum("dataset").sum("flav"),ax=ax,clear=False,  density=True)
                ax.legend(ncol=3,loc="upper right",handles=ax.get_legend_handles_labels()[0],labels=leg_label,fontsize=10)
                ax.set_xlabel(None)
                # ax.set_xticklabels(ax.get_xticklabels(), fontsize=0)
                
                
                for c in range(len(arg.compared)-1) :
                    if(c==0):
                        scale_sf=sum(hflav_nosf[notdata].integrate("dataset").integrate("flav").integrate("char").values()[()])/sum(hflav[notdata].integrate("dataset").integrate("flav").integrate("char").values()[()])
                        rax = plot.plotratio(
                                                        num=hflav[arg.compared[0][0]].sum("dataset").sum("flav"),
                                                        denom=hflav[arg.ref].sum("dataset").sum("flav"),
                                                        ax=rax,
                                                        denom_fill_opts={},
                                                        guide_opts={},
                                                        unc='num',
                                                        clear=False)
                    else:
                        plot.plotratio(
                                                         num=hflav[arg.compared[0][c]].sum("dataset").sum("flav"),
                                                         denom=hflav[arg.ref].sum("dataset").sum("flav"),
                                                         ax=rax,
                                                         error_opts=data_err_opts,
                                                         denom_fill_opts={},
                                                         guide_opts={},
                                                         unc='num',
                                                         clear=False)
                ax.set_ylabel("Events",fontsize=15)
                rax.set_ylabel('data/MC',fontsize=15)
                ax.set_ylim(0,10.)
                rax.set_ylim(0.0,2.0)
                fig.savefig("plot/%s_lin_inclusive%s.pdf" %(arg.phase, discr))
                
                

