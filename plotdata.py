import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from coffea.util import load
from coffea.hist import plot
from coffea import hist
import re
import argparse

data_err_opts = {
    'linestyle': 'none',
    'marker': '.',
    'markersize': 10.,
    'color': 'k',
    'elinewidth': 1,
}
from cycler import cycler
colors=["#F44336","#E91E63","#9C27B0","#673AB7","#3F51B5","#2196F3","#03A9F4","#00BCD4","#009688","#4CAF50","#8BC34A","#CDDC39","#FFEB3B","#FFC107","#FF9800","#FF5722","#795548","#BDBDBD","#9E9E9E","#616161","#90BED4","#607D8B","#455A64"]
from xs_scaler import scale_xs
ptbin= [25, 30, 40, 60, 80, 100,150,200,300,500]
etabin=[-2.5,-2.0,-1.5,-0.5,0.,0.5,1.5,2.0,2.5]
notdata = re.compile('(?!data)')
parser = argparse.ArgumentParser(description='hist plotter for commissioning')
parser.add_argument('--lumi',required=True,type=float,help='luminosity in /pb')
parser.add_argument('-c','--combine', default=True, type=bool,help='combined all the jets')
parser.add_argument('-p','--phase',required=True,choices=['dilep','ctag'],dest='phase',help='which phase space')
parser.add_argument('--log',type=bool,help='log on y axis')
parser.add_argument('--norm',default=True,type=bool,help='Use for reshape SF, scale to same yield as no SFs case')
parser.add_argument('-d','--discr_list',nargs='+', default=['btagDeepB', 'btagDeepC', 'btagDeepFlavB', 'btagDeepFlavC','deepcsv_CvL','deepcsv_CvB','deepflav_CvL','deepflav_CvB'],help='discriminators')
# parser.add_argument('-r' '--ref',default='data',help='referance dataset')
# parser.add_argument('--compare',nargs='+',default=['data_runB','data_runC','data_runD','data_runE','data_runF'],help='list to compare to reference')
parser.add_argument('--ext', default='data',type=str, help='addional name')

arg = parser.parse_args()
datas=re.compile('(?=%s)'%(arg.ext))
output=load('hists_%s_sf_%s_AK4.coffea' %(arg.phase,arg.phase))
lists={
    'data_runB':41.5/4.792,
    'data_runC':41.5/9.755,
    'data_runD':41.5/4.319,
    'data_runE':41.5/9.,
    'data_runF':41.5/12.5
}
events = output['sumw']
if arg.phase == 'dilep' :
    input_txt = 'dilepton'
    nj=2
elif arg.phase == 'ctag' : 
    input_txt = 'semileptonic'
    nj=4
if arg.combine : nj=1
for j in range(nj):
    for discr in arg.discr_list:
        if 'btag' in arg.discr_list or 'CvL' in arg.discr_list or 'CvB' in arg.discr_list:
            if arg.combine:
                hflav_0 = output['%sSF_0' %(discr)]
                hflav_nosf0 = output['%s_0' %(discr)]
                hflav_up0 = output['%s_up_0' %(discr)]
                hflav_dn0 = output['%s_dn_0' %(discr)]
                hflav_1 = output['%sSF_1' %(discr)]
                hflav_nosf1 = output['%s_1' %(discr)]
                hflav_up1 = output['%s_up_1' %(discr)]
                hflav_dn1 = output['%s_dn_1' %(discr)]
                hflav=hflav_0+hflav_1
                hflav_nosf=hflav_nosf0+hflav_nosf1
                hflav_up=hflav_up0+hflav_up1
                hflav_dn=hflav_dn0+hflav_dn1

                if(arg.phase=='ctag'):
                    hflav_2 = output['%sSF_2' %(discr)]
                    hflav_nosf2 = output['%s_2' %(discr)]
                    hflav_up2 = output['%s_up_2' %(discr)]
                    hflav_dn2 = output['%s_dn_2' %(discr)]
                    hflav_3 = output['%sSF_3' %(discr)]
                    hflav_nosf3 = output['%s_3' %(discr)]
                    hflav_up3 = output['%s_up_3' %(discr)]
                    hflav_dn3 = output['%s_dn_3' %(discr)]
                    hflav=hflav_0+hflav_1+hflav_2+hflav_3
                    hflav_nosf=hflav_nosf0+hflav_nosf1+hflav_nosf2+hflav_nosf3
                    hflav_up=hflav_up0+hflav_up1+hflav_up2+hflav_up3
                    hflav_dn=hflav_dn0+hflav_dn1+hflav_dn2+hflav_dn3  
                else:
                    hflav = output['%sSF_%d' %(discr,j)]     
                    hflav_nosf = output['%s_%d' %(discr,j)]     
                    hflav_up = output['%s_up_%d' %(discr,j)]
                    hflav_dn = output['%s_dn_%d' %(discr,j)]
        else: hflav = output[discr]
        hflav=hflav.rebin("flav",hist.Bin("flav", "flav", [0,1,4,5,6]))
        if 'btag' in arg.discr_list or 'CvL' in arg.discr_list or 'CvB' in arg.discr_list:scale_sf=sum(hflav_nosf[notdata].integrate("dataset").integrate("flav").integrate("ptwide").integrate("etawide").values()[()])/sum(hflav[notdata].integrate("dataset").integrate("flav").integrate("ptwide").integrate("etawide").values()[()])
        else : scale_sf=1.
        if not arg.norm : scale_sf=1.
        hflav=scale_xs(hflav,arg.lumi*scale_sf,events) 
        if 'btag' in arg.discr_list or 'deepflav' in arg.discr_list or 'deepcsv' in arg.discr_list:
            hflav_nosf=scale_xs(hflav_nosf,arg.lumi,events)  
            hflav_up=scale_xs(hflav_up,arg.lumi*scale_sf,events)   
            hflav_dn=scale_xs(hflav_dn,arg.lumi*scale_sf,events) 
        fig, ((ax),(rax)) = plt.subplots(2, 1, figsize=(8, 6), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
        fig.subplots_adjust(hspace=.07)
        if 'btag' in arg.discr_list or 'CvL' in arg.discr_list or 'CvB' in arg.discr_list:
            if('C' in discr):
                err_up=hflav_up[notdata].integrate("dataset").integrate("flav").integrate("ptwide").integrate("etawide").values()[()]
                err_dn=hflav_dn[notdata].integrate("dataset").integrate("flav").integrate("ptwide").integrate("etawide").values()[()]
            else:
                err_up=np.sqrt(np.add(np.power(np.add(hflav_up[notdata].integrate("dataset").integrate("flav").integrate("ptwide").integrate("etawide").values()[()],-1.*hflav[notdata].integrate("dataset").integrate("flav").integrate("ptwide").integrate("etawide").values()[()]),2),np.power(2*np.add(hflav_up[notdata].integrate("dataset").integrate("flav",slice(5,6)).integrate("ptwide").integrate("etawide").values()[()],-1.*hflav[notdata].integrate("dataset").integrate("flav",slice(5,6)).integrate("ptwide").integrate("etawide").values()[()]),2)))
                err_dn=np.sqrt(np.add(np.power(np.add(hflav[notdata].integrate("dataset").integrate("flav").integrate("ptwide").integrate("etawide").values()[()],-1.*hflav_dn[notdata].integrate("dataset").integrate("flav").integrate("ptwide").integrate("etawide").values()[()]),2),np.power(2*np.add(hflav[notdata].integrate("dataset").integrate("flav",slice(5,6)).integrate("ptwide").integrate("etawide").values()[()],-1.*hflav_dn[notdata].integrate("dataset").integrate("flav",slice(5,6)).integrate("ptwide").integrate("etawide").values()[()]),2)))
            data=hflav[datas].integrate("dataset").integrate("flav").integrate("ptwide").integrate("etawide").values()[()]
            # if not arg.log:
            
            
            ratio_up=np.divide(err_up,data,out=np.zeros_like(err_up), where=data!=0)
            ratio_dn=np.divide(err_dn,data,out=np.zeros_like(err_up), where=data!=0)
        else :
            data=hflav[datas].integrate("dataset").integrate("flav").values()[()]
            
        if 'btag' in arg.discr_list or 'CvL' in arg.discr_list or 'CvB' in arg.discr_list:
            ax=plot.plot1d(hflav[notdata].sum("dataset").sum("ptwide").sum("etawide").sum("flav"),ax=ax,clear=False)
            plot.plot1d(hflav[notdata].sum("dataset",sumw2=True).sum("flav").sum("ptwide").sum("etawide"), ax=ax, density=False, clear=False,error_opts={'linestyle': 'none',
                'markersize':0,
                'elinewidth': 10,
                'color':'tab:brown',
                'alpha':.3,
                'yerr':[err_up,err_dn]})
            plot.plot1d(hflav[notdata].sum("dataset").sum("flav").sum("ptwide").sum("etawide"),ax=ax,error_opts={'linestyle': 'none',
                'markersize':0,
                'elinewidth': 10,
                'color':'tab:gray',
                'alpha':.3}, clear=False,  density=False)
            
            plot.plot1d(hflav[datas].sum("dataset").sum("flav").sum("ptwide").sum("etawide"),error_opts=data_err_opts,ax=ax,clear=False,  density=False)
            
            rax = plot.plotratio(
                                            num=hflav[datas].sum("dataset").sum("flav").sum("etawide").sum("ptwide"),
                                            denom=hflav[notdata].sum("dataset").sum("flav").sum("etawide").sum("ptwide"),
                                            ax=rax,
                                            error_opts= {'linestyle': 'none','marker': '.', 'markersize': 0.,'color':'k'},
                                            denom_fill_opts={'yerr':[ratio_up,ratio_dn]},
                                            # denom_fill_opts={},
                                            guide_opts={},
                                            unc='num',
                                            clear=False)
            plot.plotratio(
                                            num=hflav[datas].sum("dataset").sum("flav").sum("etawide").sum("ptwide"),
                                            denom=hflav[notdata].sum("dataset").sum("flav").sum("etawide").sum("ptwide"),
                                            ax=rax,
                                            error_opts=data_err_opts,
                                            denom_fill_opts={},
                                            guide_opts={},
                                            unc='num',
                                            clear=False)
            hflav.scale(lists,axis="dataset")

            # print("?")
            plot.plot1d(hflav['data_runB'].sum("dataset").sum("flav").sum("ptwide").sum("etawide"),error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:orange'},ax=ax,clear=False,  density=False)
            plot.plot1d(hflav['data_runC'].sum("dataset").sum("flav").sum("ptwide").sum("etawide"),error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:green'},ax=ax,clear=False,  density=False)
            plot.plot1d(hflav['data_runD'].sum("dataset").sum("flav").sum("ptwide").sum("etawide"),error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:red'},ax=ax,clear=False,  density=False)
            plot.plot1d(hflav['data_runE'].sum("dataset").sum("flav").sum("ptwide").sum("etawide"),error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:purple'},ax=ax,clear=False,  density=False)
            plot.plot1d(hflav['data_runF'].sum("dataset").sum("flav").sum("ptwide").sum("etawide"),error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:pink'},ax=ax,clear=False,  density=False)

            ax.legend(ncol=3,loc="upper right",handles=ax.get_legend_handles_labels()[0],labels=['MC','SFs Unc.','stat. Unc.','All Data','runB','runC','runD','runE','runF'],fontsize=12)
            
            plot.plotratio(
                                            num=hflav['data_runB'].sum("dataset").sum("flav").sum("etawide").sum("ptwide"),
                                            denom=hflav[notdata].sum("dataset").sum("flav").sum("etawide").sum("ptwide"),
                                            ax=rax,
                                            error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:orange'},
                                            denom_fill_opts={},
                                            guide_opts={},
                                            unc='num',
                                            clear=False)
            plot.plotratio(
                                            num=hflav['data_runC'].sum("dataset").sum("flav").sum("etawide").sum("ptwide"),
                                            denom=hflav[notdata].sum("dataset").sum("flav").sum("etawide").sum("ptwide"),
                                            ax=rax,
                                            error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:green'},
                                            denom_fill_opts={},
                                            guide_opts={},
                                            unc='num',
                                            clear=False)
            plot.plotratio(
                                            num=hflav['data_runD'].sum("dataset").sum("flav").sum("etawide").sum("ptwide"),
                                            denom=hflav[notdata].sum("dataset").sum("flav").sum("etawide").sum("ptwide"),
                                            ax=rax,
                                            error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:red'},
                                            denom_fill_opts={},
                                            guide_opts={},
                                            unc='num',
                                            clear=False)
            plot.plotratio(
                                            num=hflav['data_runE'].sum("dataset").sum("flav").sum("etawide").sum("ptwide"),
                                            denom=hflav[notdata].sum("dataset").sum("flav").sum("etawide").sum("ptwide"),
                                            ax=rax,
                                            error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:purple'},
                                            denom_fill_opts={},
                                            guide_opts={},
                                            unc='num',
                                            clear=False)
            plot.plotratio(
                                            num=hflav['data_runF'].sum("dataset").sum("flav").sum("etawide").sum("ptwide"),
                                            denom=hflav[notdata].sum("dataset").sum("flav").sum("etawide").sum("ptwide"),
                                            ax=rax,
                                            error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:pink'},
                                            denom_fill_opts={},
                                            guide_opts={},
                                            unc='num',
                                            clear=False)
        else :
                ax=plot.plot1d(hflav[notdata].sum("dataset").sum("flav"),ax=ax,clear=False)                
                plot.plot1d(hflav[datas].sum("dataset").sum("flav"),error_opts=data_err_opts,ax=ax,clear=False,  density=False)
                rax = plot.plotratio(
                                                num=hflav[datas].sum("dataset").sum("flav"),
                                                denom=hflav[notdata].sum("dataset").sum("flav"),
                                                ax=rax,
                                                error_opts=data_err_opts,
                                                denom_fill_opts={},
                                                guide_opts={},
                                                unc='num',
                                                clear=False)
                hflav.scale(lists,axis="dataset")

                
                plot.plot1d(hflav['data_runB'].sum("dataset").sum("flav"),error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:orange'},ax=ax,clear=False,  density=False)
                plot.plot1d(hflav['data_runC'].sum("dataset").sum("flav"),error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:green'},ax=ax,clear=False,  density=False)
                plot.plot1d(hflav['data_runD'].sum("dataset").sum("flav"),error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:red'},ax=ax,clear=False,  density=False)
                plot.plot1d(hflav['data_runE'].sum("dataset").sum("flav"),error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:purple'},ax=ax,clear=False,  density=False)
                plot.plot1d(hflav['data_runF'].sum("dataset").sum("flav"),error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:pink'},ax=ax,clear=False,  density=False)

                ax.legend(ncol=3,loc="upper right",handles=ax.get_legend_handles_labels()[0],labels=['MC','All Data','runB','runC','runD','runE','runF'],fontsize=12)
                
                plot.plotratio(
                                                num=hflav['data_runB'].sum("dataset").sum("flav"),
                                                denom=hflav[notdata].sum("dataset").sum("flav"),
                                                ax=rax,
                                                error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:orange'},
                                                denom_fill_opts=None,
                                                guide_opts={},
                                                unc='num',
                                                clear=False)
                plot.plotratio(
                                                num=hflav['data_runC'].sum("dataset").sum("flav"),
                                                denom=hflav[notdata].sum("dataset").sum("flav"),
                                                ax=rax,
                                                error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:green'},
                                                denom_fill_opts=None,
                                                guide_opts={},
                                                unc='num',
                                                clear=False)
                plot.plotratio(
                                                num=hflav['data_runD'].sum("dataset").sum("flav"),
                                                denom=hflav[notdata].sum("dataset").sum("flav"),
                                                ax=rax,
                                                error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:red'},
                                                denom_fill_opts=None,
                                                guide_opts={},
                                                unc='num',
                                                clear=False)
                plot.plotratio(
                                                num=hflav['data_runE'].sum("dataset").sum("flav"),
                                                denom=hflav[notdata].sum("dataset").sum("flav"),
                                                ax=rax,
                                                error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:purple'},
                                                denom_fill_opts=None,
                                                guide_opts={},
                                                unc='num',
                                                clear=False)
                plot.plotratio(
                                                num=hflav['data_runF'].sum("dataset").sum("flav"),
                                                denom=hflav[notdata].sum("dataset").sum("flav"),
                                                ax=rax,
                                                error_opts={'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none', 'elinewidth': 1.5,'color':'tab:pink'},
                                                denom_fill_opts=None,
                                                guide_opts={},
                                                unc='num',
                                                clear=False)
        if 'btag' in arg.discr_list or 'CvL' in arg.discr_list or 'CvB' in arg.discr_list:
            maximum=max(max(max(max(max(max(max((hflav[notdata].integrate("dataset").integrate("flav").values()[()])),max(data)),max(hflav["data_runB"].integrate("dataset").integrate("flav").values()[()])),max(hflav["data_runC"].integrate("dataset").integrate("flav").values()[()])),max(hflav["data_runD"].integrate("dataset").integrate("flav").values()[()])),max(hflav["data_runE"].integrate("dataset").integrate("flav").values()[()])),max(hflav["data_runF"].integrate("dataset").integrate("flav").values()[()]))
        else:maximum=max(max(max(max(max(max(max((hflav[notdata].integrate("dataset").integrate("flav").values()[()])),max(data)),max(hflav["data_runB"].integrate("dataset").integrate("flav").values()[()])),max(hflav["data_runC"].integrate("dataset").integrate("flav").values()[()])),max(hflav["data_runD"].integrate("dataset").integrate("flav").values()[()])),max(hflav["data_runE"].integrate("dataset").integrate("flav").values()[()])),max(hflav["data_runF"].integrate("dataset").integrate("flav").values()[()]))
        if arg.log:
            if(arg.phase=="ctag"):ax.set_ylim(100,maximum*10)
            else:ax.set_ylim(1,maximum*10)
        else:ax.set_ylim(0.,maximum*1.5)
        if(arg.log):ax.semilogy()
        ax.set_ylabel("Events",fontsize=15)
        rax.set_ylabel('Data/MC',fontsize=15)

        ax.set_xlabel(None)
        # ax.set_xticklabels(ax.get_xticklabels(), fontsize=0)
        if 'CvL' in discr :discrs=discr.replace('CvL','CvB')
        elif 'CvB' in discr :discrs=discr.replace('CvB','CvL')
        else :discrs=discr
        if arg.combine :rax.set_xlabel("%s"%(discrs),fontsize=15)
        else :rax.set_xlabel("%s[%d]"%(discrs,j),fontsize=15)
        rax.set_ylim(0.5,1.5)
        # at = AnchoredText(r"semileptonic ttbar"+"\n"+
        
        at = AnchoredText(input_txt+" ttbar"+"\n"+
                                "inclusive pT, $\eta$"
                                , loc=2, prop=dict(size=15),frameon=False)
        ax.add_artist(at)
        scale=""
        if arg.norm:scale="_norm"
        if(arg.log):fig.savefig("plot/%s_%s_inclusive%s_%s_all_compare.pdf" %(arg.phase, discrs, scale, arg.ext))
        else:fig.savefig("plot/%s_lin_%s_inclusive%s_%s_all_compare.pdf" %(arg.phase, discrs, scale, arg.ext))
        
        # for i in range(len(ptbin)-1):    
        #     fig, ((ax),(rax)) = plt.subplots(2, 1, figsize=(8, 6), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
        #     fig.subplots_adjust(hspace=.07)
        #     if('C' in discr):
        #         err_up=hflav_up[notdata].integrate("dataset").integrate("flav").integrate("ptwide",slice(ptbin[i],ptbin[i+1])).integrate("etawide").values()[()]
        #         err_dn=hflav_dn[notdata].integrate("dataset").integrate("flav").integrate("ptwide",slice(ptbin[i],ptbin[i+1])).integrate("etawide").values()[()]
        #     else:
        #         err_up=np.sqrt(np.add(np.power(np.add(hflav_up[notdata].integrate("dataset").integrate("flav").integrate("ptwide").integrate("etawide").values()[()],-1.*hflav[notdata].integrate("dataset").integrate("flav").integrate("ptwide").integrate("etawide").values()[()]),2),np.power(2*np.add(hflav_up[notdata].integrate("dataset").integrate("flav",slice(5,6)).integrate("ptwide").integrate("etawide").values()[()],-1.*hflav[notdata].integrate("dataset").integrate("flav",slice(5,6)).integrate("ptwide").integrate("etawide").values()[()]),2)))
        #     err_dn=np.sqrt(np.add(np.power(np.add(hflav[notdata].integrate("dataset").integrate("flav").integrate("ptwide").integrate("etawide").values()[()],-1.*hflav_dn[notdata].integrate("dataset").integrate("flav").integrate("ptwide").integrate("etawide").values()[()]),2),np.power(2*np.add(hflav[notdata].integrate("dataset").integrate("flav",slice(5,6)).integrate("ptwide").integrate("etawide").values()[()],-1.*hflav_dn[notdata].integrate("dataset").integrate("flav",slice(5,6)).integrate("ptwide").integrate("etawide").values()[()]),2)))
        #     data=hflav[datas].integrate("dataset").integrate("flav").integrate("ptwide",slice(ptbin[i],ptbin[i+1])).integrate("etawide").values()[()]
        #     ratio_up=np.divide(err_up,data,out=np.zeros_like(err_up), where=data!=0)
        #     ratio_dn=np.divide(err_dn,data,out=np.zeros_like(err_up), where=data!=0)
        #     hflav=hflav.rebin("flav",hist.Bin("flav", "flav", [0,1,4,5,6]))
        #     ax=plot.plot1d(hflav[notdata].sum("dataset").sum("etawide").integrate("ptwide",slice(ptbin[i],ptbin[i+1])),overlay="flav",fill_opts={},error_opts=None,ax=ax,stack=True)
        #     plot.plot1d(hflav[notdata].sum("dataset",sumw2=True).sum("flav").sum("etawide").integrate("ptwide",slice(ptbin[i],ptbin[i+1])), ax=ax,  density=False, clear=False,error_opts={'linestyle': 'none',
        #     'markersize':0,
        #     'elinewidth': 10,
        #     'color':'tab:brown',
        #     'alpha':.5,
        #     'yerr':[err_up,err_dn]})
        #     plot.plot1d(hflav[notdata].sum("dataset").sum("flav").integrate("ptwide",slice(ptbin[i],ptbin[i+1])).sum("etawide"),ax=ax,error_opts={'linestyle': 'none',
        #         'markersize':0,
        #         'elinewidth': 10,
        #         'color':'tab:gray',
        #         'alpha':.5}, clear=False,  density=False)
        #     plot.plot1d(hflav[notdata].sum("dataset").sum("flav").integrate("ptwide",slice(ptbin[i],ptbin[i+1])).sum("etawide"),ax=ax,line_opts={'color':'tab:blue'},clear=False,  density=False)
        #     plot.plot1d(hflav_nosf[notdata].sum("dataset").integrate("ptwide",slice(ptbin[i],ptbin[i+1])).sum("etawide").sum("flav"),error_opts= {'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none','color' :'tab:pink' , 'elinewidth': 1.5},ax=ax,clear=False)
        #     plot.plot1d(hflav[datas].sum("dataset").sum("flav").integrate("ptwide",slice(ptbin[i],ptbin[i+1])).sum("etawide"),error_opts=data_err_opts,ax=ax,clear=False,  density=False)
        #     ax.legend(ncol=2,handles=ax.get_legend_handles_labels()[0],labels=['b','c','pileup','udsg','SFs Unc.','stat. Unc.','w/o SFs','Data'],fontsize=13)
        #     ax.set_xlabel(None)
        #     ax.set_xticklabels(ax.get_xticklabels(), fontsize=0)
        #     rax = plot.plotratio(
        #                                     num=hflav[datas].sum("dataset").sum("flav").sum("etawide").integrate("ptwide",slice(ptbin[i],ptbin[i+1])),
        #                                     denom=hflav[notdata].sum("dataset").sum("flav").sum("etawide").integrate("ptwide",slice(ptbin[i],ptbin[i+1])),
        #                                     ax=rax,
        #                                     error_opts= {'linestyle': 'none','marker': '.', 'markersize': 0.,'color': 'k'},
        #                                     denom_fill_opts={'yerr':[ratio_up,ratio_dn]},
        #                                     # denom_fill_opts={},
        #                                     guide_opts={},
        #                                     unc='num',
        #                                     clear=False)
        #     plot.plotratio(
        #                                     num=hflav[datas].sum("dataset").sum("flav").sum("etawide").integrate("ptwide",slice(ptbin[i],ptbin[i+1])),
        #                                     denom=hflav[notdata].sum("dataset").sum("flav").sum("etawide").integrate("ptwide",slice(ptbin[i],ptbin[i+1])),
        #                                     ax=rax,
        #                                     # error_opts= {'linestyle': 'none','marker': '.', 'markersize': 10.,'color': 'tab:blue',  'elinewidth': 1},
        #                                     error_opts=data_err_opts,
        #                                     # denom_fill_opts={'yerr':[ratio_up,ratio_dn]},
        #                                     denom_fill_opts={},
        #                                     guide_opts={},
        #                                     unc='num',
        #                                     clear=False)
        #     plot.plotratio(
        #                                 num=hflav_nosf[datas].sum("dataset").sum("flav").sum("etawide").integrate("ptwide",slice(ptbin[i],ptbin[i+1])),
        #                                 denom=hflav_nosf[notdata].sum("dataset").sum("flav").sum("etawide").integrate("ptwide",slice(ptbin[i],ptbin[i+1])),
        #                                 ax=rax,
        #                                 error_opts= {'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none','color' :'tab:pink' , 'elinewidth': 1.5},
        #                                 # error_opts=data_err_opts,
        #                                 # denom_fill_opts={'yerr':[ratio_up,ratio_dn]},
        #                                 denom_fill_opts={},
        #                                 guide_opts={},
        #                                 unc='num',
        #                                 clear=False)
                                
            
        #     if(arg.log):
        #         if(arg.phase=="ctag"):ax.set_ylim(1,10000000)
        #         else:ax.set_ylim(0.1,5000000)
        #     else:ax.set_ylim(0.,maximum*1.2)
        #     if(arg.log):ax.semilogy()
        #     ax.set_ylabel("Events",fontsize=15)
        #     rax.set_ylabel('Data/MC',fontsize=15)
        #     if 'CvL' in discr :discrs=discr.replace('CvL','CvB')
        #     elif 'CvB' in discr :discrs=discr.replace('CvB','CvL')
        #     else :discrs=discr
        #     rax.set_xlabel("%s"%(discrs),fontsize=15)
        #     rax.set_ylim(0.5,1.5)        
        #     at = AnchoredText(input_txt+" ttbar"+"\n"+
        #                     "%.1f$<pT<$%.1f"%(ptbin[i],ptbin[i+1])
        #                     , loc=2, frameon=False,prop=dict(size=15))
        #     ax.add_artist(at)
        #     scale=""
        #     if arg.norm:scale="_norm"
        #     if(arg.log):fig.savefig("plot/%s_unc_%s_pt%d_%d_%s_%s_all.pdf" %(arg.phase, discrs, ptbin[i],ptbin[i+1],scale, arg.ext))
        #     else:fig.savefig("plot/%s_unc_lin_%s_%s_pt%d_%d_%s_%s_all.pdf" %(arg.phase, discrs, ptbin[i],ptbin[i+1],scale, arg.ext))
            
        # # loop etabins
        # for i in range(len(etabin)-1):   
        #     fig, ((ax),(rax)) = plt.subplots(2, 1, figsize=(8, 6), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
        #     fig.subplots_adjust(hspace=.07)
        #     if('C' in discr):
        #         err_up=hflav_up[notdata].integrate("dataset").integrate("flav").integrate("etawide",slice(etabin[i],etabin[i+1])).integrate("ptwide").values()[()]
        #         err_dn=hflav_dn[notdata].integrate("dataset").integrate("flav").integrate("etawide",slice(etabin[i],etabin[i+1])).integrate("ptwide").values()[()]
        #     else:
        #         err_up=np.sqrt(np.add(np.power(np.add(hflav_up[notdata].integrate("dataset").integrate("flav").integrate("ptwide").integrate("etawide").values()[()],-1.*hflav[notdata].integrate("dataset").integrate("flav").integrate("ptwide").integrate("etawide").values()[()]),2),np.power(2*np.add(hflav_up[notdata].integrate("dataset").integrate("flav",slice(5,6)).integrate("ptwide").integrate("etawide").values()[()],-1.*hflav[notdata].integrate("dataset").integrate("flav",slice(5,6)).integrate("ptwide").integrate("etawide").values()[()]),2)))
        #         err_dn=np.sqrt(np.add(np.power(np.add(hflav[notdata].integrate("dataset").integrate("flav").integrate("ptwide").integrate("etawide").values()[()],-1.*hflav_dn[notdata].integrate("dataset").integrate("flav").integrate("ptwide").integrate("etawide").values()[()]),2),np.power(2*np.add(hflav[notdata].integrate("dataset").integrate("flav",slice(5,6)).integrate("ptwide").integrate("etawide").values()[()],-1.*hflav_dn[notdata].integrate("dataset").integrate("flav",slice(5,6)).integrate("ptwide").integrate("etawide").values()[()]),2)))
        #     data=hflav[datas].integrate("dataset").integrate("flav").integrate("etawide",slice(etabin[i],etabin[i+1])).integrate("ptwide").values()[()]
        #     ratio_up=np.divide(err_up,data,out=np.zeros_like(err_up), where=data!=0)
        #     ratio_dn=np.divide(err_dn,data,out=np.zeros_like(err_up), where=data!=0)
        #     ax=plot.plot1d(hflav[notdata].sum("dataset").sum("ptwide").integrate("etawide",slice(etabin[i],etabin[i+1])),overlay="flav",fill_opts={},error_opts=None,ax=ax,stack=True)
        #     plot.plot1d(hflav[notdata].sum("dataset",sumw2=True).sum("flav").sum("ptwide").integrate("etawide",slice(etabin[i],etabin[i+1])), ax=ax,  density=False, clear=False,error_opts={'linestyle': 'none',
        #     'markersize':0,
        #     'elinewidth': 10,
        #     'color':'tab:brown',
        #     'alpha':.5,
        #     'yerr':[err_up,err_dn]})
        #     plot.plot1d(hflav[notdata].sum("dataset").sum("flav").integrate("etawide",slice(etabin[i],etabin[i+1])).sum("ptwide"),ax=ax,error_opts={'linestyle': 'none',
        #         'markersize':0,
        #         'elinewidth': 10,
        #         'color':'tab:gray',
        #         'alpha':.5}, clear=False,  density=False)
        #     plot.plot1d(hflav_nosf[notdata].sum("dataset").sum("ptwide").integrate("etawide",slice(etabin[i],etabin[i+1])).sum("flav"),error_opts= {'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none','color' :'tab:pink' , 'elinewidth': 1.5},ax=ax,clear=False)
        #     plot.plot1d(hflav[datas].sum("dataset").sum("flav").integrate("etawide",slice(etabin[i],etabin[i+1])).sum("ptwide"),error_opts=data_err_opts,ax=ax,clear=False,  density=False)
        #     ax.legend(ncol=2,handles=ax.get_legend_handles_labels()[0],labels=['b','c','pileup','udsg','SFs Unc.','stat. Unc.','w/o SFs','Data'],fontsize=13)

        #     ax.set_xlabel(None)
        #     ax.set_xticklabels(ax.get_xticklabels(), fontsize=0)
        #     rax = plot.plotratio(
        #                                     num=hflav[datas].sum("dataset").sum("flav").sum("ptwide").integrate("etawide",slice(etabin[i],etabin[i+1])),
        #                                     denom=hflav[notdata].sum("dataset").sum("flav").sum("ptwide").integrate("etawide",slice(etabin[i],etabin[i+1])),
        #                                     ax=rax,
        #                                     error_opts= {'linestyle': 'none','marker': '.', 'markersize': 0.,'color': 'k'},
        #                                     denom_fill_opts={'yerr':[ratio_up,ratio_dn]},
        #                                     # denom_fill_opts={},
        #                                     guide_opts={},
        #                                     unc='num',
        #                                     clear=False)
        #     plot.plotratio(
        #                                     num=hflav[datas].sum("dataset").sum("flav").sum("ptwide").integrate("etawide",slice(etabin[i],etabin[i+1])),
        #                                     denom=hflav[notdata].sum("dataset").sum("flav").sum("ptwide").integrate("etawide",slice(etabin[i],etabin[i+1])),
        #                                     ax=rax,
        #                                     # error_opts= {'linestyle': 'none','marker': '.', 'markersize': 10.,'color': 'tab:blue',  'elinewidth': 1},
        #                                     error_opts=data_err_opts,
        #                                     # denom_fill_opts={'yerr':[ratio_up,ratio_dn]},
        #                                     denom_fill_opts={},
        #                                     guide_opts={},
        #                                     unc='num',
        #                                     clear=False)
        #     plot.plotratio(
        #                                     num=hflav_nosf[datas].sum("dataset").sum("flav").sum("ptwide").integrate("etawide",slice(etabin[i],etabin[i+1])),
        #                                     denom=hflav_nosf[notdata].sum("dataset").sum("flav").sum("ptwide").integrate("etawide",slice(etabin[i],etabin[i+1])),
        #                                     ax=rax,
        #                                     error_opts= {'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none','color' :'tab:pink' , 'elinewidth': 1.5},
        #                                     # denom_fill_opts={'yerr':[ratio_up,ratio_dn]},
        #                                     denom_fill_opts={},
        #                                     guide_opts={},
        #                                     unc='num',
        #                                     clear=False)
                                
            
        #     if(arg.log):
        #         if(arg.phase=="ctag"):ax.set_ylim(1,10000000)
        #         else:ax.set_ylim(0.1,5000000)
        #     else:ax.set_ylim(0.,maximum*1.2)
        #     if(arg.log):ax.semilogy()
        #     ax.set_ylabel("Events",fontsize=15)
        #     rax.set_ylabel('Data/MC',fontsize=15)
        #     if 'CvL' in discr :discrs=discr.replace('CvL','CvB')
        #     elif 'CvB' in discr :discrs=discr.replace('CvB','CvL')
        #     else :discrs=discr
        #     rax.set_xlabel("%s"%(discrs),fontsize=15)
        #     rax.set_ylim(0.5,1.5)      
        #     rax.set_xlabel("%s"%(discr),fontsize=15)
        #     at = AnchoredText(input_txt+" ttbar"+"\n"+
        #                     "%.1f$<\eta<$%.1f"%(etabin[i],etabin[i+1])
        #                     , loc=2, frameon=False,prop=dict(size=15))
        #     ax.add_artist(at)
        #     scale=""
        #     if arg.norm:scale="_norm"
        #     if(arg.log):fig.savefig("plot/%s_unc_%s_eta%.1f_%.1f_%s_%s_all.pdf" %(arg.phase, discrs, etabin[i],etabin[i+1],scale, arg.ext))
        #     else:fig.savefig("plot/%s_unc_lin_%s_%s_eta%.1f_%.1f_%s_%s_all.pdf" %(arg.phase, discrs, etabin[i],etabin[i+1],scale, arg.ext))


