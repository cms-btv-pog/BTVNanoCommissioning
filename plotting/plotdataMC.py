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
from utils.xs_scaler import scale_xs
ptbin= [25, 30, 40, 60, 80, 100,150,200,300,500]
etabin=[-2.5,-2.0,-1.5,-0.5,0.,0.5,1.5,2.0,2.5]
notdata = re.compile('(?!data)')
parser = argparse.ArgumentParser(description='hist plotter for commissioning')
parser.add_argument('--lumi',required=True,type=float,help='luminosity in /pb')
parser.add_argument('-c','--combine', type=bool,help='combined all the jets')
parser.add_argument('-p','--phase',required=True,choices=['dilep_sf','ttsemilep_sf','ctag_Wc_sf','ctag_DY_sf','ctag_ttsemilep_sf','ctag_ttdilep_sf'],dest='phase',help='which phase space')
parser.add_argument('--log',type=bool,help='log on x axis')
parser.add_argument('--norm',default=True,type=bool,help='Use for reshape SF, scale to same yield as no SFs case')
parser.add_argument('-d','--discr_list',nargs='+', default=['deepcsv_CvL','deepcsv_CvB','deepflav_CvL','deepflav_CvB','btagDeepB','btagDeepC','btagDeepFlavB','btagDeepFlavC'],help='discriminators')

parser.add_argument('--ext', type=str, default='data', help='addional name')

arg = parser.parse_args()
datas=re.compile('(?=%s)'%(arg.ext))


output=load('hists_ttsemilep_sf_ctag_AK4_runD.coffea')

events = output['sumw']
if arg.phase == 'dilep' :
    input_txt = 'dilepton ttbar'
    nj=2
elif arg.phase == 'ctag' : 
    input_txt = 'semileptonic ttbar'
    nj=4
else:
    if arg.phase == "Wc"  : input_txt = "W+c"
    elif arg.phase == "DY" : input_txt = "DY+jets"
    elif arg.phase == "ttsemilep" : input_txt  = "semileptonic ttbar"
    elif arg.phase == "ttdilep": input_txt  = "dileptonic ttbar"
    nj=1 
if arg.combine : nj=1
if 'njet' in arg.discr_list or 'nbjet' in arg.discr_list or 'mu' in arg.discr_list: nj=1


for j in range(nj):
    for discr in arg.discr_list:
        if arg.combine:
            hflav_0 = output['%sSF_0' %(discr)]
            hflav_1 = output['%sSF_1' %(discr)]
            hflav=hflav_0+hflav_1
            if 'btag' in arg.discr_list or 'CvL' in arg.discr_list or 'CvB' in arg.discr_list:
                hflav_nosf0 = output['%s_0' %(discr)]
                hflav_up0 = output['%s_up_0' %(discr)]
                hflav_dn0 = output['%s_dn_0' %(discr)]
                hflav_nosf1 = output['%s_1' %(discr)]
                hflav_up1 = output['%s_up_1' %(discr)]
                hflav_dn1 = output['%s_dn_1' %(discr)]
                hflav_nosf=hflav_nosf0+hflav_nosf1
                hflav_up=hflav_up0+hflav_up1
                hflav_dn=hflav_dn0+hflav_dn1    

                if(arg.phase=='ctag'):
                    hflav_2 = output['%sSF_2' %(discr)]
                    hflav_3 = output['%sSF_3' %(discr)]
                    hflav=hflav_0+hflav_1+hflav_2+hflav_3
                    if 'btag' in arg.discr_list or 'CvL' in arg.discr_list or 'CvB' in arg.discr_list:
                        hflav_nosf2 = output['%s_2' %(discr)]
                        hflav_up2 = output['%s_up_2' %(discr)]
                        hflav_dn2 = output['%s_dn_2' %(discr)]
                        hflav_nosf3 = output['%s_3' %(discr)]
                        hflav_up3 = output['%s_up_3' %(discr)]
                        hflav_dn3 = output['%s_dn_3' %(discr)]
                        hflav_nosf=hflav_nosf0+hflav_nosf1+hflav_nosf2+hflav_nosf3
                        hflav_up=hflav_up0+hflav_up1+hflav_up2+hflav_up3
                        hflav_dn=hflav_dn0+hflav_dn1+hflav_dn2+hflav_dn3  
        else:
            if 'btag' in discr or 'CvL' in discr or 'CvB' in discr:
                hflav = output['%sSF_%d' %(discr,j)]     
                hflav_up = output['%s_up_%d' %(discr,j)]
                hflav_dn = output['%s_dn_%d' %(discr,j)]
                hflav_nosf = output['%s_%d' %(discr,j)]  
            elif 'nbjet' in arg.discr_list:
                hflav = output[discr]
                hflav_up=   output['%s_up' %(discr)]             
                hflav_dn=   output['%s_dn' %(discr)]   
            else : hflav = output['%s' %(discr)] 

        if 'mu' not in arg.discr_list and 'njet' not in arg.discr_list and 'nbjet'  not in arg.discr_list:
            hflav=hflav.rebin("flav",hist.Bin("flav", "flav", [0,1,4,5,6]))

        if ('btag' in discr or 'CvL' in discr or 'CvB' in discr) and arg.norm : 
            if arg.phase == "Wc":
                scale_sf=sum(hflav_nosf[notdata].integrate("dataset").integrate("flav").integrate("char").values()[()])/sum(hflav[notdata].integrate("dataset").integrate("flav").integrate("char").values()[()])
                
            else :
                scale_sf=sum(hflav_nosf[notdata].integrate("dataset").integrate("flav").values()[()])/sum(hflav[notdata].integrate("dataset").integrate("flav").values()[()])
                

        
        if not arg.norm : scale_sf=1.
        hflav=scale_xs(hflav,arg.lumi*scale_sf,events) 
        if 'btag' in discr or 'CvL' in discr or 'CvB' in discr:
            hflav_nosf=scale_xs(hflav_nosf,arg.lumi,events)  
            hflav_up=scale_xs(hflav_up,arg.lumi*scale_sf,events)   
            hflav_dn=scale_xs(hflav_dn,arg.lumi*scale_sf,events) 
            if arg.phase == "Wc":
                print("============")
                print(sum(hflav_nosf[notdata].integrate("dataset").integrate("flav").integrate("char").values()[()]))
                print("b:%.3f\%",sum(hflav_nosf[notdata].integrate("dataset").integrate("flav",slice(5,6)).integrate("char").values()[()])/sum(hflav_nosf[notdata].integrate("dataset").integrate("flav").integrate("char").values()[()]))
                print("c:%.3f\%",sum(hflav_nosf[notdata].integrate("dataset").integrate("flav",slice(4,5)).integrate("char").values()[()])/sum(hflav_nosf[notdata].integrate("dataset").integrate("flav").integrate("char").values()[()]))
                print("l:%.3f\%",sum(hflav_nosf[notdata].integrate("dataset").integrate("flav",slice(0,4)).integrate("char").values()[()])/sum(hflav_nosf[notdata].integrate("dataset").integrate("flav").integrate("char").values()[()]))
                print("============")

            else :
                print("============")
                print(sum(hflav_nosf[notdata].integrate("dataset").integrate("flav").values()[()]))
                print("b:%.3f",sum(hflav_nosf[notdata].integrate("dataset").integrate("flav",slice(5,6)).values()[()])/sum(hflav_nosf[notdata].integrate("dataset").integrate("flav").values()[()]))
                print("c:%.3f",sum(hflav_nosf[notdata].integrate("dataset").integrate("flav",slice(4,5)).values()[()])/sum(hflav_nosf[notdata].integrate("dataset").integrate("flav").values()[()]))
                print("l:%.3f",sum(hflav_nosf[notdata].integrate("dataset").integrate("flav",slice(0,4)).values()[()])/sum(hflav_nosf[notdata].integrate("dataset").integrate("flav").values()[()]))
                print("============")
        elif 'nbjet' in arg.discr_list:
            hflav_up=scale_xs(hflav_up,arg.lumi,events)   
            hflav_dn=scale_xs(hflav_dn,arg.lumi,events) 
        fig, ((ax),(rax)) = plt.subplots(2, 1, figsize=(6, 6), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
        fig.subplots_adjust(hspace=.07)
        if 'btag' in discr or 'CvL' in discr or 'CvB' in discr:
            if arg.phase == "Wc" :
                if('C' in discr):
                    err_up=hflav_up[notdata].integrate("dataset").integrate("flav").integrate("char").values()[()]
                    err_dn=hflav_dn[notdata].integrate("dataset").integrate("flav").integrate("char").values()[()]
                else:
                    err_up=np.sqrt(np.add(np.power(np.add(hflav_up[notdata].integrate("dataset").integrate("flav").integrate("char").values()[()],-1.*hflav[notdata].integrate("dataset").integrate("flav").integrate("char").values()[()]),2),np.power(2*np.add(hflav_up[notdata].integrate("dataset").integrate("flav",slice(5,6)).integrate("char").values()[()],-1.*hflav[notdata].integrate("dataset").integrate("flav",slice(5,6)).integrate("char").values()[()]),2)))
                    err_dn=np.sqrt(np.add(np.power(np.add(hflav[notdata].integrate("dataset").integrate("flav").integrate("char").values()[()],-1.*hflav_dn[notdata].integrate("dataset").integrate("flav").integrate("char").values()[()]),2),np.power(2*np.add(hflav[notdata].integrate("dataset").integrate("flav",slice(5,6)).integrate("char").values()[()],-1.*hflav_dn[notdata].integrate("dataset").integrate("flav",slice(5,6)).integrate("char").values()[()]),2)))
                data=hflav[datas].integrate("dataset").integrate("flav").integrate("char").values()[()]
                maximum=max(max((hflav[notdata].integrate("dataset").integrate("flav").integrate("char").values()[()]+err_up)),max(data))
            else :
                if('C' in discr):
                    err_up=hflav_up[notdata].integrate("dataset").integrate("flav").values()[()]
                    err_dn=hflav_dn[notdata].integrate("dataset").integrate("flav").values()[()]
                else:
                    err_up=np.sqrt(np.add(np.power(np.add(hflav_up[notdata].integrate("dataset").integrate("flav").values()[()],-1.*hflav[notdata].integrate("dataset").integrate("flav").values()[()]),2),np.power(2*np.add(hflav_up[notdata].integrate("dataset").integrate("flav",slice(5,6)).values()[()],-1.*hflav[notdata].integrate("dataset").integrate("flav",slice(5,6)).values()[()]),2)))
                    err_dn=np.sqrt(np.add(np.power(np.add(hflav[notdata].integrate("dataset").integrate("flav").values()[()],-1.*hflav_dn[notdata].integrate("dataset").integrate("flav").values()[()]),2),np.power(2*np.add(hflav[notdata].integrate("dataset").integrate("flav",slice(5,6)).values()[()],-1.*hflav_dn[notdata].integrate("dataset").integrate("flav",slice(5,6)).values()[()]),2)))
                data=hflav[datas].integrate("dataset").integrate("flav").values()[()]
                maximum=max(max((hflav[notdata].integrate("dataset").integrate("flav").values()[()]+err_up)),max(data))
            ratio_up=np.divide(err_up,data,out=np.zeros_like(err_up), where=data!=0)
            ratio_dn=np.divide(err_dn,data,out=np.zeros_like(err_up), where=data!=0)
            
            if arg.phase == "Wc":
                ax=plot.plot1d(hflav[notdata].sum("dataset").integrate("char"),overlay="flav",fill_opts={},error_opts=None,ax=ax,stack=True)
                plot.plot1d(hflav[notdata].sum("dataset",sumw2=True).sum("flav").sum("char"), ax=ax, density=False, clear=False,error_opts={'linestyle': 'none',
                    'markersize':0,
                    'elinewidth': 10,
                    'color':'tab:brown',
                    'alpha':.3,
                    'yerr':[err_up,err_dn]})
                plot.plot1d(hflav[notdata].sum("dataset").sum("flav").sum("char"),ax=ax,error_opts={'linestyle': 'none',
                    'markersize':0,
                    'elinewidth': 10,
                    'color':'tab:gray',
                    'alpha':.3}, clear=False,  density=False)
                
                plot.plot1d(hflav_nosf[notdata].sum("dataset").integrate("char").sum("flav"),error_opts= {'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none','color' :'tab:pink' , 'elinewidth': 1.5},ax=ax,clear=False)
                plot.plot1d(hflav[datas].sum("dataset").sum("flav").sum("char"),error_opts=data_err_opts,ax=ax,clear=False,  density=False)
                ax.legend(ncol=2,loc="upper right",handles=ax.get_legend_handles_labels()[0],labels=['b','c','pileup','udsg','SFs Unc.','stat. Unc.','w/o SFs','%s'%(arg.ext)],fontsize=13)
                ax.set_xlabel(None)

                ax.set_xticklabels(ax.get_xticklabels(), fontsize=0)
                rax = plot.plotratio(
                                                num=hflav[datas].sum("dataset").sum("flav").sum("char"),
                                                denom=hflav[notdata].sum("dataset").sum("flav").sum("char"),
                                                ax=rax,
                                                error_opts= {'linestyle': 'none','marker': '.', 'markersize': 0.,'color':'k'},
                                                denom_fill_opts={'yerr':[ratio_up,ratio_dn]},
                                                # denom_fill_opts={},
                                                guide_opts={},
                                                unc='num',
                                                clear=False)
                plot.plotratio(
                                                num=hflav[datas].sum("dataset").sum("flav").sum("char"),
                                                denom=hflav[notdata].sum("dataset").sum("flav").sum("char"),
                                                ax=rax,
                                                error_opts=data_err_opts,
                                                denom_fill_opts={},
                                                guide_opts={},
                                                unc='num',
                                                clear=False)
                plot.plotratio(
                                                num=hflav_nosf[datas].sum("dataset").sum("flav").sum("char"),
                                                denom=hflav_nosf[notdata].sum("dataset").sum("flav").sum("char"),
                                                ax=rax,
                                                error_opts= {'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none','color' :'tab:pink' , 'elinewidth': 1.5},
                                                denom_fill_opts={},
                                                guide_opts={},
                                                unc='num',
                                                clear=False)
            else :
                ax=plot.plot1d(hflav[notdata].sum("dataset"),overlay="flav",fill_opts={},error_opts=None,ax=ax,stack=True)

                plot.plot1d(hflav[notdata].sum("dataset",sumw2=True).sum("flav"), ax=ax, density=False, clear=False,error_opts={'linestyle': 'none',
                    'markersize':0,
                    'elinewidth': 10,
                    'color':'tab:brown',
                    'alpha':.3,
                    'yerr':[err_up,err_dn]})
                plot.plot1d(hflav[notdata].sum("dataset").sum("flav"),ax=ax,error_opts={'linestyle': 'none',
                    'markersize':0,
                    'elinewidth': 10,
                    'color':'tab:gray',
                    'alpha':.3}, clear=False,  density=False)
                
                plot.plot1d(hflav_nosf[notdata].sum("dataset").sum("flav"),error_opts= {'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none','color' :'tab:pink' , 'elinewidth': 1.5},ax=ax,clear=False)
                plot.plot1d(hflav[datas].sum("dataset").sum("flav"),error_opts=data_err_opts,ax=ax,clear=False,  density=False)
                ax.legend(ncol=2,loc="upper right",handles=ax.get_legend_handles_labels()[0],labels=['b','c','pileup','udsg','SFs Unc.','stat. Unc.','w/o SFs','%s'%(arg.ext)],fontsize=13)

                ax.set_xlabel(None)

                ax.set_xticklabels(ax.get_xticklabels(), fontsize=0)
                rax = plot.plotratio(
                                                num=hflav[datas].sum("dataset").sum("flav"),
                                                denom=hflav[notdata].sum("dataset").sum("flav"),
                                                ax=rax,
                                                error_opts= {'linestyle': 'none','marker': '.', 'markersize': 0.,'color':'k'},
                                                denom_fill_opts={'yerr':[ratio_up,ratio_dn]},
                                                # denom_fill_opts={},
                                                guide_opts={},
                                                unc='num',
                                                clear=False)
                plot.plotratio(
                                                num=hflav[datas].sum("dataset").sum("flav"),
                                                denom=hflav[notdata].sum("dataset").sum("flav"),
                                                ax=rax,
                                                error_opts=data_err_opts,
                                                denom_fill_opts={},
                                                guide_opts={},
                                                unc='num',
                                                clear=False)
                plot.plotratio(
                                                num=hflav_nosf[datas].sum("dataset").sum("flav"),
                                                denom=hflav_nosf[notdata].sum("dataset").sum("flav"),
                                                ax=rax,
                                                error_opts= {'linestyle': 'none','marker': 'o', 'markersize': 5.,'mfc': 'none','color' :'tab:pink' , 'elinewidth': 1.5},
                                                denom_fill_opts={},
                                                guide_opts={},
                                                unc='num',
                                                clear=False)                                    
        elif 'nbjet' in arg.discr_list :
            
            err_up=np.add(hflav_up[notdata].integrate("dataset").values()[()],-1.*hflav[notdata].integrate("dataset").values()[()])
            err_dn=np.add(hflav[notdata].integrate("dataset").values()[()],-1.*hflav_dn[notdata].integrate("dataset").values()[()])
            data=hflav[datas].integrate("dataset").values()[()]
            if not arg.log:maximum=max(max((hflav[notdata].integrate("dataset").values()[()]+err_up)),max(data))
            ratio_up=np.divide(err_up,data,out=np.zeros_like(err_up), where=data!=0)
            ratio_dn=np.divide(err_dn,data,out=np.zeros_like(err_up), where=data!=0)
            
            
            ax=plot.plot1d(hflav[notdata].sum("dataset"),fill_opts={},error_opts=None,ax=ax)
            plot.plot1d(hflav[notdata].sum("dataset",sumw2=True), ax=ax, density=False, clear=False,error_opts={'linestyle': 'none',
                'markersize':0,
                'elinewidth': 10,
                'color':'tab:brown',
                'alpha':.3,
                'yerr':[err_up,err_dn]})
            plot.plot1d(hflav[notdata].sum("dataset"),ax=ax,error_opts={'linestyle': 'none',
                'markersize':0,
                'elinewidth': 10,
                'color':'tab:gray',
                'alpha':.3}, clear=False,  density=False)
            
            plot.plot1d(hflav[datas].sum("dataset"),error_opts=data_err_opts,ax=ax,clear=False,  density=False)
            ax.legend(ncol=2,loc="upper right",handles=ax.get_legend_handles_labels()[0],labels=['MC','SFs Unc.','stat. Unc.','%s'%(arg.ext)],fontsize=13)
            ax.set_xlabel(None)

            rax = plot.plotratio(
                                            num=hflav[datas].sum("dataset"),
                                            denom=hflav[notdata].sum("dataset"),
                                            ax=rax,
                                            error_opts= {'linestyle': 'none','marker': '.', 'markersize': 0.,'color':'k'},
                                            denom_fill_opts={'yerr':[ratio_up,ratio_dn]},
                                            # denom_fill_opts={},
                                            guide_opts={},
                                            unc='num',
                                            clear=False)
            plot.plotratio(
                                            num=hflav[datas].sum("dataset"),
                                            denom=hflav[notdata].sum("dataset"),
                                            ax=rax,
                                            error_opts=data_err_opts,
                                            denom_fill_opts={},
                                            guide_opts={},
                                            unc='num',
                                            clear=False)
  
        
        elif 'mu' not in arg.discr_list and 'njet' is not arg.discr_list :
            data=hflav[datas].integrate("dataset").values()[()]
            if not arg.log:maximum=max(max((hflav[notdata].integrate("dataset").values()[()])),max(data))
          
            ax=plot.plot1d(hflav[notdata].sum("dataset"),error_opts=None,ax=ax)            
            plot.plot1d(hflav[datas].sum("dataset"),error_opts=data_err_opts,ax=ax,clear=False,  density=False)
            ax.legend(ncol=2,loc="upper right",handles=ax.get_legend_handles_labels()[0],labels=['MC','%s'%(arg.ext)],fontsize=13)
            ax.set_xlabel(None)
            # ax.set_xticklabels(ax.get_xticklabels(), fontsize=0)
            
            rax = plot.plotratio(
                                            num=hflav[datas].sum("dataset"),
                                            denom=hflav[notdata].sum("dataset"),
                                            ax=rax,
                                            error_opts=data_err_opts,
                                            denom_fill_opts={},
                                            guide_opts={},
                                            unc='num',
                                            clear=False)                                    
        else:
            if not arg.log:maximum=max(max((hflav[notdata].integrate("dataset").integrate("flav").values()[()])),max(data))
          
            ax=plot.plot1d(hflav[notdata].sum("dataset"),overlay="flav",fill_opts={},error_opts=None,ax=ax,stack=True)            
            plot.plot1d(hflav[datas].sum("dataset").sum("flav"),error_opts=data_err_opts,ax=ax,clear=False,  density=False)
            ax.legend(ncol=2,loc="upper right",handles=ax.get_legend_handles_labels()[0],labels=['b','c','pileup','udsg','%s'%(arg.ext)],fontsize=13)
            ax.set_xlabel(None)
            # ax.set_xticklabels(ax.get_xticklabels(), fontsize=0)
            
            rax = plot.plotratio(
                                            num=hflav[datas].sum("dataset").sum("flav"),
                                            denom=hflav[notdata].sum("dataset").sum("flav"),
                                            ax=rax,
                                            error_opts=data_err_opts,
                                            denom_fill_opts={},
                                            guide_opts={},
                                            unc='num',
                                            clear=False)
        if arg.log: 
            ax.set_ylim(1,maximum*100)
            ax.semilogy()
        else:ax.set_ylim(0,maximum*1.8)


        
        ax.set_ylabel("Events",fontsize=15)
        rax.set_ylabel('Data/MC',fontsize=15)
        if 'CvL' in discr :discrs=discr.replace('CvL','CvB')
        elif 'CvB' in discr :discrs=discr.replace('CvB','CvL')
        else :discrs=discr
        if arg.combine :rax.set_xlabel("%s"%(discrs),fontsize=15)
        else :rax.set_xlabel("%s[%d]"%(discrs,j),fontsize=15)
        rax.set_ylim(0.5,1.5)
        
        at = AnchoredText(input_txt+"\n"
                                #+ "inclusive pT, $\eta$"
                                , loc=2, prop=dict(size=15),frameon=False)
        ax.add_artist(at)
        scale=""
        if arg.norm:scale="_norm"
        name="all"
        if not arg.combine:name=str(j)
        if(arg.log):fig.savefig("plot/%s_unc_%s_inclusive%s_%s_%s.pdf" %(arg.phase, discrs, scale, arg.ext,name))
        else:fig.savefig("plot/%s_unc_lin_%s_inclusive%s_%s_%s.pdf" %(arg.phase, discrs, scale, arg.ext,name))
