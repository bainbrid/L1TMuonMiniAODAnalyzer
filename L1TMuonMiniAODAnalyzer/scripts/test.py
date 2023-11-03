################################################################################
# Imports
from __future__ import print_function
import builtins
import future
from future.utils import raise_with_traceback
from argparse import ArgumentParser
import pandas as pd
import numpy as np
import awkward as ak
import ROOT as r
import uproot
import math
from setTDRStyle import setTDRStyle

################################################################################
# Utility

def phi_phasewrap(phi) :
    return (phi + np.pi) % (2 * np.pi) - np.pi

def deltaR(obj1,obj2) :
    deta = obj1.eta - obj2.eta
    dphi = phi_phasewrap(obj1.phi - obj2.phi)
    return np.sqrt( deta*deta + dphi*dphi )

################################################################################
# Configure

parser = ArgumentParser()
parser.add_argument('--verbosity',default=0,type=int)
parser.add_argument('--nevents',default=0,type=int)
args = parser.parse_args()
print("Command line args:",vars(args))

columns = [
    'event',
    'run',
    'ls',
    'bx',
    'nPV',
    #
    'gen_idx',
    'gen_acc',
    'gen_pt',
    'gen_eta',
    'gen_phi',
    'gen_charge',
    'gen_dxy',
    'gen_pdgid',
    'gen_muon_dr',
    'gen_muon_match',
    'gen_muon_idx',
    #
    'muon_pt',
    'muon_eta',
    'muon_phi',
    'muon_charge',
    'muon_etaAtSt1',
    'muon_etaAtSt2',
    'muon_phiAtSt1',
    'muon_phiAtSt2',
    'muon_dz',
    'muon_dzError',
    'muon_dxy',
    'muon_dxyError',
    'muon_3dIP',
    'muon_3dIPError',
    'muon_PassTightID',
    'muon_PassLooseID',
    'muon_isSAMuon',
    'muon_trg_dr',
    'muon_trg_match',
    'muon_trg_idx',
    #
    'trg_pt',
    'trg_pt_dxy',
    'trg_eta',
    'trg_phi',
    'trg_charge',
    'trg_etaAtVtx',
    'trg_phiAtVtx',
    'trg_qual',
    'trg_dxy',
    'trg_tfIdx',
    'trg_bx',
    ]
columns = list(set(columns))

################################################################################
# Get DataFrame
def get_data(files,columns) :
    if args.verbosity > 0 : print('Getting files:\n', '\n'.join(files))
    dfs = [ uproot.open(i)['L1TMuonMiniAODAnalyzerFlat/tree'].arrays(columns,library="pd") for i in files ]
    if args.verbosity > 0 : print('Extracted branches: ',columns)
    df = pd.concat(dfs)

    if args.verbosity > 0 : print('Available branches: ',df.keys())
    if args.nevents > 0 : 
        print("Considering only first {:.0f} events ...".format(args.nevents))
        df = df.head(args.nevents)
    return df

files = [
    #"./L1TMuonNtuple_run3_mc.root",
    #"./L1TMuonNtuple_run3_mc_effs.root",
    "./L1TMuonNtuple_run3_mc_purities.root",
    ]
data = get_data(files,columns)
if args.verbosity > 2 :
    print(data.columns)
    print(data.dtypes)

################################################################################
# Engineered features

################################################################################
# Print summary of DataFrame
if args.verbosity > 1 :
   pd.options.display.max_columns=None
   pd.options.display.width=None
   print(data.describe().T)
   print(data.info())
   #pretty=lambda df:tabulate(df,headers='keys',tablefmt='psql') # 'html'
   #print(pretty(data.describe().T)))

################################################################################
# Plotting, utility

def plot1d(histos,data=None):
    for signal,name,title,binning,vals in histos:
        if vals is None and data is not None : vals = data[name]
        c = r.TCanvas()
        prefix = 'signal' if signal else 'bkgd'
        his = r.TH1F(prefix+'_'+name,'',*binning)
        for val in ak.ravel(vals) : his.Fill(val if val < binning[2] else binning[2]-1.e-6)
        his.GetXaxis().SetTitle(title)
        his.SetLineWidth(2)
        his.SetLineColor(r.kGreen+3 if signal else r.kRed)
        his.Draw('')
        c.SaveAs('plots/'+prefix+'_'+name+'.pdf')

def plot2d(histos):
    for signal,name,title,binning,xvals,yvals in histos:
        c = r.TCanvas()
        prefix = 'signal' if signal else 'bkgd'
        his = r.TH2F(prefix+'_'+name,'',*binning)
        for xval,yval in zip(np.ravel(xvals),np.ravel(yvals)) :
            his.Fill(
                xval if xval < binning[2] else binning[2]-1.e-6,
                yval if yval < binning[5] else binning[5]-1.e-6,
                )
        his.GetXaxis().SetTitle(title.split(':')[1])
        his.GetYaxis().SetTitle(title.split(':')[2])
        his.Draw('colz text')
        r.gStyle.SetOptStat(0)
        c.SaveAs('plots/'+prefix+'_'+name+'.pdf')

################################################################################
# Plotting, simple 1D

def plot1d_gen(data,purity=False):
    histos = [
        (True,'gen_idx','GEN particle index',(501,-0.5,500.5),data.gen_idx),
        (True,'gen_acc','GEN muon in acceptance',(2,-0.5,0.5),data.gen_acc),
        (True,'gen_pt','GEN muon p_{T} [GeV]',(100,0.,10.),data.gen_pt),
        (True,'gen_eta','GEN muon #eta',(100,-5.,5.),data.gen_eta),
        (True,'gen_phi','GEN muon #phi',(64,-3.2,3.2),data.gen_phi),
        (True,'gen_charge','GEN muon charge',(3,-1.5,1.5),data.gen_charge),
        (True,'gen_dxy','GEN muon #Deltaxy',(100,0.,1.),data.gen_dxy),
        (True,'gen_pdgid','GEN particle PDG ID',(41,-20.5,20.5),data.gen_pdgid),
        (True,'gen_muon_dr','#DeltaR(GEN,muon)',(100,0.,0.1),data.gen_muon_dr),
        (True,'gen_muon_match','GEN muon matched to muon',(2,-0.5,1.5),data.gen_muon_match),
        (True,'gen_muon_idx','GEN matched to muon idx',(22,-1.5,20.5),data.gen_muon_idx),
        ] if purity==False else [
        (True,'gen_idx','GEN particle index',(501,-0.5,500.5),data.gen_idx[data.gen_muon_match>0]),
        (True,'gen_acc','GEN muon in acceptance',(2,-0.5,0.5),data.gen_acc[data.gen_muon_match>0]),
        (True,'gen_pt','GEN muon p_{T} [GeV]',(100,0.,10.),data.gen_pt[data.gen_muon_match>0]),
        (True,'gen_eta','GEN muon #eta',(100,-5.,5.),data.gen_eta[data.gen_muon_match>0]),
        (True,'gen_phi','GEN muon #phi',(64,-3.2,3.2),data.gen_phi[data.gen_muon_match>0]),
        (True,'gen_charge','GEN muon charge',(3,-1.5,1.5),data.gen_charge[data.gen_muon_match>0]),
        (True,'gen_dxy','GEN muon #Deltaxy',(100,0.,1.),data.gen_dxy[data.gen_muon_match>0]),
        (True,'gen_pdgid','GEN particle PDG ID',(41,-20.5,20.5),data.gen_pdgid[data.gen_muon_match>0]),
        (True,'gen_muon_dr','#DeltaR(GEN,muon)',(100,0.,0.1),data.gen_muon_dr[data.gen_muon_match>0]),
        (True,'gen_muon_match','GEN muon matched to muon',(2,-0.5,1.5),data.gen_muon_match[data.gen_muon_match>0]),
        (True,'gen_muon_idx','GEN matched to muon idx',(22,-1.5,20.5),data.gen_muon_idx[data.gen_muon_match>0]),
        ]
    plot1d(histos)

def plot1d_muon(data,purity=False):    
    histos = [
        (True,'muon_pt','Muon p_{T} [GeV]',(100,0.,10.),data.muon_pt[data.gen_muon_match>0]),
        (True,'muon_eta','Muon #eta',(100,-5.,5.),data.muon_eta[data.gen_muon_match>0]),
        (True,'muon_phi','Muon #phi',(64,-3.2,3.2),data.muon_phi[data.gen_muon_match>0]),
        (True,'muon_charge','Muon charge',(3,-1.5,1.5),data.muon_charge[data.gen_muon_match>0]),
        (True,'muon_dz','Muon #Deltaz',(100,0.,10.),data.muon_dz[data.gen_muon_match>0]),
        (True,'muon_dzError','Muon #Deltaz error',(100,0.,1.),data.muon_dzError[data.gen_muon_match>0]),
        (True,'muon_dxy','Muon #Deltaxy',(100,0.,1.),data.muon_dxy[data.gen_muon_match>0]),
        (True,'muon_dxyError','Muon #Deltaxy error',(100,0.,0.1),data.muon_dxyError[data.gen_muon_match>0]),
        (True,'muon_3dIP','Muon 3D IP',(100,0.,10.),data.muon_3dIP[data.gen_muon_match>0]),
        (True,'muon_3dIPError','Muon 3D IP error',(100,0.,0.1),data.muon_3dIPError[data.gen_muon_match>0]),
        (True,'muon_PassTightID','Muon Tight ID',(2,-0.5,1.5),data.muon_PassTightID[data.gen_muon_match>0]),
        (True,'muon_PassLooseID','Muon Loose ID',(2,-0.5,1.5),data.muon_PassLooseID[data.gen_muon_match>0]),
        (True,'muon_isSAMuon','Muon standalone',(2,-0.5,1.5),data.muon_isSAMuon[data.gen_muon_match>0]),
        (True,'muon_trg_dr','#DeltaR(TRG@Vtx,muon)',(100,0.,1.),data.muon_trg_dr[data.gen_muon_match>0]),
        (True,'muon_trg_match','TRG muon matched to muon',(2,-0.5,1.5),data.muon_trg_match[data.gen_muon_match>0]),
        (True,'muon_trg_idx','TRG matched to muon idx',(22,-1.5,20.5),data.muon_trg_idx[data.gen_muon_match>0]),
        ] if purity==False else [
        (True,'muon_pt','Muon p_{T} [GeV]',(100,0.,10.),data.muon_pt[data.muon_trg_match>0]),
        (True,'muon_eta','Muon #eta',(100,-5.,5.),data.muon_eta[data.muon_trg_match>0]),
        (True,'muon_phi','Muon #phi',(64,-3.2,3.2),data.muon_phi[data.muon_trg_match>0]),
        (True,'muon_charge','Muon charge',(3,-1.5,1.5),data.muon_charge[data.muon_trg_match>0]),
        (True,'muon_dz','Muon #Deltaz',(100,0.,10.),data.muon_dz[data.muon_trg_match>0]),
        (True,'muon_dzError','Muon #Deltaz error',(100,0.,1.),data.muon_dzError[data.muon_trg_match>0]),
        (True,'muon_dxy','Muon #Deltaxy',(100,0.,1.),data.muon_dxy[data.muon_trg_match>0]),
        (True,'muon_dxyError','Muon #Deltaxy error',(100,0.,0.1),data.muon_dxyError[data.muon_trg_match>0]),
        (True,'muon_3dIP','Muon 3D IP',(100,0.,10.),data.muon_3dIP[data.muon_trg_match>0]),
        (True,'muon_3dIPError','Muon 3D IP error',(100,0.,0.1),data.muon_3dIPError[data.muon_trg_match>0]),
        (True,'muon_PassTightID','Muon Tight ID',(2,-0.5,1.5),data.muon_PassTightID[data.muon_trg_match>0]),
        (True,'muon_PassLooseID','Muon Loose ID',(2,-0.5,1.5),data.muon_PassLooseID[data.muon_trg_match>0]),
        (True,'muon_isSAMuon','Muon standalone',(2,-0.5,1.5),data.muon_isSAMuon[data.muon_trg_match>0]),
        (True,'muon_trg_dr','#DeltaR(TRG@Vtx,muon)',(100,0.,1.),data.muon_trg_dr[data.muon_trg_match>0]),
        (True,'muon_trg_match','TRG muon matched to muon',(2,-0.5,1.5),data.muon_trg_match[data.muon_trg_match>0]),
        (True,'muon_trg_idx','TRG matched to muon idx',(22,-1.5,20.5),data.muon_trg_idx[data.muon_trg_match>0]),
        ]
    plot1d(histos)
    
def plot1d_trg(data,purity=False):
    histos = [
        (True,'trg_pt','L1 muon p_{T} [GeV]',(41,-0.25,20.25),data.trg_pt[data.muon_trg_match>0]),
        (True,'trg_pt_dxy','L1 muon displaced p_{T} [GeV]',(21,-0.5,20.5),data.trg_pt_dxy[data.muon_trg_match>0]),
        (True,'trg_eta','L1 muon #eta',(100,-5.,5.),data.trg_eta[data.muon_trg_match>0]),
        (True,'trg_phi','L1 muon #phi',(64,-3.2,3.2),data.trg_phi[data.muon_trg_match>0]),
        (True,'trg_charge','L1 muon charge',(3,-1.5,1.5),data.trg_charge[data.muon_trg_match>0]),
        (True,'trg_etaAtVtx','L1 muon #eta at vertex',(100,-5.,5.),data.trg_etaAtVtx[data.muon_trg_match>0]),
        (True,'trg_phiAtVtx','L1 muon #phi at vertex',(64,-3.2,3.2),data.trg_phiAtVtx[data.muon_trg_match>0]),
        (True,'trg_dxy','L1 muon #Deltaxy',(4,-0.5,3.5),data.trg_dxy[data.muon_trg_match>0]),
        (True,'trg_tfIdx','L1 muon TF index',(128,-0.5,127.5),data.trg_tfIdx[data.muon_trg_match>0]),
        (True,'trg_qual','L1 muon quality',(16,-0.5,15.5),data.trg_qual[data.muon_trg_match>0]),
        (True,'trg_bx','L1 muon bunch crossing',(21,-10.5,10.5),data.trg_bx[data.muon_trg_match>0]),
        ] if purity==False else [
        (True,'trg_pt','L1 muon p_{T} [GeV]',(41,-0.25,20.25),data.trg_pt),
        (True,'trg_pt_dxy','L1 muon displaced p_{T} [GeV]',(21,-0.5,20.5),data.trg_pt_dxy),
        (True,'trg_eta','L1 muon #eta',(100,-5.,5.),data.trg_eta),
        (True,'trg_phi','L1 muon #phi',(64,-3.2,3.2),data.trg_phi),
        (True,'trg_charge','L1 muon charge',(3,-1.5,1.5),data.trg_charge),
        (True,'trg_etaAtVtx','L1 muon #eta at vertex',(100,-5.,5.),data.trg_etaAtVtx),
        (True,'trg_phiAtVtx','L1 muon #phi at vertex',(64,-3.2,3.2),data.trg_phiAtVtx),
        (True,'trg_dxy','L1 muon #Deltaxy',(4,-0.5,3.5),data.trg_dxy),
        (True,'trg_tfIdx','L1 muon TF index',(128,-0.5,127.5),data.trg_tfIdx),
        (True,'trg_qual','L1 muon quality',(16,-0.5,15.5),data.trg_qual),
        (True,'trg_bx','L1 muon bunch crossing',(21,-10.5,10.5),data.trg_bx),
        ]
    plot1d(histos)
    
################################################################################
#
def plot_eff(data,option="pT",cumu=False):

    setTDRStyle()
    c1 = r.TCanvas()
    r.gStyle.SetOptStat(1111)
    
    xmin = 0.14
    ymax = 0.88

    legend = r.TLegend(xmin,ymax-0.035*9,xmin+0.15,ymax)
    legend.SetTextSize(0.04)
    legend.SetTextFont(42)
    legend.SetTextAngle(0)
    legend.SetTextColor(r.kBlack)
    legend.SetTextSize(0.03)
    legend.SetTextAlign(12)
    legend.SetBorderSize(0)

    trg_nbins = 9
    trg_min = 7.
    trg_max = 16.

    pt_nbins = 20
    pt_min = 0.
    pt_max = 20.

    effs = {}
    for idx in range(trg_nbins) :
        trg_pt_thr = trg_min + idx

        name = f"trg_{idx}"
        title = "L1 p_{T}"f" > {trg_pt_thr:.0f} GeV"

        presel = (data.gen_acc>0) & (data.gen_muon_match>0) #& (data.muon_pt>trg_min)
        passed = presel & (data.muon_trg_match>0)

        # Which trigger pT variable to use?
        if option == "pT"      : passed &= (data.trg_pt>trg_pt_thr)
        elif option == "Upt"   : passed &= (data.trg_pt_dxy>trg_pt_thr)
        elif option == "pTUpt" : passed &= (data.trg_pt>trg_pt_thr) & (data.trg_pt_dxy>trg_pt_thr)
        else                   : passed = False
        values = data.muon_pt[presel]

        numer = r.TH1F("numer","",pt_nbins,pt_min,pt_max)
        denom = r.TH1F("denom","",pt_nbins,pt_min,pt_max)

        for val in values:         denom.Fill(min(val,pt_max-1.e-6))
        for val in values[passed]: numer.Fill(min(val,pt_max-1.e-6))

        if cumu==False : eff = r.TEfficiency(numer,denom)
        else : eff = r.TEfficiency(numer.GetCumulative(False),denom.GetCumulative(False))

        effs[name] = eff
        eff.Draw("") if idx==0 else eff.Draw("same") 
        c1.Pad().Update()
        if idx==0:
            eff.Draw("")
            temp = ";Muon p_{T} threshold [GeV];Cumulative efficiency" if cumu else ";Muon p_{T} [GeV];Differential efficiency"
            eff.SetTitle(temp)
            eff.GetPaintedGraph().GetXaxis().SetLimits(pt_min,pt_max)
            eff.GetPaintedGraph().GetXaxis().SetNdivisions(505)
            eff.GetPaintedGraph().GetYaxis().SetRangeUser(0.,1.)
            eff.GetPaintedGraph().GetYaxis().SetNdivisions(505)
        eff.SetMarkerColor(1+idx)
        eff.SetMarkerStyle(20)
        eff.SetMarkerSize(1.5)
        eff.SetLineColor(1+idx)
        legend.AddEntry(eff,title,"pe")
    
        #eff = effs[name]
        #c1 = r.TCanvas()
        #his_denom = eff.GetCopyTotalHisto().GetCumulative(False).Draw()
        #c1.SaveAs(f"plots/temp_denom_trg_pt{trg_pt_thr:.0f}_{option}.pdf")
        #his_numer = eff.GetCopyPassedHisto().GetCumulative(False).Draw()
        #c1.SaveAs(f"plots/temp_numer_trg_pt{trg_pt_thr:.0f}_{option}.pdf")
        #eff.Draw()
        #c1.SaveAs(f"plots/temp_eff_{option}.pdf")

    r.gStyle.SetOptStat(0)
    legend.Draw()
    c1.Update()
    suffix = "_cumu" if cumu else ""
    c1.SaveAs(f"plots/eff_{option}{suffix}.pdf")

################################################################################
#
def cumu(data,option="pT",fixed=True,purity=False):
    trg_nbins = 9
    trg_min = 7.
    trg_max = 16.
    his_numer = r.TH1F("numer","",trg_nbins,trg_min-0.5,trg_max-0.5)
    his_denom = r.TH1F("denom","",trg_nbins,trg_min-0.5,trg_max-0.5)
    for idx in range(trg_nbins):
        trg_pt_thr = trg_min + idx
        thr = 4. if fixed else trg_pt_thr+1.
        denom = None
        numer = None
        if purity==True:
            if option == "pT"      : denom = (data.trg_pt>trg_pt_thr)
            elif option == "Upt"   : denom = (data.trg_pt_dxy>trg_pt_thr)
            elif option == "pTUpt" : denom = (data.trg_pt>trg_pt_thr) & (data.trg_pt_dxy>trg_pt_thr)
            else                   : denom = False
            numer = denom & (data.gen_acc>0) & (data.gen_muon_match>0) & (data.muon_pt>thr) & (data.muon_trg_match>0)
        elif purity==False:
            denom = (data.gen_acc>0) & (data.gen_muon_match>0) & (data.muon_pt>thr)
            numer = denom & (data.muon_trg_match>0)
            if option == "pT"      : numer &= (data.trg_pt>trg_pt_thr)
            elif option == "Upt"   : numer &= (data.trg_pt_dxy>trg_pt_thr)
            elif option == "pTUpt" : numer &= (data.trg_pt>trg_pt_thr) & (data.trg_pt_dxy>trg_pt_thr)
            else                   : numer = False
        his_numer.SetBinContent(idx+1,np.sum(numer))
        his_denom.SetBinContent(idx+1,np.sum(denom))
    his_numer.Sumw2()
    his_denom.Sumw2()
    eff = r.TEfficiency(his_numer.GetCumulative(False),his_denom.GetCumulative(False))
    gr = eff.GetPaintedGraph()
    return eff.Clone()

################################################################################
#
def plot_roc(data,fixed=True,quantity="eff"):
    setTDRStyle()
    #r.gStyle.SetOptStat(1111)

    trg_nbins = 9
    trg_min = 7.
    trg_max = 16.
    
    rates = {}
    rates["pT"] = [(66516,902.171),(43977,733.564),(30419,610.098),(21915,517.842),
                   (16458,448.756),(12811,395.934),(9911 ,348.25 ),(8076 ,314.356),(6461 ,281.168),]
    rates["Upt"] = [(26577,570.268),(22478,524.45 ),(19113,483.604),(19590,489.604),
                    (14500,421.219),(12787,395.555),(11453,374.358),(10217,353.584),(9361 ,338.439),]
    rates["pTUpt"] = [(19003,482.208),(14622,422.993),(11123,368.919),(8504 ,322.583),
                      (6877 ,290.08 ),(5323 ,255.208),(4209 ,226.949),(3414 ,204.386),(2802 ,185.168),]

    effs = {}
    effs["pT"] = cumu(data,option="pT",fixed=fixed,purity=(quantity=="purity"))
    effs["Upt"] = cumu(data,option="Upt",fixed=fixed,purity=(quantity=="purity"))
    effs["pTUpt"] = cumu(data,option="pTUpt",fixed=fixed,purity=(quantity=="purity"))

    c1 = r.TCanvas()
    c1.SetLogx()
    c1.SetLogy()
    c1.SetGrid()

    xmin = 0.14
    ymax = 0.88
    legend = r.TLegend(xmin,ymax-0.06*3,xmin+0.34,ymax)
    legend.SetTextSize(0.04)
    legend.SetTextFont(42)
    legend.SetTextAngle(0)
    legend.SetTextColor(r.kBlack)
    legend.SetTextSize(0.04)
    legend.SetTextAlign(12)
    legend.SetBorderSize(0)
    
    graphs = {}
    pretty = ["p_{T} (7#minus15 GeV)","Up_{T} (7#minus15 GeV)","p_{T} & Up_{T} (7#minus15 GeV)"]
    colours = [r.kRed,r.kOrange+1,r.kGreen+1]
    for idx,option in enumerate(["pT","Upt","pTUpt"]) :

        obj = effs[option] if quantity=="eff" or quantity=="purity" else effs[option].GetCopyPassedHisto()

        #c1 = r.TCanvas()
        #his_denom = obj.GetCopyTotalHisto().Draw()
        #his_denom.Draw()
        #c1.SaveAs(f"plots/temp_denom_{option}.pdf")
        #his_numer = obj.GetCopyPassedHisto().Draw()
        #his_numer.Draw()
        #c1.SaveAs(f"plots/temp_numer_{option}.pdf")
        #obj.Draw()
        #c1.SaveAs(f"plots/temp_eff_{option}.pdf")

        _x = None
        _ex = None
        if quantity=="eff" or quantity=="purity" :
            _x  = np.array([ float(obj.GetEfficiency(idx+1)) for idx in range(trg_nbins) ])
            _ex = np.array([ float(obj.GetEfficiencyErrorUp(idx+1)) for idx in range(trg_nbins) ])
        else: 
            _x  = np.array([ float(obj.GetBinContent(x+1)) for x in range(obj.GetNbinsX()) ])
            #ex = np.array([ float(obj.GetBinError(x+1)) for x in range(obj.GetNbinsX()) ])
            _ex = np.array([ np.sqrt(x) for x in _x ])

        _y  = np.array([ float(x) for x in list(list(zip(*rates[option]))[0])])
        _ey = np.array([ float(x) for x in list(list(zip(*rates[option]))[1])])
        #print(_x,_ex,_y,_ey)

        graphs[option] = r.TGraphErrors(trg_nbins,_x,_y,_ex,_ey)
        graphs[option].SetLineStyle(1)
        graphs[option].SetLineWidth(2)
        graphs[option].SetLineColor(colours[idx])
        graphs[option].SetMarkerStyle(20)
        graphs[option].SetMarkerColor(colours[idx])
        if idx==0:
            r.gStyle.SetTitleFont(42)
            r.gStyle.SetTitleAlign(23)
            r.gStyle.SetTitleX(0.5)
            temp = "Muon p_{T} threshold = 4 GeV" if fixed else "Muon p_{T} threshold = L1 threshold + 1 GeV"
            graphs[option].SetTitle(temp)
            if quantity=="eff":
                graphs[option].GetXaxis().SetTitle("L1 trigger efficiency")
                limits = (0.01,0.2) if fixed else (0.1,1.)
                graphs[option].GetXaxis().SetLimits(*limits)
            elif quantity=="purity":
                graphs[option].GetXaxis().SetTitle("L1 trigger purity")
                limits = (0.01,1.) if fixed else (0.01,1.)
                graphs[option].GetXaxis().SetLimits(*limits)
            elif quantity=="number":
                graphs[option].GetXaxis().SetTitle("Number of offline muons matched to trigger")
                graphs[option].GetXaxis().SetLimits(500.,100000.)
            else:
                graphs[option].GetXaxis().SetTitle("Unknown quantity!")
            
            graphs[option].GetYaxis().SetTitle("L1 trigger rate [Hz]")
            graphs[option].GetYaxis().SetRangeUser(1000.,100000.)
            graphs[option].Draw("ALPZ")
        else:
            graphs[option].Draw("LPZsame")
        legend.AddEntry(graphs[option],pretty[idx],"LPZ")

    legend.Draw("same")
    suffix = "fixed" if fixed else "dynamic"
    c1.SaveAs(f"plots/roc_{suffix}_{quantity}.pdf")
    
################################################################################
# Steer plotting

plot1d_gen(data,purity=True)
plot1d_muon(data,purity=True)
plot1d_trg(data,purity=True)

plot_eff(data,option="pT",cumu=False)
plot_eff(data,option="Upt",cumu=False)
plot_eff(data,option="pTUpt",cumu=False)
plot_eff(data,option="pT",cumu=True)
plot_eff(data,option="Upt",cumu=True)
plot_eff(data,option="pTUpt",cumu=True)

#plot_roc(data,fixed=True,quantity="eff")
#plot_roc(data,fixed=True,quantity="number")
plot_roc(data,fixed=True,quantity="purity")
#plot_roc(data,fixed=False,quantity="eff")
#plot_roc(data,fixed=False,quantity="number")
plot_roc(data,fixed=False,quantity="purity")
    
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
