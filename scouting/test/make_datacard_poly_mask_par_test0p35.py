# import ROOT in batch mode                                                                                                          
import os
import sys

#import sys                                                                                                                          
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

import numpy as np
from array import array

from ROOT import TH1F, TH1D, TH2D, TFile, TLorentzVector, TVector3, TChain, TProfile, TTree

# load FWLite C++ libraries                                                                                                          
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()

# load FWlite python libraries                                                                                                       
from DataFormats.FWLite import Handle, Events

mass = "0p35"

tree_muMC = ROOT.TChain('events')
tree_muMC.Add("/cms/routray/muon_ntuples_sl7/CMSSW_10_2_5/src/MuonAnalysis/Scouting/condor/flat_dimuon_tree_ggPhi_mass{}_ct5_new.root".format(mass))
tree_mudata = ROOT.TChain('events')
tree_mudata.Add("/cms/routray/muon_ntuples_sl7/CMSSW_10_2_5/src/MuonAnalysis/Scouting/condor/flatdimuontree_2017and2018_v1.root")

# lxybins = np.array([[0,0.1], [0.1,1], [1,3], [3,7], [7,11]])
 
lxybins = np.array([[0.0,0.1]])
#print lxybins[0,0], lxybins[0,1]

ggphipoly = open("ggphipolynew.csv", "a")
ggphipoly.write(" mass\tlxy bin\tpoly order\tchi2\tndof\tExpected 50.0%: r < \tExpected 16.0%: r < \tExpected 84.0%: r < \tExpected 2.5%: r < \tExpected 97.5%: r < \tObserved Limit\n")

# outfile = TFile("simple-shapes-TH1_mass0p35_ctau5_Lxy0_1.root", "recreate")
# outfile.cd()

binwidth = 0.001
ndecimal = 4

h3 = ROOT.TH1F("h3","h3", int(round(10/binwidth)), 0, 10)
tree_muMC.Draw('dimuon_mass>>h3','','')
fit_x = ROOT.TF1("fit_x", "gaus", 0, 10)
h3.Fit("fit_x", "R", "SAME")
mu =  h3.GetFunction("fit_x").GetParameter("Mean")
sig =  h3.GetFunction("fit_x").GetParameter("Sigma")

print float(binwidth/2.0)

xsigup = h3.GetBinCenter(h3.GetXaxis().FindBin(mu + 3*sig)) + float(binwidth/2.0)
xsigdown = h3.GetBinCenter(h3.GetXaxis().FindBin(mu - 3*sig)) - float(binwidth/2.0)

xfitup = h3.GetBinCenter(h3.GetXaxis().FindBin(mu + 10*sig)) + float(binwidth/2.0)
xfitdown = h3.GetBinCenter(h3.GetXaxis().FindBin(mu - 10*sig)) - float(binwidth/2.0)

bins = int(round((xfitup-xfitdown)/binwidth))

print xsigup, xsigdown, xfitup, xfitdown, bins
                                                                                   
signal2 = ROOT.TH1F("signal2", "Histogram of signal__x", int(bins), float(xfitdown), float(xfitup))
tree_muMC.Draw('dimuon_mass>>signal2','','') 

signal_sigmaUp = ROOT.TH1F("signal_sigmaUp", "Histogram of signal__x", int(bins), float(xfitdown), float(xfitup))
tree_muMC.Draw('dimuon_mass>>signal_sigmaUp','','')

signal_sigmaDown = ROOT.TH1F("signal_sigmaDown", "Histogram of signal__x", int(bins), float(xfitdown), float(xfitup))
tree_muMC.Draw('dimuon_mass>>signal_sigmaDown','','')
 

# gauss = ROOT.TF1("gauss", "({}/({}*sqrt(2*pi)))*exp((-0.5*pow(x-{},2))/(pow({},2)))".format(ns,sig,mu,sig), float(xfitdown), float(xfitup))
# for i in range(signal.GetNbinsX()):
#                 if(signal.GetBinCenter(i+1)>float(xfitdown)):
#                         print signal.GetBinCenter(i+1), gauss.Eval(signal.GetBinCenter(i+1)), gauss.Integral(signal.GetXaxis().GetBinLowEdge(i+1),signal.GetXaxis().GetBinUpEdge(i+1))                   
#                         signal.SetBinContent(i+1,gauss.Integral(signal.GetXaxis().GetBinLowEdge(i+1),signal.GetXaxis().GetBinUpEdge(i+1)))
#                         signal_sigmaUp.SetBinContent(i+1,gauss.Integral(signal.GetXaxis().GetBinLowEdge(i+1),signal.GetXaxis().GetBinUpEdge(i+1)))
#                         signal_sigmaDown.SetBinContent(i+1,gauss.Integral(signal.GetXaxis().GetBinLowEdge(i+1),signal.GetXaxis().GetBinUpEdge(i+1)))
# print signal.Integral()
#print len(lxybins)

h2 = []
h1 = []

for j in range(len(lxybins)):

        print "Looking at lxy bin----------",lxybins[j,0], "-", lxybins[j,1], "----------------"   


        signal1 = ROOT.TH1F("signal1", "signal1", int(bins), float(xfitdown), float(xfitup))                                                                                                               
        tree_muMC.Draw('dimuon_mass>>signal1',"lxy > {} && lxy < {} && muon1_trkiso < 0.1 && muon2_trkiso < 0.1 && dRmuon1jet > 0.3 && dRmuon2jet> 0.3 && abs(dphidimudv) < 0.02 && log(abs(deta\
mumu/dphimumu)) < 1.25 && distPixel > 0.05 && dimuon_mass > {} && dimuon_mass < {}".format(lxybins[j,0], lxybins[j,1], xfitdown, xfitup), '')                                                                                                                   
        print "signal events in this bin", signal1.Integral()

        accepatance = signal1.Integral()/100000

        # ns = signal1.Integral()
        ns = 100

        signal2.Scale(ns/signal2.Integral())
        signal_sigmaUp.Scale(ns/signal_sigmaUp.Integral())
        signal_sigmaDown.Scale(ns/signal_sigmaDown.Integral())

        h2.append(ROOT.TH1F("h2[{}]".format(j),"h2[{}]".format(j), int(bins), float(xfitdown), float(xfitup)))

        tree_mudata.Draw('dimuon_mass>>h2[{}]'.format(j),"lxy > {} && lxy < {} && muon1_trkiso < 0.1 && muon2_trkiso < 0.1 && dRmuon1jet > 0.3 && dRmuon2jet> 0.3 && abs(dphidimudv) < 0.02 && log(abs(detamumu/dphimumu)) < 1.25 && distPixel > 0.05 && dimuon_mass > {} && dimuon_mass < {}".format(lxybins[j,0], lxybins[j,1], xfitdown, xfitup),'')

        print float(xfitdown), float(xfitup)

        x_unmasked = []
        y_unmasked = []
        print h2[j].GetNbinsX()
        for i in range(h2[j].GetNbinsX()):
                if(h2[j].GetBinCenter(i+1)>float(xfitdown)): 
                        #print h2[j].GetBinCenter(i+1)
                        x_unmasked.append(round(h2[j].GetBinCenter(i+1),ndecimal))
                        y_unmasked.append(h2[j].GetBinContent(i+1))

        print "x_unmasked", x_unmasked
        print "y_unmasked", y_unmasked
        print "y_unmasked_error", np.sqrt(y_unmasked)

        data = ROOT.TH1F("data","Histogram of data_obs__x", int(bins), float(xfitdown), float(xfitup))
        for i in range(h2[j].GetNbinsX()):
                # print i, x_unmasked[i],y_unmasked[i]
                data.SetBinContent(i+1,y_unmasked[i])

        h1.append(ROOT.TH1F("h1[{}]".format(j),"h1[{}]".format(j), int(bins), float(xfitdown), float(xfitup)))

        tree_mudata.Draw('dimuon_mass>>h1[{}]'.format(j),"lxy > {} && lxy < {} && muon1_trkiso < 0.1 && muon2_trkiso < 0.1 && dRmuon1jet > 0.3 && dRmuon2jet> 0.3 && abs(dphidimudv) < 0.02 && log(abs(detamumu/dphimumu)) < 1.25 && distPixel > 0.05 && dimuon_mass > {} && dimuon_mass < {} && (dimuon_mass < {} || dimuon_mass > {})".format(lxybins[j,0], lxybins[j,1], xfitdown, xfitup, xsigdown, xsigup),'')

        x_masked = []
        y_masked = []
        print h1[j].GetNbinsX()
        for i in range(h1[j].GetNbinsX()):
                # if(h1[j].GetBinCenter(i+1)>float(xfitdown) and h1[j].GetBinContent(i+1)>0):
                #         x_masked.append(round(h1[j].GetBinCenter(i+1),3))
                #         y_masked.append(h1[j].GetBinContent(i+1))

                if(h1[j].GetBinCenter(i+1)>float(xfitdown) and h1[j].GetBinCenter(i+1)<float(xsigdown)):
                        x_masked.append(round(h1[j].GetBinCenter(i+1),ndecimal))
                        y_masked.append(h1[j].GetBinContent(i+1))

                # if(h1[j].GetBinCenter(i+1)>=float(xsigdown) and h1[j].GetBinCenter(i+1)<=float(xsigup) and h1[j].GetBinContent(i+1)>0):
                #         x_masked.append(round(h1[j].GetBinCenter(i+1),3))
                #         y_masked.append(h1[j].GetBinContent(i+1))
 
                if(h1[j].GetBinCenter(i+1)>float(xsigup) and h1[j].GetBinCenter(i+1)<float(xfitup)):
                        x_masked.append(round(h1[j].GetBinCenter(i+1),ndecimal))
                        y_masked.append(h1[j].GetBinContent(i+1))


        print "x_masked", x_masked
        print "y_masked", y_masked
        print "y_masked_error", np.sqrt(y_masked)


        x_sigdata = [] 
        y_sigdata = []
        for i in range(h2[j].GetNbinsX()):
                if(h2[j].GetBinCenter(i+1)>float(xsigdown) and h2[j].GetBinCenter(i+1)<float(xsigup)):
                        x_sigdata.append(round(h2[j].GetBinCenter(i+1),ndecimal))
                        y_sigdata.append(h2[j].GetBinContent(i+1))

############################################################################################################

        x = ROOT.RooRealVar("x","x",float(xfitdown),float(xfitup))

        x.setRange("R1",float(xfitdown),float(xsigdown))
        x.setRange("R2",float(xsigup),float(xfitup))

        l = ROOT.RooArgList(x)

        data_obs = ROOT.RooDataHist("data_obs", "data_obs", l, data)

	mean = ROOT.RooRealVar("mean","Mean of Gaussian",mu)
        sigma = ROOT.RooRealVar("sigma","Width of Gaussian",sig)
        signal = ROOT.RooGaussian("signal","signal",x,mean,sigma)

	nS = ns
	sig_norm = ROOT.RooRealVar("sig_norm","sig_norm",nS,0,10*nS)

	p0 =  ROOT.RooRealVar("p0","p0",1,-10000000,10000000)
        p1 =  ROOT.RooRealVar("p1","p1",1,-10000000,10000000)
        p2 =  ROOT.RooRealVar("p2","p2",1,-10000000,10000000)
        p3 =  ROOT.RooRealVar("p3","p3",1,-10000000,10000000)
        # p4 =  ROOT.RooRealVar("p4","p4",1,-10000000,10000000)
        # p5 =  ROOT.RooRealVar("p5","p5",1,-10000000,10000000)
        # p6 =  ROOT.RooRealVar("p6","p6",1,-10000000,10000000)


        # p0 =  ROOT.RooRealVar("p0","p0",1000,-10000000,10000000)                                                                                                                                         
        # p1 =  ROOT.RooRealVar("p1","p1",-1000,-10000000,10000000)                                                                                                                                        
        # p2 =  ROOT.RooRealVar("p2","p2",-1000,-10000000,10000000)                                                                                                                                        
        # p3 =  ROOT.RooRealVar("p3","p3",-1000,-10000000,10000000)                                                                                                                                        
        # p4 =  ROOT.RooRealVar("p4","p4",-1000,-10000000,10000000)                                                                                                                                        
        # p5 =  ROOT.RooRealVar("p5","p5",20000,-10000000,10000000)                                                                                                                                        
        # p6 =  ROOT.RooRealVar("p6","p6",100000,-10000000,10000000)                                                                                                                                       



        # background = ROOT.RooGenericPdf("background","background","p0 + p1*pow(x,1) + p2*pow(x,2) + p3*pow(x,3) + p4*pow(x,4) + p5*pow(x,5) + p6*pow(x,6)",ROOT.RooArgSet(x,p0,p1,p2,p3,p4,p5,p6))
        # background = ROOT.RooPolynomial("background","background", x, ROOT.RooArgList(p0,p1,p2,p3,p4))
	# background = ROOT.RooChebychev("background","background", x, ROOT.RooArgList(p0,p1,p2,p3,p4,p5,p6))
	background = ROOT.RooBernstein("background","background", x, ROOT.RooArgList(p0,p1,p2,p3))
  
	nB = data.Integral()
	background_norm = ROOT.RooRealVar("background_norm","background_norm",nB,0.9*nB,1.1*nB)
	# background.fitTo(data_obs, ROOT.RooFit.Range("R1,R2"))
 	
        model = ROOT.RooAddPdf("model","model",ROOT.RooArgList(background),ROOT.RooArgList(background_norm))
	result = ROOT.RooFitResult(model.fitTo(data_obs, ROOT.RooFit.Range("R1,R2"), ROOT.RooFit.Save(ROOT.kTRUE), ROOT.RooFit.Minimizer("Minuit2","Migrad")))
        model.fitTo(data_obs,ROOT.RooFit.Range("R1,R2"))

	# chi2 = ROOT.RooChi2Var("chi2","chi2",model,data,ROOT.RooFit.Range("R1,R2"))
	# m = ROOT.RooMinuit(chi2) 
	# m.migrad() 
	# m.improve()
	# m.hesse() 
        # result = ROOT.RooFitResult(m.save())

        bkg_component = ROOT.RooArgSet(background)

        xframe = x.frame(ROOT.RooFit.Title("Data Fit"))
        data_obs.plotOn(xframe, ROOT.RooFit.Name("data"))
        # model.plotOn(xframe,ROOT.RooFit.LineColor(3),ROOT.RooFit.Components(bkg_component),ROOT.RooFit.Name("bkg"), ROOT.RooFit.LineStyle(2), ROOT.RooFit.Range("R1,R2"), ROOT.RooFit.NormRange("R1,R2"))
        model.plotOn(xframe,ROOT.RooFit.LineColor(3),ROOT.RooFit.Components(bkg_component),ROOT.RooFit.Name("bkg"), ROOT.RooFit.LineStyle(2), ROOT.RooFit.Range("Full"), ROOT.RooFit.NormRange("Full"))
        print xframe.chiSquare(7)
	print result.minNll()
	model.plotOn(xframe,ROOT.RooFit.VisualizeError(result,1,ROOT.kFALSE), ROOT.RooFit.DrawOption("F"), ROOT.RooFit.Components(bkg_component), ROOT.RooFit.Name("errorband"), ROOT.RooFit.FillColor(ROOT.kOrange), ROOT.RooFit.Range("Full"), ROOT.RooFit.NormRange("Full")) 
        data_obs.plotOn(xframe, ROOT.RooFit.Name("data"))
        model.plotOn(xframe,ROOT.RooFit.LineColor(3),ROOT.RooFit.Components(bkg_component), ROOT.RooFit.Name("bkg"), ROOT.RooFit.LineStyle(2), ROOT.RooFit.Range("Full"), ROOT.RooFit.NormRange("Full"))

        c1 = ROOT.TCanvas("c1","c1") 
        c1.cd()

        leg1 = ROOT.TLegend(1.75,5100,1.85,5300)
        leg1.SetLineColor(0)
        leg1.SetFillColor(0)
        leg1.AddEntry(xframe.findObject("data"), "Data", "pe")
	leg1.AddEntry(xframe.findObject("bkg"), "#color[2]{Poly Fit}", "l")                                                                                          
	leg1.AddEntry(xframe.findObject("errorband"), "Fit Error", "")
        leg1.SetTextFont(42)
        leg1.SetBorderSize(0)
        leg1.Draw()

        # xframe.GetYaxis().SetRangeUser(14000,16000)
        xframe.Draw()

        c1.BuildLegend()
        c1.Draw()
        c1.SaveAs("testmass{}_lxy{}_{}_poly.png".format(mass, lxybins[j,0], lxybins[j,1]))

        #######################################

        c2 = ROOT.TCanvas("c2","c2")
        pad1 = ROOT.TPad("pad1", "The pad 80% of the height", 0.0, 0.2, 1.0, 1.0, 0)
        pad2 = ROOT.TPad("pad2", "The pad 20% of the height", 0.0, 0.0, 1.0, 0.2, 0)
        c2.cd()

        pad1.Draw()
        pad2.Draw()

        pad1.cd()
        pad1.SetTickx()
        pad1.SetTicky()
        pad1.SetBottomMargin(0.01)

        # ROOT.gStyle.SetEndErrorSize(0)
        xframe2 = x.frame(ROOT.RooFit.Title(" "))
        data_obs.plotOn(xframe2, ROOT.RooFit.Name("data"))
        model.plotOn(xframe2,ROOT.RooFit.LineColor(3),ROOT.RooFit.Components(bkg_component), ROOT.RooFit.LineStyle(2), ROOT.RooFit.Name("bkg"), ROOT.RooFit.Range("Full"), ROOT.RooFit.NormRange("Full"))
	model.plotOn(xframe2,ROOT.RooFit.VisualizeError(result,1,ROOT.kFALSE), ROOT.RooFit.DrawOption("F"), ROOT.RooFit.Components(bkg_component), ROOT.RooFit.Name("errorband"), ROOT.RooFit.FillColor(ROOT.kOrange), ROOT.RooFit.Range("Full"), ROOT.RooFit.NormRange("Full"))
        data_obs.plotOn(xframe2, ROOT.RooFit.Name("data"))
        model.plotOn(xframe2,ROOT.RooFit.LineColor(2),ROOT.RooFit.Components(bkg_component), ROOT.RooFit.Name("bkg"), ROOT.RooFit.LineStyle(1), ROOT.RooFit.Range("Full"), ROOT.RooFit.NormRange("Full"))

        xframe2.Draw()
        xframe2.GetYaxis().SetRangeUser(-20,120)                                                                                                                                                       
        xframe2.GetYaxis().SetTitle("Events/ 0.01 GeV")
        xframe2.GetYaxis().SetTitleSize(0.05)
        xframe2.GetYaxis().SetLabelSize(0.045)
        xframe2.GetYaxis().SetTitleOffset(0.95)

	box = ROOT.TBox(float(xsigdown),0,float(xsigup),500) 
	box.SetFillColorAlpha(7,0.35)
	box.SetFillStyle(1001)
	box.Draw()

        leg1 = ROOT.TLegend(0.1,0.6,0.4,0.9)
        leg1.SetLineColor(0)
        leg1.SetFillColor(0)
        leg1.SetFillStyle(0)
        leg1.AddEntry(xframe2.findObject("data"), "Data [{} < lxy < {}]".format(lxybins[j,0], lxybins[j,1]), "pe") 
        leg1.AddEntry(xframe2.findObject("bkg"), "#color[2]{Poly Fit}", "l")
        leg1.AddEntry(xframe2.findObject("errorband"), "Fit Error", "f")
        leg1.SetTextFont(42)
        leg1.SetBorderSize(0)
        leg1.Draw()

        pull = ROOT.RooHist(xframe2.pullHist("data","bkg"))
        pull.SetFillColor(ROOT.kRed)
        pull.SetLineWidth(0)

        xframe3 = x.frame(ROOT.RooFit.Title(" "))
        
        xframe3.addPlotable(pull,"B X")
        xframe3.GetXaxis().SetLabelSize(0.17)
        xframe3.GetYaxis().SetLabelSize(0.15)
        xframe3.GetXaxis().SetTitleSize(0.21)
        xframe3.GetYaxis().SetTitleSize(0.15)
        xframe3.GetXaxis().SetTitleOffset(0.85)
        xframe3.GetYaxis().SetTitleOffset(0.28)
        xframe3.GetXaxis().SetTitle("Dimuon Mass [GeV]")
        xframe3.GetYaxis().SetTitle("#scale[1.3]{#frac{data - fit}{#sigma_{data}}}")
        # xframe3.GetYaxis().SetTitle("Pull")
        xframe3.GetYaxis().SetLabelSize(0.15)

        pad2.cd()
        pad2.SetTickx()
        pad2.SetTicky()
        #pad2.SetGridy()                                                                                                                                                                                 
        pad2.SetTopMargin(0.0)
        pad2.SetBottomMargin(0.4)
        xframe3.Draw()

        # c2.BuildLegend()
        c2.Draw()
        c2.SaveAs("mass{}_lxy{}_{}_poly.png".format(mass, lxybins[j,0], lxybins[j,1]))
        
        ###########################################################################################
        
	myWS = ROOT.RooWorkspace("myWS", "workspace")
        getattr(myWS,'import')(data_obs)
        getattr(myWS,'import')(background)
        getattr(myWS,'import')(signal)
	# getattr(myWS, 'import')(model,ROOT.RooCmdArg())                                                                                                                                                  
        # getattr(myWS, 'import')(model)   



        
        datacard = open("simple-shapes-TH1_mass{}_ctau5_Lxy{}_{}_poly.txt".format(mass, lxybins[j,0],lxybins[j,1]), "w")

        datacard.write("imax 1  number of channels\n")
        datacard.write("jmax 1  number of backgrounds\n")
        datacard.write("kmax *  number of nuisance parameters (sources of systematical uncertainties)\n")
        datacard.write("------------------------------------\n")
        datacard.write("shapes * * simple-shapes-TH1_mass{}_ctau5_Lxy{}_{}_poly.root myWS:$PROCESS\n".format(mass, lxybins[j,0],lxybins[j,1]))
        datacard.write("------------------------------------\n")
        datacard.write("bin bin1\n")
        datacard.write("observation -1\n")
        datacard.write("------------------------------------\n")
        datacard.write("bin bin1 bin1\n")
        datacard.write("process signal background\n")
        datacard.write("process 0 1\n")
        datacard.write("rate {} {}\n".format(nS, nB))
        datacard.write("------------------------------------\n")
        datacard.write("lumi lnN 1.025 1.0\n")
        # datacard.write("bgnorm lnN 1.00 {}\n".format(background_alphaUp.Integral()/background.Integral()))
        # datacard.write("bgnorm lnN 1.00 1.2\n") 
        # datacard.write("alpha shapeN2 - 1 uncertainty on background shape and normalization\n")
        datacard.close() 

	myWS.writeToFile("simple-shapes-TH1_mass{}_ctau5_Lxy{}_{}_poly.root".format(mass, lxybins[j,0],lxybins[j,1]))
        myWS.Print()
        print "RooWorkspace made"
        ROOT.gDirectory.Add(myWS)
        
        os.system('combine -M  AsymptoticLimits --rAbsAcc=0.0001 --rRelAcc=0.001 simple-shapes-TH1_mass{}_ctau5_Lxy{}_{}_poly.txt > com.out'.format(mass, lxybins[j,0],lxybins[j,1]))

        os.system('cat com.out')
        com_out = open('com.out','r')

        for line in com_out:
                if line[:15] == 'Observed Limit:':
                        coml_obs = float(line[19:])

                elif line[:15] == 'Expected  2.5%:':
                        coml_2sd = float(line[19:])

                elif line[:15] == 'Expected 16.0%:':
                        coml_1sd = float(line[19:])
        
                elif line[:15] == 'Expected 50.0%:':
                        coml_exp = float(line[19:])

                elif line[:15] == 'Expected 84.0%:':
                        coml_1su = float(line[19:])

                elif line[:15] == 'Expected 97.5%:':
                        coml_2su = float(line[19:])


        ggphipoly.write(" {}\t{} - {}\t6\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(mass, lxybins[j,0], lxybins[j,1], xframe.chiSquare(7), int(0), coml_exp, coml_1sd, coml_1su, coml_2sd, coml_2su, coml_obs))

        
