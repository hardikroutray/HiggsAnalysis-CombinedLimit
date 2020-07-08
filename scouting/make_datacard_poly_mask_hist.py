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

mass = "2"

tree_muMC = ROOT.TChain('events')
tree_muMC.Add("/uscms/home/hroutray/nobackup/CMSSW_10_2_13/src/flat_dimuon_tree_ggPhi_mass{}_ct5_new.root".format(mass))
tree_mudata = ROOT.TChain('events')
tree_mudata.Add("/uscms/home/hroutray/nobackup/CMSSW_10_2_13/src/flatdimuontree_2017and2018_v1.root")

# lxybins = np.array([[0,0.1], [0.1,1], [1,3], [3,7], [7,11]])
 
lxybins = np.array([[0.0,0.1]])
#print lxybins[0,0], lxybins[0,1]

ggphipoly = open("ggphipolynew.csv", "a")
ggphipoly.write(" mass\tlxy bin\tpoly order\tchi2\tndof\tsqrt(sum signalpull**2)\tsqrt(sum fiterr**2)\tExpected 50.0%: r < \tExpected 16.0%: r < \tExpected 84.0%: r < \tExpected 2.5%: r < \tExpected 97.5%: r < \tObserved Limit\n")

# outfile = TFile("simple-shapes-TH1_mass0p35_ctau5_Lxy0_1.root", "recreate")
# outfile.cd()

binwidth = 0.01
ndecimal = 3

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

# print bins
# print int(round(bins))
# print (int(round(bins)) % 2.0)

# if int(round(bins)) % 2.0 == 0:
#         bins = bins
# else:
#         bins = bins+1

print xsigup, xsigdown, xfitup, xfitdown, bins

# signal1 = ROOT.TH1F("signal1", "signal1", int(bins), float(xfitdown), float(xfitup))
# tree_muMC.Draw('dimuon_mass>>signal1','','')
                                                                                   
signal = ROOT.TH1F("signal", "Histogram of signal__x", int(bins), float(xfitdown), float(xfitup))
tree_muMC.Draw('dimuon_mass>>signal','','') 

signal_sigmaUp = ROOT.TH1F("signal_sigmaUp", "Histogram of signal__x", int(bins), float(xfitdown), float(xfitup))
tree_muMC.Draw('dimuon_mass>>signal_sigmaUp','','')

signal_sigmaDown = ROOT.TH1F("signal_sigmaDown", "Histogram of signal__x", int(bins), float(xfitdown), float(xfitup))
tree_muMC.Draw('dimuon_mass>>signal_sigmaDown','','')
 
# print signal1.Integral()
# ns = 1000

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

        signal.Scale(ns/signal.Integral())
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

        data_obs = ROOT.TH1F("data_obs","Histogram of data_obs__x", int(bins), float(xfitdown), float(xfitup))
        for i in range(h2[j].GetNbinsX()):
                # print i, x_unmasked[i],y_unmasked[i]
                data_obs.SetBinContent(i+1,y_unmasked[i])

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


##########################################################################################################

	'''
	pol = ROOT.TF1("pol","[0] + [1]*x + [2]*x**2 + [3]*x**3 + [4]*x**4 + [5]*x**5 + [6]*x**6", float(xfitdown) , float(xfitup))
        pol.SetParameters(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0)

        r = h2[j].Fit("pol", "SL", "SAME")
        f = h2[j].GetFunction("pol")
 
        print f.GetChisquare()
        print f.GetNDF()

        for p in range(pol.GetNpar()):
                print pol.GetParameter(p), pol.GetParameter(p) - pol.GetParError(p), pol.GetParameter(p) + pol.GetParError(p), "\n",

	'''

############################################################################################################

        x = ROOT.RooRealVar("x","x",float(xfitdown),float(xfitup))

        x.setRange("R1",float(xfitdown),float(xsigdown))
        x.setRange("R2",float(xsigup),float(xfitup))


        l = ROOT.RooArgList(x)

        data = ROOT.RooDataHist("data", "data", l, data_obs)

        # mean = ROOT.RooRealVar("mean","Mean of Gaussian",mu)
        # sigma = ROOT.RooRealVar("sigma","Width of Gaussian",sig)
        # sign = ROOT.RooGaussian("sign","sign",x,mean,sigma)

	# p0 =  ROOT.RooRealVar("p0","p0",504087,-1.68909e+06,2.69726e+06)                                                                                                                                 
        # p1 =  ROOT.RooRealVar("p1","p1",-228890,-1.34128e+06,883500)                                                                                                                                     
        # p2 =  ROOT.RooRealVar("p2","p2",-35819.9,-443971,372331)                                                                                                                                         
        # p3 =  ROOT.RooRealVar("p3","p3",14137.2,-158933,187207)                                                                                                                                          
        # p4 =  ROOT.RooRealVar("p4","p4",12272.3,-64110.8,88655.3)                                                                                                                                        
	# p5 =  ROOT.RooRealVar("p5","p5",9013.25,-23359.3,41385.8)                                                                                                                                        
	# p6 =  ROOT.RooRealVar("p6","p6",-4882.69,-16952.9,7187.51)       


	# p0 =  ROOT.RooRealVar("p0","p0",pol.GetParameter(0), pol.GetParameter(0) - 1000*pol.GetParError(0), pol.GetParameter(0) + 1000*pol.GetParError(0))
        # p1 =  ROOT.RooRealVar("p1","p1",pol.GetParameter(1), pol.GetParameter(1) - 1000*pol.GetParError(1), pol.GetParameter(1) + 1000*pol.GetParError(1))
        # p2 =  ROOT.RooRealVar("p2","p2",pol.GetParameter(2), pol.GetParameter(2) - 1000*pol.GetParError(2), pol.GetParameter(2) + 1000*pol.GetParError(2))
        # p3 =  ROOT.RooRealVar("p3","p3",pol.GetParameter(3), pol.GetParameter(3) - 1000*pol.GetParError(3), pol.GetParameter(3) + 1000*pol.GetParError(3))
        # p4 =  ROOT.RooRealVar("p4","p4",pol.GetParameter(4), pol.GetParameter(4) - 1000*pol.GetParError(4), pol.GetParameter(4) + 1000*pol.GetParError(4))
        # p5 =  ROOT.RooRealVar("p5","p5",pol.GetParameter(5), pol.GetParameter(5) - 1000*pol.GetParError(5), pol.GetParameter(5) + 1000*pol.GetParError(5))
        # p6 =  ROOT.RooRealVar("p6","p6",pol.GetParameter(6), pol.GetParameter(6) - 1000*pol.GetParError(6), pol.GetParameter(6) + 1000*pol.GetParError(6))

        # p0 =  ROOT.RooRealVar("p0","p0",1000,-10000000,10000000)                                                                                                                                  
        # p1 =  ROOT.RooRealVar("p1","p1",-1000,-10000000,10000000)                                                                                                                                      
        # p2 =  ROOT.RooRealVar("p2","p2",-1000,-10000000,10000000)                                                                                                                                        
        # p3 =  ROOT.RooRealVar("p3","p3",-1000,-10000000,10000000)                                                                                                                                        
        # p4 =  ROOT.RooRealVar("p4","p4",-1000,-10000000,10000000)                                                                                                                                        
        # p5 =  ROOT.RooRealVar("p5","p5",20000,-10000000,10000000)                                                                                                                                        
        # p6 =  ROOT.RooRealVar("p6","p6",100000,-10000000,10000000)                                                                                                                                       
 	# p7 = ROOT.RooRealVar("p7","p7",1000,-10000000,10000000)


	p0 =  ROOT.RooRealVar("p0","p0",1,-10000000,10000000)
        p1 =  ROOT.RooRealVar("p1","p1",1,-10000000,10000000)
        p2 =  ROOT.RooRealVar("p2","p2",1,-10000000,10000000)
        p3 =  ROOT.RooRealVar("p3","p3",1,-10000000,10000000)
        p4 =  ROOT.RooRealVar("p4","p4",1,-10000000,10000000)
        p5 =  ROOT.RooRealVar("p5","p5",1,-10000000,10000000)
        p6 =  ROOT.RooRealVar("p6","p6",1,-10000000,10000000)

	# p0 =  ROOT.RooRealVar("p0","p0",100,-10000000,10000000)
        # p1 =  ROOT.RooRealVar("p1","p1",1,-10000000,10000000)
        # p2 =  ROOT.RooRealVar("p2","p2",1,-10000000,10000000)
        # p3 =  ROOT.RooRealVar("p3","p3",1,-10000000,10000000)
        # p4 =  ROOT.RooRealVar("p4","p4",100,-10000000,10000000)
        # p5 =  ROOT.RooRealVar("p5","p5",1,-10000000,10000000)
        # p6 =  ROOT.RooRealVar("p6","p6",100,-10000000,10000000)
        # p7 = ROOT.RooRealVar("p7","p7",100,-10000000,10000000)




        # background = ROOT.RooGenericPdf("background","background","p0 + p1*pow(x,1) + p2*pow(x,2) + p3*pow(x,3) + p4*pow(x,4) + p5*pow(x,5) + p6*pow(x,6)",ROOT.RooArgSet(x,p0,p1,p2,p3,p4,p5,p6))
        background = ROOT.RooPolynomial("background","background", x, ROOT.RooArgList(p0,p1,p2,p3,p4,p5,p6))
	# background = ROOT.RooChebychev("background","background", x, ROOT.RooArgList(p0,p1,p2,p3,p4,p5,p6))
	# background = ROOT.RooBernstein("background","background", x, ROOT.RooArgList(p0,p1,p2,p3,p4,p5,p6))


        # background.fitTo(data)    
        # data_obs2 = ROOT.RooDataHist("data_obs2","data_obs2",x,data_obs) 

        # nB = data_obs.Integral()
        # nS = ns
        # background_norm = ROOT.RooRealVar("background_norm","background_norm",nB,0.9*nB,1.1*nB)
        # sig_norm = ROOT.RooRealVar("sig_norm","sig_norm",nS,0,10*nS)

	nS = ns
	nB = data_obs.Integral()
	background_norm = ROOT.RooRealVar("background_norm","background_norm",nB,0.9*nB,1.1*nB)

        # model = ROOT.RooAddPdf("model","model",ROOT.RooArgList(background,sign),ROOT.RooArgList(background_norm, sig_norm))
        model = ROOT.RooAddPdf("model","model",ROOT.RooArgList(background),ROOT.RooArgList(background_norm))


        # result1 = ROOT.RooFitResult(background.fitTo(data, ROOT.RooFit.Save(ROOT.kTRUE), ROOT.RooFit.Strategy(2), ROOT.RooFit.SumW2Error(ROOT.kFALSE), ROOT.RooFit.Minimizer("Miniuit2")))
        # result1 = ROOT.RooFitResult(background.fitTo(data, ROOT.RooFit.Save(ROOT.kTRUE), ROOT.RooFit.Strategy(2), ROOT.RooFit.SumW2Error(ROOT.kFALSE)))
        # background.fitTo(data)
        # background.fitTo(data, ROOT.RooFit.Range("R1,R2"))
 
        # result = ROOT.RooFitResult(model.fitTo(data, ROOT.RooFit.Save(ROOT.kTRUE), ROOT.RooFit.Strategy(2), ROOT.RooFit.SumW2Error(ROOT.kFALSE), ROOT.RooFit.Minimizer("Miniuit2")))
        # result = ROOT.RooFitResult(model.fitTo(data, ROOT.RooFit.Save(ROOT.kTRUE), ROOT.RooFit.Strategy(2), ROOT.RooFit.SumW2Error(ROOT.kFALSE)))
        # result = ROOT.RooFitResult(model.fitTo(data, ROOT.RooFit.Range("R1,R2"), ROOT.RooFit.Save(ROOT.kTRUE), ROOT.RooFit.Minimizer("Minuit2","Migrad")))  
	
	result = ROOT.RooFitResult(model.fitTo(data, ROOT.RooFit.Range("R1,R2"), ROOT.RooFit.Save(ROOT.kTRUE), ROOT.RooFit.Minimizer("Minuit2","Migrad")))
        # model.fitTo(data)
        model.fitTo(data,ROOT.RooFit.Range("R1,R2"))

	'''
	print p0.getVal(), p0.getVal() - 10*p0.getError(), p0.getVal() + 10*p0.getError(), "\n"
	print p1.getVal(), p1.getVal() - 10*p1.getError(), p1.getVal() + 10*p1.getError(), "\n"
	print p2.getVal(), p2.getVal() - 10*p2.getError(), p2.getVal() + 10*p2.getError(), "\n"
	print p3.getVal(), p3.getVal() - 10*p3.getError(), p3.getVal() + 10*p3.getError(), "\n"
	print p4.getVal(), p4.getVal() - 10*p4.getError(), p4.getVal() + 10*p4.getError(), "\n"
	print p5.getVal(), p5.getVal() - 10*p5.getError(), p5.getVal() + 10*p5.getError(), "\n"
	print p6.getVal(), p6.getVal() - 10*p6.getError(), p6.getVal() + 10*p6.getError(), "\n"

	
	p0 =  ROOT.RooRealVar("p0","p0", p0.getVal(), p0.getVal() - 3*p0.getError(), p0.getVal() + 3*p0.getError())
        p1 =  ROOT.RooRealVar("p1","p1",p1.getVal(), p1.getVal() - 3*p1.getError(), p1.getVal() + 3*p1.getError())
        p2 =  ROOT.RooRealVar("p2","p2",p2.getVal(), p2.getVal() - 3*p2.getError(), p2.getVal() + 3*p2.getError())
        p3 =  ROOT.RooRealVar("p3","p3",p3.getVal(), p3.getVal() - 3*p3.getError(), p3.getVal() + 3*p3.getError())
        p4 =  ROOT.RooRealVar("p4","p4",p4.getVal(), p4.getVal() - 3*p4.getError(), p4.getVal() + 3*p4.getError())
        p5 =  ROOT.RooRealVar("p5","p5",p5.getVal(), p5.getVal() - 3*p5.getError(), p5.getVal() + 3*p5.getError())
        p6 =  ROOT.RooRealVar("p6","p6",p6.getVal(), p6.getVal() - 3*p6.getError(), p6.getVal() + 3*p6.getError())

	background = ROOT.RooPolynomial("background","background", x, ROOT.RooArgList(p0,p1,p2,p3,p4,p5,p6))
	nB = data_obs.Integral()
        nS = ns
        background_norm = ROOT.RooRealVar("background_norm","background_norm",nB,0.9*nB,1.1*nB)
        sig_norm = ROOT.RooRealVar("sig_norm","sig_norm",nS,0,10*nS)
	model = ROOT.RooAddPdf("model","model",ROOT.RooArgList(background),ROOT.RooArgList(background_norm))

	result = ROOT.RooFitResult(model.fitTo(data, ROOT.RooFit.Range("R1,R2"), ROOT.RooFit.Save(ROOT.kTRUE), ROOT.RooFit.Minimizer("Miniuit", "Migrad")))
	model.fitTo(data,ROOT.RooFit.Range("R1,R2"))
	'''

	# chi2 = ROOT.RooChi2Var("chi2","chi2",model,data,ROOT.RooFit.Range("R1,R2"))
	# m = ROOT.RooMinuit(chi2) 
	# m.migrad() 
	# m.improve()
	# m.hesse() 
        # result = ROOT.RooFitResult(m.save())

        bkg_component = ROOT.RooArgSet(background)
        # sig_component = ROOT.RooArgSet(sign)

        xframe = x.frame(ROOT.RooFit.Title("data generated from composite"))
        data.plotOn(xframe, ROOT.RooFit.Name("data"))
        # model.plotOn(xframe,ROOT.RooFit.LineColor(2))
        # model.plotOn(xframe,ROOT.RooFit.LineColor(3),ROOT.RooFit.Components(bkg_component),ROOT.RooFit.Name("bkg"), ROOT.RooFit.LineStyle(2), ROOT.RooFit.Range("R1,R2"), ROOT.RooFit.NormRange("R1,R2"))
        model.plotOn(xframe,ROOT.RooFit.LineColor(3),ROOT.RooFit.Components(bkg_component),ROOT.RooFit.Name("bkg"), ROOT.RooFit.LineStyle(2), ROOT.RooFit.Range("Full"), ROOT.RooFit.NormRange("Full"))

        # model.plotOn(xframe,ROOT.RooFit.LineColor(6),ROOT.RooFit.Name("bkg"), ROOT.RooFit.LineStyle(1)) 
        # background.plotOn(xframe,ROOT.RooFit.LineColor(3),ROOT.RooFit.Name("bkg"), ROOT.RooFit.LineStyle(2))
        # model.plotOn(xframe,ROOT.RooFit.LineColor(6),ROOT.RooFit.Components(sig_component),ROOT.RooFit.LineStyle(2))
        # func = ROOT.TF1(background.asTF(l, ROOT.RooArgList(p0,p1,p2,p3,p4,p5,p6), ROOT.RooArgSet(x)))
        print xframe.chiSquare()
	print result.minNll()
        # model.plotOn(xframe,ROOT.RooFit.VisualizeError(result,1,ROOT.kFALSE), ROOT.RooFit.DrawOption("L"), ROOT.RooFit.Components(bkg_component), ROOT.RooFit.Name("errorband"), ROOT.RooFit.FillColor(ROOT.kOrange))
        # model.plotOn(xframe,ROOT.RooFit.Components(bkg_component), ROOT.RooFit.VisualizeError(result,ROOT.RooArgSet(sig_norm),2), ROOT.RooFit.Name("errorband"), ROOT.RooFit.FillColor(ROOT.kOrange))
        # model.plotOn(xframe,ROOT.RooFit.Components(bkg_component), ROOT.RooFit.VisualizeError(result,1), ROOT.RooFit.Name("errorband"), ROOT.RooFit.FillColor(ROOT.kOrange), ROOT.RooFit.Range("R1,R2"), ROOT.RooFit.NormRange("R1,R2"))

        # model.plotOn(xframe, ROOT.RooFit.Components(bkg_component), ROOT.RooFit.VisualizeError(result,1), ROOT.RooFit.Name("errorband"), ROOT.RooFit.FillColor(ROOT.kOrange), ROOT.RooFit.Range("Full"), ROOT.RooFit.NormRange("Full"))

	model.plotOn(xframe,ROOT.RooFit.VisualizeError(result,1,ROOT.kFALSE), ROOT.RooFit.DrawOption("F"), ROOT.RooFit.Components(bkg_component), ROOT.RooFit.Name("errorband"), ROOT.RooFit.FillColor(ROOT.kOrange), ROOT.RooFit.Range("Full"), ROOT.RooFit.NormRange("Full")) 

        # background.plotOn(xframe,ROOT.RooFit.VisualizeError(result1,1), ROOT.RooFit.Name("errorband"), ROOT.RooFit.FillColor(ROOT.kOrange), ROOT.RooFit.Range("Full"))
 
        data.plotOn(xframe, ROOT.RooFit.Name("data"))
        # model.plotOn(xframe,ROOT.RooFit.LineColor(6),ROOT.RooFit.Name("bkg"), ROOT.RooFit.LineStyle(1))
        # model.plotOn(xframe,ROOT.RooFit.LineColor(3),ROOT.RooFit.Components(bkg_component), ROOT.RooFit.Name("bkg"), ROOT.RooFit.LineStyle(2), ROOT.RooFit.Range("R1,R2"), ROOT.RooFit.NormRange("R1,R2"))
        model.plotOn(xframe,ROOT.RooFit.LineColor(3),ROOT.RooFit.Components(bkg_component), ROOT.RooFit.Name("bkg"), ROOT.RooFit.LineStyle(2), ROOT.RooFit.Range("Full"), ROOT.RooFit.NormRange("Full"))
        # model.plotOnWithErrorBand(xframe,ROOT.RooFit.LineColor(3),ROOT.RooFit.Components(bkg_component), ROOT.RooFit.Name("bkg"), ROOT.RooFit.LineStyle(2), ROOT.RooFit.Range("Full"), ROOT.RooFit.NormRange("Full"))

        # background.plotOn(xframe,ROOT.RooFit.LineColor(3),ROOT.RooFit.Name("bkg"), ROOT.RooFit.LineStyle(2), ROOT.RooFit.Range("Full"))
                
        # scale_factor = data_obs.Integral()*data_obs.GetBinWidth(1)
        scale_factor = background_norm.getVal()*data_obs.GetBinWidth(1)
        # scaled_background = ROOT.RooFormulaVar("scaled_background", "@0 * @1", ROOT.RooArgList(background, ))
        # formula = ROOT.TString("{} * background".format(scale_factor))
        scaled_background = ROOT.RooFormulaVar("scaled_background", "{} * background".format(scale_factor), ROOT.RooArgList(background))
        scaledfunc = ROOT.TF1(scaled_background.asTF(l, ROOT.RooArgList(p0,p1,p2,p3,p4,p5,p6), ROOT.RooArgSet(x)))

        xframe.Print("v")

        fitcurve = ROOT.RooCurve(xframe.getCurve("bkg"))
        fitgraph = ROOT.TGraph(fitcurve.GetN())

        for i in range(fitcurve.GetN()):
                fitgraph.SetPoint(i, fitcurve.GetX()[i], fitcurve.GetY()[i])

        central = ROOT.RooCurve(xframe.getCurve("errorband"))
        curve = ROOT.RooCurve(xframe.getCurve("errorband"))
        upBound = ROOT.TGraph(central.GetN())
        loBound = ROOT.TGraph(central.GetN())


        for i in range(curve.GetN()):
                if i < central.GetN():
                        upBound.SetPoint(i, curve.GetX()[i], curve.GetY()[i])
                else:
                        loBound.SetPoint(i, curve.GetX()[i], curve.GetY()[i])

                
        y_pred = []
        sigmaerr = []

        print h2[j].GetNbinsX()
        for i in range(h2[j].GetNbinsX()):                                                                                         
                if(h2[j].GetBinCenter(i+1)>float(xfitdown)):                                                                       
                        #print h2[j].GetBinCenter(i+1)                 
                        dataCount  = h2[j].GetBinContent(i+1)                                                                      
                        fitValue = scaledfunc.Eval(h2[j].GetBinCenter(i+1))
                        upperValue = upBound.Eval(h2[j].GetBinCenter(i+1))
                        lowerValue = loBound.Eval(h2[j].GetBinCenter(i+1))
                        fitValue2 = fitgraph.Eval(h2[j].GetBinCenter(i+1))
                   
                        y_pred.append(fitValue2)
                        sigmaerr.append(upperValue-fitValue2)

                        print h2[j].GetBinCenter(i+1), "----", dataCount, "----", fitValue, "----", upperValue, "-----", fitValue2, "----", fitValue2 - (upperValue - fitValue2), "----", background_norm.getVal(), "----" , upperValue-fitValue2, "-----", lowerValue


        # scaledfunc = ROOT.TF1(scaled_background.asTF(l))

        # binnedbkg = background.generateBinned(ROOT.RooArgSet(x),background_norm.getVal())
        # binnedbkg.plotOn(xframe,ROOT.RooFit.LineColor(2))

        # var = ROOT.RooRealVar(model.getObservables(ROOT.RooArgSet(x)).first())
        # # var = ROOT.RooRealVar(model.getObservables(ROOT.RooArgSet(x)))

        # for i in range(h2[j].GetNbinsX()):
        #         if(h2[j].GetBinCenter(i+1)>float(xfitdown)):
        #                 #print h2[j].GetBinCenter(i+1)                                                                              
        #                 dataCount  = h2[j].GetBinContent(i+1)
        #                 var.setVal(h2[j].GetBinCenter(i+1))
                        
        #                 print var.getVal(), "----", dataCount, model.getVal(ROOT.RooArgSet(var))*background_norm.getVal()*binwidth, " ---- ", model.getVal(), "----"
        
        c1 = ROOT.TCanvas("c1","c1") 
        c1.cd()


        leg1 = ROOT.TLegend(1.75,5100,1.85,5300)
        leg1.SetLineColor(0)
        leg1.SetFillColor(0)
        leg1.AddEntry(xframe.findObject("data"), "Data", "pe")
 
        # TString chi2_bbf4_signal="#color[3]{chi^2 bbf4} = " + std::to_string(xframe2->chiSquare("background","data_obs2"));                                                                             
        # TString chi2_bbf4="#scale[1]{#color[2]{pol6}} #scale[1.0]{#color[2]{ #chi^{2}}} = " + std::to_string(signalchi);//+"/"+std::to_string(n_bins-3);                                                          
        # lxybin = ROOT.TString("#scale[1]{#color[4]{0 < lxy < 0.1}}")
        # leg1.AddEntry(xframe.findObject("errorband"), "Fit Error", "")
        leg1.AddEntry(xframe.findObject("bkg"), "#color[2]{Pol6 + Gaus Fit}", "l")                                                                                          
        leg1.AddEntry(xframe.findObject("errorband"), "Fit Error", "")
        leg1.SetTextFont(42)
        leg1.SetBorderSize(0)
        leg1.Draw()

        # xframe.GetYaxis().SetRangeUser(14000,16000)
        xframe.Draw()
        # func.Draw("SAME")
        # func.SetLineColor(2)

        # scaledfunc.Draw("SAME")
        # scaledfunc.SetLineColor(2)



	# fafill = ROOT.TF1("fafill","masked region",float(xsigdown),float(xsigup))



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
        data.plotOn(xframe2, ROOT.RooFit.Name("data"))
        model.plotOn(xframe2,ROOT.RooFit.LineColor(3),ROOT.RooFit.Components(bkg_component), ROOT.RooFit.LineStyle(2), ROOT.RooFit.Name("bkg"), ROOT.RooFit.Range("Full"), ROOT.RooFit.NormRange("Full"))
        # background.plotOn(xframe,ROOT.RooFit.LineColor(3),ROOT.RooFit.Name("bkg"), ROOT.RooFit.LineStyle(2), ROOT.RooFit.Range("Full"), ROOT.RooFit.NormRange("Full")) 
        # model.plotOn(xframe2,ROOT.RooFit.Components(bkg_component), ROOT.RooFit.VisualizeError(result,1), ROOT.RooFit.Name("errorband"), ROOT.RooFit.FillColor(ROOT.kOrange), ROOT.RooFit.Range("Full"), ROOT.RooFit.NormRange("Full"))
	model.plotOn(xframe2,ROOT.RooFit.VisualizeError(result,1,ROOT.kFALSE), ROOT.RooFit.DrawOption("F"), ROOT.RooFit.Components(bkg_component), ROOT.RooFit.Name("errorband"), ROOT.RooFit.FillColor(ROOT.kOrange), ROOT.RooFit.Range("Full"), ROOT.RooFit.NormRange("Full"))
        # background.plotOn(xframe2,ROOT.RooFit.VisualizeError(result1,1), ROOT.RooFit.Name("errorband"), ROOT.RooFit.FillColor(ROOT.kOrange)) 
        data.plotOn(xframe2, ROOT.RooFit.Name("data"))
        model.plotOn(xframe2,ROOT.RooFit.LineColor(2),ROOT.RooFit.Components(bkg_component), ROOT.RooFit.Name("bkg"), ROOT.RooFit.LineStyle(1), ROOT.RooFit.Range("Full"), ROOT.RooFit.NormRange("Full"))

        xframe2.Draw()
        xframe2.GetYaxis().SetRangeUser(5400,6400)                                                                                                                                                       
        xframe2.GetYaxis().SetTitle("Events/ 0.01 GeV")
        xframe2.GetYaxis().SetTitleSize(0.05)
        xframe2.GetYaxis().SetLabelSize(0.045)
        xframe2.GetYaxis().SetTitleOffset(0.95)

	box = ROOT.TBox(float(xsigdown),650,float(xsigup),900) 
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

        #######################################



        # pol = ROOT.TF1("pol","[0] + [1]*x + [2]*x**2 + [3]*x**3 + [4]*x**4 + [5]*x**5 + [6]*x**6", float(xfitdown) , float(xfitup))
        # pol.SetParameters(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0)

        # pol = ROOT.TF1("pol","[0] + [1]*x + [2]*x**2", float(xfitdown) , float(xfitup))
        # pol.SetParameters(2.0, 1.0, 1.0)

 
        # r = h1[j].Fit("pol", "S", "SAME") 
        # f = h1[j].GetFunction("pol")
        # print f.GetChisquare()
        # print f.GetNDF()


        # y_pred = []
        # for i in range(h1[j].GetNbinsX()):                                                                                       
          
        #         if(h1[j].GetBinCenter(i+1) > float(xfitdown)):
        #                 y_pred.append(h1[j].GetFunction("pol").Eval(h1[j].GetBinCenter(i+1)))
                        # print(h1.GetBinCenter(i)),",",                                                                            
                        # print(h1.GetFunction("p3").Eval(h1.GetBinCenter(i))),",",                                                 

                
###########################################################################################

#         import numpy as np
#         from matplotlib import pyplot as plt

#         from sklearn.gaussian_process import GaussianProcessRegressor
#         from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C, Matern, DotProduct, WhiteKernel

#         np.random.seed(1)

#         x_masked = np.atleast_2d(x_masked).T
#         y_masked = np.atleast_2d(y_masked).T 

#         print x_masked.shape, y_masked.shape

#         #logX = np.log(X)
#         #logy = np.log(y)

#         # Mesh the input space for evaluations of the real function, the prediction and
#         # its MSE
#         x = np.atleast_2d(np.linspace(float(xfitdown),float(xfitup),int(bins))).T

#         # Instantiate a Gaussian Process model
#         # kernel = RBF(length_scale=0.001, length_scale_bounds=(0.001, 100000.0))
#         # kernel = C(0.1, (1e-10, 1e10)) * RBF(length_scale=2, length_scale_bounds=(0.001, 1e10)) #change Kernal and retry
#         kernel = C(0.1, (1, 1e10)) * RBF(length_scale=2, length_scale_bounds=(6*float(sig), 1e10))
#         # kernel = C(0.1, (20, 1e10)) * RBF(length_scale=2, length_scale_bounds=(0.2, 1e10)) #change Kernal and retry           

#         gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=40, alpha=np.add(y_masked.reshape(-1),1))

#         # print(gp.kernel) 

#         # for hyperparameter in kernel.hyperparameters: print(hyperparameter)
#         # params = kernel.get_params()
#         # for key in sorted(params): print("%s : %s" % (key, params[key]))

#         # Fit to data using Maximum Likelihood Estimation of the parameters
#         gp.fit(x_masked, y_masked)

#         print(gp.kernel_)

#         for hyperparameter in kernel.hyperparameters: print(hyperparameter)
#         params = kernel.get_params()
#         for key in sorted(params): print("%s : %s" % (key, params[key]))

#         # Make the prediction on the meshed x-axis (ask for MSE as well)
#         y_pred, sigma = gp.predict(x, return_std=True)

#         # print hyperparameter
#         # for hyperparameter in kernel.hyperparameters:
#         #     print hyperparameter

#         # # Plot the function, the prediction and the 95% confidence interval based on
#         # # the MSE
#         # plt.figure(figsize=(8, 8))
#         # # plt.plot(x, f(x), 'r:', label=r'$f(x) = x\,\sin(x)$')
#         # #plt.plot(X, y, 'r.', markersize=10, label='Observations')
#         # plt.errorbar(X, y, yerr=1*np.sqrt(y.reshape(-1)), fmt='o', color = 'red', label = "Data")   

#         # plt.plot(x, y_pred, 'b-', label='Prediction')
#         # plt.fill(np.concatenate([x, x[::-1]]),
#         #          np.concatenate([y_pred - 1.00 * sigma,
#         #                         (y_pred + 1.00 * sigma)[::-1]]),
#         #          alpha=100, fc='b', ec='None')
#         # plt.xlabel('$x$')
#         # plt.ylabel('$f(x)$')
#         # #plt.xscale("log")
#         # #plt.yscale("log")
#         # # plt.ylim(-10, 20)
#         # plt.legend(loc='lower right')
#         # plt.show()

# ##########################################################################################
#         y_pred = y_pred.reshape(-1)
#         y_pred = y_pred.flatten()
#         print "y_pred", y_pred
        
#         x_pulls = [k - float(binwidth/2.0) for k in x_unmasked]

#         pulls = []
#         for i in range(h2[j].GetNbinsX()):
#                 if(h2[j].GetBinCenter(i+1)> float(xfitdown) and h2[j].GetBinCenter(i+1)< float(xfitup)):
#                         pulls.append((y_unmasked[i] - y_pred[i])/(np.sqrt(y_unmasked[i]+1)))
#         print "pulls", pulls
#         pull_width = binwidth

#         from matplotlib.gridspec import GridSpec
#         import matplotlib.patches as mpatches
#         # import mplhep as hep
#         # from matplotlib import rc
#         # plt.rc('text', usetex=True)
#         # plt.rc('font', family='serif')
#         heights = [40, 7]
#         gs_kw = dict(height_ratios=heights)
#         print(gs_kw)
#         fig, f_axes = plt.subplots(ncols=1, nrows=2, gridspec_kw=gs_kw, sharex = True)
#         plt.subplots_adjust(wspace=0, hspace=0.05)
#         #fig.set_size_inches(18.5, 10.5)
#         f_axes[0].errorbar(x_masked.reshape(-1), y_masked.reshape(-1), yerr=1*np.sqrt(np.add(y_masked.reshape(-1),1)), fmt='o', color = 'black', label = "data")   
#         line1, = f_axes[0].plot(x,y_pred , lw=2, label = "GP fit", color='red')
#         f_axes[0].fill(np.concatenate([x, x[::-1]]),                                                                                
#                  np.concatenate([y_pred - 1.00 * sigma,                                                                           
#                                 (y_pred + 1.00 * sigma)[::-1]]),                                                                  
#                        alpha=0.2, fc='b', ec='None', color = "red", label = "Fit Error")    
#         f_axes[0].errorbar(x_sigdata, y_sigdata, yerr=1*np.sqrt(np.add(y_sigdata,1)), fmt='o', color = 'darkcyan')
#         # line1, = f_axes[0].plot(x_sig,y_fit , lw=2, linestyle ="--", label = "signal", color='green')
#         f_axes[0].set_ylabel("Events/{} GeV".format(binwidth))
#         # f_axes[0].legend(loc="upper left")
#         #f_axes[0].set_ylim(9000, 18000) 
#         # f_axes[0].set_ylim(bottom=0)
#         f_axes[0].set_xlim(xfitdown-float(binwidth/2.0), xfitup+float(binwidth/2.0))
#         f_axes[0].axvspan(xsigdown, xsigup, alpha=0.15, color='cyan', label = "masked region")
#         cyan_patch = mpatches.Patch(color='cyan',alpha=0.15, label='masked region')
#         white_patch = mpatches.Patch(color='white',alpha=0.15, label='{}-{}'.format(lxybins[j,0],lxybins[j,1]))
#         f_axes[0].legend(loc="best")
        

#         f_axes[1].bar(x_pulls, pulls, pull_width, color = 'red', label="pull")
#         f_axes[1].set_ylim([-5, 5])
#         f_axes[1].set_ylabel('Pull')
#         f_axes[1].axvspan(xsigdown, xsigup, alpha=0.15, color='cyan')
    
#         #f_axes[2].errorbar(X.reshape(-1), y.reshape(-1), yerr=1*np.sqrt(y.reshape(-1)), fmt='o', color = 'black', label = "data")  
#         #line, = f_axes[2].plot(x, y_pred, lw=2, label = "fit", color='blue')
#         #f_axes[2].set_xlabel("Dimuon Mass(GeV)")
#         #f_axes[2].set_ylabel("Events/0.01 GeV")
#         #plt.yscale("log")
#         f_axes[1].set_xlabel("Dimuon Mass(GeV)") 
#         plt.text(float(xfitdown)+float(binwidth/2.0), 40, "{} < lxy < {} cm".format(lxybins[j,0],lxybins[j,1]))
#         plt.show()
#         fig.savefig('mass{}_lxy{}_{}.png'.format(mass, lxybins[j,0], lxybins[j,1]), dpi=100)

#########################################################################################
 
        print "y_pred", y_pred

        background = ROOT.TH1F("background","Histogram of background__x", int(bins), float(xfitdown), float(xfitup))
        for i in range(h1[j].GetNbinsX()):
                background.SetBinContent(i+1,y_pred[i])
                background.SetBinError(i+1,sigmaerr[i])

        x_unmasked = np.array(x_unmasked)
        # sigma = np.zeros(h1[j].GetNbinsX())
        #sigma = np.array([0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0])
        # values = r.GetConfidenceIntervals(int(bins), 1, 1, x_unmasked, sigma, 0.683, False), ","

        print "sigma", sigmaerr

        background_alphaUp = ROOT.TH1F("background_alphaUp","Histogram of background__x", int(bins), float(xfitdown), float(xfitup))
        print float(xfitdown), float(xfitup)
        for i in range(h1[j].GetNbinsX()):
                # print x[i],y_pred[i],sigma[i]
                # background_alphaUp.SetBinContent(i+1,y_pred[i]+sigmaerr[i])
                background_alphaUp.SetBinContent(i+1,y_pred[i])
                background_alphaUp.SetBinError(i+1,sigmaerr[i])

        background_alphaDown = ROOT.TH1F("background_alphaDown","Histogram of background__x", int(bins), float(xfitdown), float(xfitup))
        for i in range(h1[j].GetNbinsX()):
                #print x[i],y_pred[i],sigma[i]
                # background_alphaDown.SetBinContent(i+1,y_pred[i]-sigmaerr[i])
                background_alphaDown.SetBinContent(i+1,y_pred[i])
                background_alphaDown.SetBinError(i+1,sigmaerr[i])

        signalpulls = []
        for i in range(h2[j].GetNbinsX()):
                if(h2[j].GetBinCenter(i+1)> float(xsigdown) and h2[j].GetBinCenter(i+1)< float(xsigup)):
                        signalpulls.append((y_unmasked[i] - y_pred[i])/(np.sqrt(y_unmasked[i])))
        print "signal pulls", signalpulls

        print "signal integral", signal.Integral()
        print "bkg integral", background.Integral()
        if background.Integral() == 0:
                bkgint = background.Integral() + 0.0001
        else:
                bkgint = background.Integral()
        # print "bkgup/bkg", background_alphaUp.Integral()+0.000001/background.Integral()+0.000001
        print "sumsq_pulls", np.sqrt(sum(i*i for i in signalpulls))
        print "sumsq_fiterr", np.sqrt(sum(i*i for i in sigmaerr))

        
        datacard = open("simple-shapes-TH1_mass{}_ctau5_Lxy{}_{}_poly.txt".format(mass, lxybins[j,0],lxybins[j,1]), "w")

        datacard.write("imax 1  number of channels\n")
        datacard.write("jmax 1  number of backgrounds\n")
        datacard.write("kmax *  number of nuisance parameters (sources of systematical uncertainties)\n")
        datacard.write("------------------------------------\n")
        datacard.write("shapes * * simple-shapes-TH1_mass{}_ctau5_Lxy{}_{}_poly.root $PROCESS $PROCESS_$SYSTEMATIC\n".format(mass, lxybins[j,0],lxybins[j,1]))
        datacard.write("------------------------------------\n")
        datacard.write("bin bin1\n")
        datacard.write("observation -1\n")
        datacard.write("------------------------------------\n")
        datacard.write("bin bin1 bin1\n")
        datacard.write("process signal background\n")
        datacard.write("process 0 1\n")
        datacard.write("rate {} {}\n".format(nS, bkgint))
        # datacard.write("rate {} {}\n".format(signal.Integral(), background_norm.getVal()))
        datacard.write("------------------------------------\n")
        datacard.write("lumi lnN 1.025 1.0\n")
        # datacard.write("bgnorm lnN 1.00 {}\n".format(background_alphaUp.Integral()/background.Integral()))
        # datacard.write("bgnorm lnN 1.00 1.2\n") 
        datacard.write("alpha shapeN2 - 1 uncertainty on background shape and normalization\n")
        datacard.close() 

        outfile = TFile("simple-shapes-TH1_mass{}_ctau5_Lxy{}_{}_poly.root".format(mass, lxybins[j,0],lxybins[j,1]), "recreate")                                                       
        outfile.cd()                                                                                                                      
        signal.Write()
        signal_sigmaUp.Write()
        signal_sigmaDown.Write()
        data_obs.Write()
        background.Write()
        background_alphaUp.Write()
        background_alphaDown.Write()

        outfile.Write()
        outfile.Close()

        
        os.system('combine -M  AsymptoticLimits simple-shapes-TH1_mass{}_ctau5_Lxy{}_{}_poly.txt > com.out'.format(mass, lxybins[j,0],lxybins[j,1]))

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


        ggphipoly.write(" {}\t{} - {}\t6\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(mass, lxybins[j,0], lxybins[j,1], int(0), int(0), np.sqrt(sum(i*i for i in signalpulls)), np.sqrt(sum(i*i for i in sigmaerr)), coml_exp, coml_1sd, coml_1su, coml_2sd, coml_2su, coml_obs))

        

        # combine -M AsymptoticLimits simple-shapes-TH1_mass0p35_ctau5_Lxy{}_{}.root.format(lxybins[i,0],lxybins[i,1])
        data_obs.SetDirectory(0)
        background.SetDirectory(0)
        background_alphaUp.SetDirectory(0)
        background_alphaDown.SetDirectory(0)


#os.system('combineCards.py simple-shapes-TH1_mass0p35_ctau5_Lxy{}_{}.txt.format(lxybins[0,0],lxybins[0,1]) simple-shapes-TH1_mass0p35_ctau5_Lxy{}_{}.txt.format(lxybins[1,0],lxybins[1,1]) simple-shapes-TH1_mass0p35_ctau5_Lxy{}_{}.txt.format(lxybins[2,0],lxybins[2,1]) simple-shapes-TH1_mass0p35_ctau5_Lxy{}_{}.txt.format(lxybins[3,0],lxybins[3,1]) simple-shapes-TH1_mass0p35_ctau5_Lxy{}_{}.txt.format(lxybins[4,0],lxybins[4,1])) 
