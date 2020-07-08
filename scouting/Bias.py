import PyFunctions
from PyFunctions import *
import ROOT
import math


from array import array
import re
import json
import types
import os

ntoys = 100
#years = [2017, 2018]
#masses =[0.35, 1.5, 2, 6]
name = "analysis"

cardname = "simple-shapes-TH1_mass2_ctau5_Lxy0.0_0.1_poly.txt" 
# Signal Injection
#INJ = [0., 1., 2., 5., 10.]
INJ = [5.]
for i in INJ:
    os.system("combine %s -M GenerateOnly -t %d -m 2 --saveToys --toysFrequentist --expectSignal %f -n %s%f --bypassFrequentistFit" %
              (cardname, ntoys, i, name, i))
    os.system("combine -M FitDiagnostics -d %s -m 2 --bypassFrequentistFit --skipBOnlyFit -t %d --toysFile higgsCombine%s%f.GenerateOnly.mH1600.123456.root --rMin -5 --rMax %f --saveWorkspace -n %s%f" %
              (cardname, ntoys, name, i, max(i*5, 5), name, i))

    # Plot
    F = ROOT.TFile("fitDiagnostics%s%f.root" % (name, i))
    T = F.Get("tree_fit_sb")
    # H = ROOT.TH1F("Bias Test, injected r="+str(int(i)),
    #               "%s;(#mu_{measured} - #mu_{injected})/#sigma_{#mu};toys", 40, -20., 20.)
    # T.Draw("(r-%f)/rErr>>Bias Test, injected r=%d" %
    #        (i, int(i)))

    H = ROOT.TH1F("Bias Test, injected r="+str(int(i)),
                  "Signal Injection Test;r_{measured};toys", 100, -50., 50.)
    T.Draw("r>>Bias Test, injected r=%d" %                                                                                                                                                      
           (int(i)))

    G = ROOT.TF1("f"+name+str(i), "gaus(0)", 5., 15.)
    G.SetParLimits(0, 1, 2500)
    G.SetParLimits(1, -5, 5)
    # G.SetParLimits(1, -20, 20)
    H.Fit(G)
    ROOT.gStyle.SetOptFit(1111)
    C_B = ROOT.TCanvas()
    C_B.cd()
    
    H.SetLineWidth(2)
    H.Draw("e0")
    
#C_B.Print("../results/"+NAME+"/BIAS"+str(int(i))+"_"+cardname+".root")
    C_B.Print("plots/"+"biastest_"+cardname+"_r_"+str(int(i))+".png")
#C_B.Close()
#C_B.Print(cardname+".png")
# os.system("rm *.out")
# os.system("rm *.root")
