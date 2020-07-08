import PyFunctions
from PyFunctions import *
import ROOT
import math
from array import array
import re
import json
import types
import os
import numpy as np

ntoys = 250

# masses = [1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800]
masses = [2]

for mass in masses:
    channel = "%d" % (mass)
    cardname = "simple-shapes-TH1_mass2_ctau5_Lxy1.0_3.0_poly.txt"
    print("Channel: %s" % (channel))

    # Run without toys (?)
    print("\nWithout Toys: ")
    os.system(
        "combine %s -M GoodnessOfFit --algo=saturated -m " % (cardname) + str(mass))
    KS_Fs = ROOT.TFile("higgsCombineTest.GoodnessOfFit.mH" + str(mass) + ".root")
    KS_Ts = KS_Fs.Get("limit")
    KS_Vs = []


    for i in range(0, KS_Ts.GetEntries()):
        KS_Ts.GetEntry(i)
        if (KS_Ts.limit < 10000):
            KS_Vs.append(KS_Ts.limit)

    # Run with toys
    print("\nWith Toys: ")
    os.system(
        "combine %s -M GoodnessOfFit --algo=saturated -m " % (cardname) + str(mass) + " -t %d" % (ntoys))
    KS_F = ROOT.TFile(
        "higgsCombineTest.GoodnessOfFit.mH" + str(mass) + ".123456.root")
    KS_T = KS_F.Get("limit")
    KS_V = []
    for i in range(0, KS_T.GetEntries()):
        KS_T.GetEntry(i)
        if (KS_T.limit  <10000):
            KS_V.append(KS_T.limit)

    # Plot
    minKS = min(min(KS_V), min(KS_Vs))
    maxKS = max(max(KS_V), max(KS_Vs))
    rangeKS = maxKS - minKS
    KS_plot = ROOT.TH1F("KS_plot", "%s;Goodness Of Fit Statistic (Saturated);toys" % ("Goodness of Fit"),
                        50, minKS-(rangeKS/10), maxKS+(rangeKS/10))
    KS_plot.SetStats(0)
    for i in KS_V:
        KS_plot.Fill(i)
    GoodPlotFormat(KS_plot, "markers", ROOT.kBlack, 20)
    KS_mk = ROOT.TLine(KS_Vs[0], 0., KS_Vs[0],
                        KS_plot.GetMaximum())
    KS_mk.SetLineColor(ROOT.kRed)
    KS_mk.SetLineWidth(3)

    integral = KS_plot.Integral(1, KS_plot.FindBin(KS_Vs[0]))

    # Legend
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.AddEntry(KS_plot, "Toy Models", "pe")
    legend.AddEntry(KS_mk, "Bg, p = %.3f" % (integral / ntoys), "l")

    C_KS = ROOT.TCanvas()
    C_KS.cd()
    KS_plot.Draw("e")
    KS_mk.Draw("same")
    legend.Draw()
    ROOT.gPad.SetTicks(1, 1)
    ROOT.gPad.RedrawAxis()
    # AddCMSLumi(ROOT.gPad, plot_lumi, cmsextra)

    C_KS.Print("plots/gop_%s.png" % (channel))
    # os.system("rm *.out")
    # os.system("rm *.root")
