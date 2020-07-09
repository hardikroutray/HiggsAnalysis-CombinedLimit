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

mass = "1p5"

tree_muMC = ROOT.TChain('events')
tree_muMC.Add("/cms/routray/muon_ntuples_sl7/CMSSW_10_2_5/src/MuonAnalysis/Scouting/condor/flat_dimuon_tree_ggPhi_mass{}_ct50_new.root".format(mass))
tree_mudata = ROOT.TChain('events')
tree_mudata.Add("/cms/routray/muon_ntuples_sl7/CMSSW_10_2_5/src/MuonAnalysis/Scouting/condor/flatdimuontree_2017and2018_v1.root")

lxybins = np.array([[0,0.1], [0.1,1], [1,3], [3,7], [7,11]]) 
#lxybins = np.array([[0.0,0.1]])
#print lxybins[0,0], lxybins[0,1]

ggphipoly = open("ggphigp_lim.csv", "a")
ggphipoly.write(" mass\tlxy bin\tgp_kernel_param \tsqrt(sum signalpull**2)\tsqrt(sum fiterr**2)\tExpected 50.0%: r < \tExpected 16.0%: r < \tExpected 84.0%: r < \tExpected 2.5%: r < \tExpected 97.5%: r < \tObserved Limit\n")

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

# print bins
# print int(round(bins))
# print (int(round(bins)) % 2.0)

# if int(round(bins)) % 2.0 == 0:
#         bins = bins
# else:
#         bins = bins+1

print xsigup, xsigdown, xfitup, xfitdown, bins

signal = ROOT.TH1F("signal", "Histogram of signal__x", int(bins), float(xfitdown), float(xfitup))
tree_muMC.Draw('dimuon_mass>>signal','','') 

signal_sigmaUp = ROOT.TH1F("signal_sigmaUp", "Histogram of signal__x", int(bins), float(xfitdown), float(xfitup))
tree_muMC.Draw('dimuon_mass>>signal_sigmaUp','','')

signal_sigmaDown = ROOT.TH1F("signal_sigmaDown", "Histogram of signal__x", int(bins), float(xfitdown), float(xfitup))
tree_muMC.Draw('dimuon_mass>>signal_sigmaDown','','')

#print len(lxybins)
h2 = []
h1 = []

for j in range(len(lxybins)):

        print "Looking at lxy bin----------",lxybins[j,0], "-", lxybins[j,1], "----------------"   

        signal1 = ROOT.TH1F("signal1", "signal1", int(bins), float(xfitdown), float(xfitup))

	# tree_muMC.Draw('dimuon_mass>>signal1',"lxy > {} && lxy < {} && muon1_trkiso < 0.1 && muon2_trkiso < 0.1 && dRmuon1jet > 0.3 && dRmuon2jet> 0.3 && abs(dphidimudv) < 0.02 && log(abs(detamumu/dphimumu)) < 1.25 && distPixel > 0.05 && dimuon_mass > {} && dimuon_mass < {}".format(lxybins[j,0], lxybins[j,1], xfitdown, xfitup), '')                                                                                 
	tree_muMC.Draw('dimuon_mass>>signal1',"lxy > {} && lxy < {} && muon1_trkiso < 0.1 && muon2_trkiso < 0.1 && dRmuon1jet > 0.3 && dRmuon2jet> 0.3 && abs(dphidimudv) < 0.02 && log(abs(detamumu/dphimumu)) < 1.25 && distPixel > 0.05 && dimuon_mass > {} && dimuon_mass < {}".format(0.1, 11, xfitdown, xfitup), '')      
        
	print "signal events in this bin", signal1.Integral(), "acceptance", signal1.Integral()/100000

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



        # pol = ROOT.TF1("pol","[0] + [1]*x + [2]*x**2 + [3]*x**3 + [4]*x**4 + [5]*x**5 + [6]*x**6", float(xfitdown) , float(xfitup))
        # pol.SetParameters(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0)

        # pol = ROOT.TF1("pol","[0] + [1]*x + [2]*x**2", float(xfitdown) , float(xfitup))
        # pol.SetParameters(2.0, 1.0, 1.0)

 
        # r = h1[j].Fit("pol", "S", "SAME") 
        # f = h1[j].GetFunction("pol")
        # print f.GetChisquare()
        # print f.GetNDF()
 
        y_pred = []
        # for i in range(h1[j].GetNbinsX()):                                                                                       
          
        #         if(h1[j].GetBinCenter(i+1) > float(xfitdown)):
        #                 y_pred.append(h1[j].GetFunction("pol").Eval(h1[j].GetBinCenter(i+1)))
                        # print(h1.GetBinCenter(i)),",",                                                                            
                        # print(h1.GetFunction("p3").Eval(h1.GetBinCenter(i))),",",                                                                 
###########################################################################################

        import numpy as np
        from matplotlib import pyplot as plt

        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C, Matern, DotProduct, WhiteKernel

        np.random.seed(1)

        x_masked = np.atleast_2d(x_masked).T
        y_masked = np.atleast_2d(y_masked).T 

        print x_masked.shape, y_masked.shape

        #logX = np.log(X)
        #logy = np.log(y)

        # Mesh the input space for evaluations of the real function, the prediction and
        # its MSE
        x = np.atleast_2d(np.linspace(float(xfitdown),float(xfitup),int(bins))).T

        # Instantiate a Gaussian Process model
        # kernel = RBF(length_scale=0.001, length_scale_bounds=(0.001, 100000.0))
        # kernel = C(0.1, (1e-10, 1e10)) * RBF(length_scale=2, length_scale_bounds=(0.001, 1e10)) #change Kernal and retry
        kernel = C(0.1, (1, 1e10)) * RBF(length_scale=2, length_scale_bounds=(6*float(sig), 1e10))
        # kernel = C(0.1, (20, 1e10)) * RBF(length_scale=2, length_scale_bounds=(0.2, 1e10)) #change Kernal and retry           

        gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=40, alpha=np.add(y_masked.reshape(-1),1))

        # print(gp.kernel) 

        # for hyperparameter in kernel.hyperparameters: print(hyperparameter)
        # params = kernel.get_params()
        # for key in sorted(params): print("%s : %s" % (key, params[key]))

        # Fit to data using Maximum Likelihood Estimation of the parameters
        gp.fit(x_masked, y_masked)

        print(gp.kernel_)

        gpkernel = gp.kernel_

        for hyperparameter in kernel.hyperparameters: print(hyperparameter)
        params = kernel.get_params()
        for key in sorted(params): print("%s : %s" % (key, params[key]))

        # Make the prediction on the meshed x-axis (ask for MSE as well)
        y_pred, sigma = gp.predict(x, return_std=True)

        # print hyperparameter
        # for hyperparameter in kernel.hyperparameters:
        #     print hyperparameter

        # # Plot the function, the prediction and the 95% confidence interval based on
        # # the MSE
        # plt.figure(figsize=(8, 8))
        # # plt.plot(x, f(x), 'r:', label=r'$f(x) = x\,\sin(x)$')
        # #plt.plot(X, y, 'r.', markersize=10, label='Observations')
        # plt.errorbar(X, y, yerr=1*np.sqrt(y.reshape(-1)), fmt='o', color = 'red', label = "Data")   

        # plt.plot(x, y_pred, 'b-', label='Prediction')
        # plt.fill(np.concatenate([x, x[::-1]]),
        #          np.concatenate([y_pred - 1.00 * sigma,
        #                         (y_pred + 1.00 * sigma)[::-1]]),
        #          alpha=100, fc='b', ec='None')
        # plt.xlabel('$x$')
        # plt.ylabel('$f(x)$')
        # #plt.xscale("log")
        # #plt.yscale("log")
        # # plt.ylim(-10, 20)
        # plt.legend(loc='lower right')
        # plt.show()

##########################################################################################
        y_pred = y_pred.reshape(-1)
        y_pred = y_pred.flatten()
        print "y_pred", y_pred
        
        x_pulls = [k - float(binwidth/2.0) for k in x_unmasked]

        pulls = []
        for i in range(h2[j].GetNbinsX()):
                if(h2[j].GetBinCenter(i+1)> float(xfitdown) and h2[j].GetBinCenter(i+1)< float(xfitup)):
                        pulls.append((y_unmasked[i] - y_pred[i])/(np.sqrt(y_unmasked[i]+1)))
        print "pulls", pulls
        pull_width = binwidth

        from matplotlib.gridspec import GridSpec
        import matplotlib.patches as mpatches
        # import mplhep as hep
        # from matplotlib import rc
        # plt.rc('text', usetex=True)
        # plt.rc('font', family='serif')
        heights = [40, 7]
        gs_kw = dict(height_ratios=heights)
        print(gs_kw)
        fig, f_axes = plt.subplots(ncols=1, nrows=2, gridspec_kw=gs_kw, sharex = True)
        plt.subplots_adjust(wspace=0, hspace=0.05)
        #fig.set_size_inches(18.5, 10.5)
        f_axes[0].errorbar(x_masked.reshape(-1), y_masked.reshape(-1), yerr=1*np.sqrt(np.add(y_masked.reshape(-1),1)), fmt='o', color = 'black', label = 'Data [{} - {}]'.format(lxybins[j,0],lxybins[j,1]))
        line1, = f_axes[0].plot(x,y_pred , lw=2, label = "GP fit", color='red')
        f_axes[0].fill(np.concatenate([x, x[::-1]]),                                                                                
                 np.concatenate([y_pred - 1.00 * sigma,                                                                           
                                (y_pred + 1.00 * sigma)[::-1]]),                                                                  
                       alpha=0.2, fc='b', ec='None', color = "red", label = "Fit Error")    
        f_axes[0].errorbar(x_sigdata, y_sigdata, yerr=1*np.sqrt(np.add(y_sigdata,1)), fmt='o', color = 'darkcyan')
        # line1, = f_axes[0].plot(x_sig,y_fit , lw=2, linestyle ="--", label = "signal", color='green')
        f_axes[0].set_ylabel("Events/{} GeV".format(binwidth))
        # f_axes[0].legend(loc="upper left")
        # f_axes[0].set_ylim(5400, 6400) 
        # f_axes[0].set_ylim(bottom=0)
        f_axes[0].set_xlim(xfitdown-float(binwidth/2.0), xfitup+float(binwidth/2.0))
        f_axes[0].axvspan(xsigdown, xsigup, alpha=0.15, color='cyan', label = "masked region")
        # f_axes[0].axvspan(xfitdown, xfitdown, alpha=0.0, color='white', label = '{} cm - {} cm'.format(lxybins[j,0],lxybins[j,1]))
        f_axes[0].axvspan(xfitdown, xfitdown, alpha=0.0, color='white', label = '{}'.format(gp.kernel_))
        # cyan_patch = mpatches.Patch(color='cyan',alpha=0.15, label='masked region')
        # white_patch = mpatches.Patch(color='white',alpha=0.15, label='{}-{}'.format(lxybins[j,0],lxybins[j,1]))        
        # white_patch1 = mpatches.Patch(color='white',alpha=0.15, label='{}'.format(gp.kernel_))
        f_axes[0].legend(loc="best", fontsize = 'small')
        # legend = f_axes[0].legend(loc="best", fontsize = 'small')

        # text = legend.get_texts()[4]
        # text.set_fontsize(10)

        f_axes[1].bar(x_pulls, pulls, pull_width, color = 'red', label="pull")
        f_axes[1].set_ylim([-5, 5])
        f_axes[1].set_ylabel('Pull')
        f_axes[1].axvspan(xsigdown, xsigup, alpha=0.15, color='cyan')
    
        #f_axes[2].errorbar(X.reshape(-1), y.reshape(-1), yerr=1*np.sqrt(y.reshape(-1)), fmt='o', color = 'black', label = "data")  
        #line, = f_axes[2].plot(x, y_pred, lw=2, label = "fit", color='blue')
        #f_axes[2].set_xlabel("Dimuon Mass(GeV)")
        #f_axes[2].set_ylabel("Events/0.01 GeV")
        #plt.yscale("log")
        f_axes[1].set_xlabel("Dimuon Mass(GeV)") 
        # plt.text(float(xfitdown)+float(binwidth/2.0), 50, "{} < lxy < {} cm".format(lxybins[j,0],lxybins[j,1]))
        # plt.text(float(xfitdown)+float(binwidth/2.0), 51, "gp kernel param = {}".format(gp.kernel_))
        # plt.show()
        fig.savefig('mass{}_lxy{}_{}_gp.png'.format(mass, lxybins[j,0], lxybins[j,1]), dpi=100)

#########################################################################################
 
        print "y_pred", y_pred

        background = ROOT.TH1F("background","Histogram of background__x", int(bins), float(xfitdown), float(xfitup))
        for i in range(h1[j].GetNbinsX()):
                background.SetBinContent(i+1,y_pred[i])
                background.SetBinError(i+1,sigma[i])

        x_unmasked = np.array(x_unmasked)
        # sigma = np.zeros(h1[j].GetNbinsX())
        #sigma = np.array([0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0])
        # values = r.GetConfidenceIntervals(int(bins), 1, 1, x_unmasked, sigma, 0.683, False), ","

        print "sigma", sigma

        background_alphaUp = ROOT.TH1F("background_alphaUp","Histogram of background__x", int(bins), float(xfitdown), float(xfitup))
        print float(xfitdown), float(xfitup)
        for i in range(h1[j].GetNbinsX()):
                # print x[i],y_pred[i],sigma[i]
                background_alphaUp.SetBinContent(i+1,y_pred[i])
                background_alphaUp.SetBinError(i+1,sigma[i])

        background_alphaDown = ROOT.TH1F("background_alphaDown","Histogram of background__x", int(bins), float(xfitdown), float(xfitup))
        for i in range(h1[j].GetNbinsX()):
                #print x[i],y_pred[i],sigma[i]
                background_alphaDown.SetBinContent(i+1,y_pred[i])
                background_alphaDown.SetBinError(i+1,sigma[i])

        signalpulls = []
        for i in range(h2[j].GetNbinsX()):
                if(h2[j].GetBinCenter(i+1)> float(xsigdown) and h2[j].GetBinCenter(i+1)< float(xsigup)):
                        signalpulls.append((y_unmasked[i] - y_pred[i])/(np.sqrt(y_unmasked[i])))
        print "signal pulls", signalpulls

        print "signal integral", signal.Integral()
        print "bkg integral", background.Integral()
        if background.Integral() == 0:
                bkgint = background.Integral() + 0.1
        else:
                bkgint = background.Integral()
        # print "bkgup/bkg", background_alphaUp.Integral()+0.000001/background.Integral()+0.000001
        print "sumsq_pulls", np.sqrt(sum(i*i for i in signalpulls))
        print "sumsq_fiterr", np.sqrt(sum(i*i for i in sigma))

        
        datacard = open("simple-shapes-TH1_mass{}_ctau5_Lxy{}_{}_gp.txt".format(mass, lxybins[j,0],lxybins[j,1]), "w")

        datacard.write("imax 1  number of channels\n")
        datacard.write("jmax 1  number of backgrounds\n")
        datacard.write("kmax *  number of nuisance parameters (sources of systematical uncertainties)\n")
        datacard.write("------------------------------------\n")
        datacard.write("shapes * * simple-shapes-TH1_mass{}_ctau5_Lxy{}_{}_gp.root $PROCESS $PROCESS_$SYSTEMATIC\n".format(mass, lxybins[j,0],lxybins[j,1]))
        datacard.write("------------------------------------\n")
        datacard.write("bin bin1\n")
        datacard.write("observation -1\n")
        datacard.write("------------------------------------\n")
        datacard.write("bin bin1 bin1\n")
        datacard.write("process signal background\n")
        datacard.write("process 0 1\n")
        datacard.write("rate {} {}\n".format(ns, bkgint))
        # datacard.write("rate {} 0.0001\n".format(signal.Integral()))
        datacard.write("------------------------------------\n")
        datacard.write("lumi lnN 1.025 1.0\n")
        # datacard.write("bgnorm lnN 1.00 {}\n".format(background_alphaUp.Integral()/background.Integral()))
        # datacard.write("bgnorm lnN 1.00 1.2\n") 
        datacard.write("alpha shapeN2 - 1 uncertainty on background shape and normalization\n")
        datacard.close() 

        outfile = TFile("simple-shapes-TH1_mass{}_ctau5_Lxy{}_{}_gp.root".format(mass, lxybins[j,0],lxybins[j,1]), "recreate")                                                       
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

        
        os.system('combine -M  AsymptoticLimits --rAbsAcc=0.0001 --rRelAcc=0.001 simple-shapes-TH1_mass{}_ctau5_Lxy{}_{}_gp.txt > com.out'.format(mass, lxybins[j,0],lxybins[j,1]))

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


        ggphipoly.write(" {}\t{} - {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(mass, lxybins[j,0], lxybins[j,1], gp.kernel_,np.sqrt(sum(i*i for i in signalpulls)), np.sqrt(sum(i*i for i in sigma)), coml_exp, coml_1sd, coml_1su, coml_2sd, coml_2su, coml_obs))

        

        # combine -M AsymptoticLimits simple-shapes-TH1_mass0p35_ctau5_Lxy{}_{}.root.format(lxybins[i,0],lxybins[i,1])
        data_obs.SetDirectory(0)
        background.SetDirectory(0)
        background_alphaUp.SetDirectory(0)
        background_alphaDown.SetDirectory(0)


#os.system('combineCards.py simple-shapes-TH1_mass0p35_ctau5_Lxy{}_{}.txt.format(lxybins[0,0],lxybins[0,1]) simple-shapes-TH1_mass0p35_ctau5_Lxy{}_{}.txt.format(lxybins[1,0],lxybins[1,1]) simple-shapes-TH1_mass0p35_ctau5_Lxy{}_{}.txt.format(lxybins[2,0],lxybins[2,1]) simple-shapes-TH1_mass0p35_ctau5_Lxy{}_{}.txt.format(lxybins[3,0],lxybins[3,1]) simple-shapes-TH1_mass0p35_ctau5_Lxy{}_{}.txt.format(lxybins[4,0],lxybins[4,1])) 
