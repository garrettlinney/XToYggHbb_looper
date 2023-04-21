import ROOT
import argparse
from collections import OrderedDict
import copy
from datetime import date    
import glob
import numpy
import os
import plotUtils


sampleFillColor=dict()
sampleFillColor["Data"]              = None
sampleFillColor["GJets"]             = ROOT.kCyan-7
sampleFillColor["Diphoton"]          = ROOT.kSpring+5
sampleFillColor["tt+X"]              = ROOT.kPink-4
sampleFillColor["DY"]                = ROOT.kOrange+3
sampleFillColor["VG"]                = ROOT.kMagenta+2
sampleFillColor["Diboson"]           = ROOT.kPink
sampleFillColor["Single H"]          = ROOT.kGreen+2
sampleFillColor["ttH_M125"]          = ROOT.kOrange+7
sampleFillColor["HHbbgg"]            = ROOT.kBlue-4
sampleFillColor["NMSSM_XToYHTo2G2B"] = None

sampleLineColor=dict()
sampleLineColor["Data"]              = ROOT.kBlack
sampleLineColor["GJets"]             = None
sampleLineColor["Diphoton"]          = None
sampleLineColor["tt+X"]              = None
sampleLineColor["DY"]                = None
sampleLineColor["VG"]                = None
sampleLineColor["Diboson"]           = None
sampleLineColor["Single H"]          = None
sampleLineColor["ttH_M125"]          = None
sampleLineColor["HHbbgg"]            = None
sampleLineColor["NMSSM_XToYHTo2G2B"] = ROOT.kYellow

sampleLineWidth=dict()
sampleLineWidth["Data"]              = 1
sampleLineWidth["GJets"]             = 0
sampleLineWidth["Diphoton"]          = 0
sampleLineWidth["tt+X"]              = 0
sampleLineWidth["DY"]                = 0
sampleLineWidth["VG"]                = 0
sampleLineWidth["Diboson"]           = 0
sampleLineWidth["Single H"]          = 0
sampleLineWidth["ttH_M125"]          = 0
sampleLineWidth["HHbbgg"]            = 0
sampleLineWidth["NMSSM_XToYHTo2G2B"] = 2

sampleMarkerStyle=dict()
sampleMarkerStyle["Data"]              = 20
sampleMarkerStyle["GJets"]             = None
sampleMarkerStyle["Diphoton"]          = None
sampleMarkerStyle["tt+X"]              = None
sampleMarkerStyle["DY"]                = None
sampleMarkerStyle["VG"]                = None
sampleMarkerStyle["Diboson"]           = None
sampleMarkerStyle["Single H"]          = None
sampleMarkerStyle["ttH_M125"]          = None
sampleMarkerStyle["HHbbgg"]            = None
sampleMarkerStyle["NMSSM_XToYHTo2G2B"] = None

sampleMarkerSize=dict()
sampleMarkerSize["Data"]              = 1.2
sampleMarkerSize["GJets"]             = None
sampleMarkerSize["Diphoton"]          = None
sampleMarkerSize["tt+X"]              = None
sampleMarkerSize["DY"]                = None
sampleMarkerSize["VG"]                = None
sampleMarkerSize["Diboson"]           = None
sampleMarkerSize["Single H"]          = None
sampleMarkerSize["ttH_M125"]          = None
sampleMarkerSize["HHbbgg"]            = None
sampleMarkerSize["NMSSM_XToYHTo2G2B"] = None

sampleLegend=dict()
sampleLegend["Data"]              = "Data"
sampleLegend["GJets"]             = "GJets"
sampleLegend["Diphoton"]          = "Diphoton"
sampleLegend["tt+X"]              = "tt+X"
sampleLegend["DY"]                = "DY"
sampleLegend["VG"]                = "VG"
sampleLegend["Diboson"]           = "Diboson"
sampleLegend["Single H"]          = "Single H"
sampleLegend["ttH_M125"]          = "ttH"
sampleLegend["HHbbgg"]            = "ggHH"
sampleLegend["NMSSM_XToYHTo2G2B"] = "XYH"

epsilon = 1e-6

def BTagSF(tree):
  tree.Draw("weight_beforeBTagSF>>hBefore","","goff")
  hBefore = ROOT.gDirectory.Get("hBefore")
  wBefore = hBefore.GetEntries()*hBefore.GetMean()
  tree.Draw("weight_afterBTagSF>>hAfter","","goff")
  hAfter = ROOT.gDirectory.Get("hAfter")
  wAfter = hAfter.GetEntries()*hAfter.GetMean()
  SF = wBefore / wAfter
  return SF

def get_plots(samples, year, plotname, cut, plotBins, plotXTitles):

  htempDict=OrderedDict()
  nBins = plotBins[plotname][0]; lowBin = plotBins[plotname][1]; highBin = plotBins[plotname][2]
  htemp = ROOT.TH1F(plotname,"",nBins,lowBin,highBin);

  if year!="all" and year!="2016":
    years=[year]
  elif year=="2016":
    years=["2016nonAPV","2016APV"]
  else:
    years=["2016nonAPV","2016APV","2017","2018"]
  for tyear in years:
    for i,sample in enumerate(samples):
      if "NMSSM_XToYHTo2G2B" in sample:
        for mass in args.signalMass:
          infile = ROOT.TFile(args.inDir+"output_"+sample+"_"+mass+"_"+tyear+".root")
          tree = infile.Get("tout")
          if tree.GetEntries() == 0:
            print("0 entries for sample %s%s, skipping..."%(sample,mass))
            continue
          tree.Draw(plotname+">>htemp("+str(nBins)+","+str(lowBin)+","+str(highBin)+")",cut,"goff")
          if (sample+"_"+mass) not in htempDict.keys():
            htempDict[sample+"_"+mass]=[]
          htemp = ROOT.gDirectory.Get("htemp")
          htemp.GetYaxis().SetTitle("Events")
          htemp.GetXaxis().SetTitle(plotXTitles[plotname])
          htempDict[sample+"_"+mass].append(copy.deepcopy(htemp))
          htempDict[sample+"_"+mass][-1].Scale(BTagSF(tree))
      elif sample=="GJets":
        for m1,m2 in zip(["40","100","200","400","600"],["100","200","400","600","Inf"]): 
          infile = ROOT.TFile(args.inDir+"output_GJets_HT-"+m1+"To"+m2+"_"+tyear+".root")
          tree = infile.Get("tout")
          tree.Draw(plotname+">>htemp("+str(nBins)+","+str(lowBin)+","+str(highBin)+")",cut,"goff")
          if sample not in htempDict.keys():
            htempDict[sample]=[]
          htemp = ROOT.gDirectory.Get("htemp")
          htemp.GetYaxis().SetTitle("Events")
          htemp.GetXaxis().SetTitle(plotXTitles[plotname])
          htempDict[sample].append(copy.deepcopy(htemp))
          htempDict[sample][-1].Scale(BTagSF(tree))
      else:
        infile = ROOT.TFile(args.inDir+"output_"+sample+"_"+tyear+".root")
        tree = infile.Get("tout")
        tree.Draw(plotname+">>htemp("+str(nBins)+","+str(lowBin)+","+str(highBin)+")",cut,"goff")
        if sample not in htempDict.keys():
          htempDict[sample]=[]
        htemp = ROOT.gDirectory.Get("htemp")
        htemp.GetYaxis().SetTitle("Events")
        htemp.GetXaxis().SetTitle(plotXTitles[plotname])
        htempDict[sample].append(copy.deepcopy(htemp))
        htempDict[sample][-1].Scale(BTagSF(tree))

  plotDict=OrderedDict()
  groupedSamples = OrderedDict()
  tempGroups = OrderedDict()
  tempGroups["Diphoton"] = ["DiPhotonLow","DiPhoton"]
  tempGroups["tt+X"] = ["TTGG","TTGJets","TTJets"]
  tempGroups["VG"] = ["WG","ZG"]
  tempGroups["Diboson"]   = ["WW","WZ","ZZ"]
  tempGroups["Single H"]   = ["ggHToDiPhoM125","VBFH_M125","VH_M125"]
  for sample in htempDict.keys():
    if sample in tempGroups["Diphoton"]:
      if "Diphoton" not in groupedSamples.keys():
        groupedSamples["Diphoton"]=[]
      groupedSamples["Diphoton"].append(sample)
    elif sample in tempGroups["tt+X"]:
      if "tt+X" not in groupedSamples.keys():
        groupedSamples["tt+X"]=[]
      groupedSamples["tt+X"].append(sample)
    elif sample in tempGroups["VG"]:
      if "VG" not in groupedSamples.keys():
        groupedSamples["VG"]=[]
      groupedSamples["VG"].append(sample)
    elif sample in tempGroups["Diboson"]:
      if "Diboson" not in groupedSamples.keys():
        groupedSamples["Diboson"]=[]
      groupedSamples["Diboson"].append(sample)
    elif sample in tempGroups["Single H"]:
      if "Single H" not in groupedSamples.keys():
        groupedSamples["Single H"]=[]
      groupedSamples["Single H"].append(sample)
    else:
      groupedSamples[sample] = [sample]

  for gsample in groupedSamples.keys():
    tplot=None
    for sample in groupedSamples[gsample]:
      for tsample in htempDict.keys():
        if not tsample==sample:
          continue
        for htemp in htempDict[tsample]:
          if not tplot:
            tplot = copy.deepcopy(htemp)
          else:
            tplot.Add(htemp)

    #for b in range(0, tplot.GetNbinsX()+2):
    #  if tplot.GetBinContent(b)<0.0 or numpy.isnan(tplot.GetBinContent(b)) or not numpy.isfinite(tplot.GetBinContent(b)):
    #    tplot.SetBinContent(b,0.0)
    #    tplot.SetBinError(b,0.0)

    plotDict[gsample] = tplot

  return plotDict

def customize_plot(sample, plot, fillColor, lineColor, lineWidth, markerStyle, markerSize):

  error = ROOT.TMath.Sqrt(plot.GetBinError(0)*plot.GetBinError(0)+plot.GetBinError(1)*plot.GetBinError(1))
  plot.SetBinContent(1, plot.GetBinContent(1) + plot.GetBinContent(0))
  plot.SetBinError(1, error)
  plot.SetBinContent(0, 0.0)
  plot.SetBinError(0, 0.0)

  error = ROOT.TMath.Sqrt(plot.GetBinError(plot.GetNbinsX()+1)*plot.GetBinError(plot.GetNbinsX()+1)+plot.GetBinError(plot.GetNbinsX())*plot.GetBinError(plot.GetNbinsX()))
  plot.SetBinContent(plot.GetNbinsX(), plot.GetBinContent(plot.GetNbinsX()+1) + plot.GetBinContent(plot.GetNbinsX()))
  plot.SetBinError(plot.GetNbinsX(), error)
  plot.SetBinContent(plot.GetNbinsX()+1, 0.0)
  plot.SetBinError(plot.GetNbinsX()+1, 0.0)

  if fillColor: 
    plot.SetFillColor(fillColor)
    plot.SetLineColor(fillColor)
    plot.SetMarkerColor(fillColor)
  if lineColor: 
    plot.SetLineColor(lineColor)
    plot.SetMarkerColor(lineColor)
  if lineWidth:
    plot.SetLineWidth(lineWidth)
  if markerStyle:
    plot.SetMarkerStyle(markerStyle)
  if markerSize:
    plot.SetMarkerSize(markerSize)
  #plot.Sumw2()

  ### Remove spikes
  if sample!="Data" and not "met_pt" in plot.GetName():
    for b in range(1, plot.GetNbinsX()+1):
      if plot.GetBinContent(b)>0 and plot.GetBinError(b)/plot.GetBinContent(b)>0.75:
        plot.SetBinContent(b,0.0)
        plot.SetBinError(b,0.0)

  return plot

def get_yields(plotDict, plotname, lumi, year, plotData=False):
  curYields=OrderedDict()
  curErrors=OrderedDict()
  totalSMYield = None
  totalSMError = None

  for i,sample in enumerate(plotDict.keys()):
    # Signal
    if "NMSSM_XToYHTo2G2B" in sample:
      curYields[sample] = plotDict[sample].GetBinContent(1)
      curErrors[sample] = plotDict[sample].GetBinError(1)
      print(sample.replace(" ","_")+":\t%.2f +/- %.2f"%(curYields[sample],curErrors[sample]))
    # Data
    elif sample=="Data": 
      if plotData:
        curYields[sample] = plotDict[sample].GetBinContent(1)
        curErrors[sample] = plotDict[sample].GetBinError(1)
        print(sample.replace(" ","_")+":\t%.2f +/- %.2f"%(curYields[sample],curErrors[sample]))
    # Bkg
    else:
      curYields[sample] = plotDict[sample].GetBinContent(1)
      curErrors[sample] = plotDict[sample].GetBinError(1)
      print(sample.replace(" ","_")+":\t%.2f +/- %.2f"%(curYields[sample],curErrors[sample]))
      if not totalSMYield:
        totalSMYield = curYields[sample]
        totalSMError = curErrors[sample]*curErrors[sample]
      else:
        totalSMYield = totalSMYield + curYields[sample]
        totalSMError = totalSMError + curErrors[sample]
  totalSMError = ROOT.TMath.Sqrt(totalSMError)
  print("Total_SM:\t%.2f +/- %.2f"%(totalSMYield,totalSMError))

def draw_plot(plotDict, plotname, lumi, year, logY=True, logX=False, plotData=False, doRatio=True):

  # Labels
  latex = ROOT.TLatex()
  latex.SetTextFont(42)
  latex.SetTextAlign(31)
  latex.SetTextSize(0.04)
  latex.SetNDC(True)

  latexCMS = ROOT.TLatex()
  latexCMS.SetTextFont(61)
  latexCMS.SetTextSize(0.04)
  latexCMS.SetNDC(True)

  latexCMSExtra = ROOT.TLatex()
  latexCMSExtra.SetTextFont(52)
  latexCMSExtra.SetTextSize(0.04)
  latexCMSExtra.SetNDC(True)

  legoffset = 0.0
  latexSel = ROOT. TLatex()
  latexSel.SetTextAlign(11)
  latexSel.SetTextFont(42)
  latexSel.SetTextSize(0.03)
  latexSel.SetNDC(True)

  yearenergy=""
  if year!="all" or lumi<100.0:
    if year!="all":
      yearenergy="%.1f fb^{-1} (%s, 13 TeV)"%(lumi,year)
    else:
      yearenergy="%.1f fb^{-1} (2016-2018, 13 TeV)"%(lumi)
  else:
    yearenergy="%.0f fb^{-1} (13 TeV)"%(lumi)
  if plotData:
    cmsExtra="Preliminary"
  else:
    cmsExtra="Simulation"

  curPlots=OrderedDict()

  totalSM = None
  lowToHighBinsCumulative = True
  for i,sample in enumerate(plotDict.keys()):
    # Signal
    if "NMSSM_XToYHTo2G2B" in sample:
      model = sample.split("_M")[0]
      mass = sample.split("B_")[1]
      curPlots[sample] = copy.deepcopy(customize_plot(sample,plotDict[sample],sampleFillColor[model],sampleLineColor[model]+i%len(args.signalMass),sampleLineWidth[model],sampleMarkerStyle[model],sampleMarkerSize[model]))
      if args.shape and curPlots[sample].Integral(0,-1)>0.0:
        curPlots[sample].Scale(1.0/curPlots[sample].Integral(0,-1))
        if args.cumulative:
          curPlots[sample] = plotUtils.GetCumulative(curPlots[sample],lowToHighBinsCumulative)
    # Data
    elif sample=="Data": 
      if plotData:
        curPlots[sample] = copy.deepcopy(customize_plot(sample,plotDict[sample],sampleFillColor[sample],sampleLineColor[sample],sampleLineWidth[sample],sampleMarkerStyle[sample],sampleMarkerSize[sample]))
        if args.shape and curPlots[sample].Integral(0,-1)>0.0:
          curPlots[sample].Scale(1.0/curPlots[sample].Integral(0,-1))
        if args.cumulative:
          curPlots[sample] = plotUtils.GetCumulative(curPlots[sample],lowToHighBinsCumulative)
    # Bkg
    else:
      curPlots[sample] = copy.deepcopy(customize_plot(sample,plotDict[sample],sampleFillColor[sample],sampleLineColor[sample],sampleLineWidth[sample],sampleMarkerStyle[sample],sampleMarkerSize[sample]))
      if not totalSM:
        totalSM = curPlots[sample].Clone("totalSM")
      else:
        totalSM.Add(curPlots[sample])

  if args.dataOnly:
    totalSM = curPlots["Data"].Clone("totalSM")
  totalScale = totalSM.Integral(0,-1)
  if args.cumulative:
    totalSM = plotUtils.GetCumulative(totalSM,lowToHighBinsCumulative)
  if args.shape and totalScale>0.0:
    totalSM.Scale(1.0/totalScale)


  # Build stack
  stack = ROOT.THStack("stack","")
  if not args.dataOnly:
    for i,sample in enumerate(reversed(plotDict.keys())):
      # Bkg
      if not ("NMSSM_XToYHTo2G2B" in sample or sample=="Data"):
        if args.shape and totalScale>0.0:
          curPlots[sample].Scale(1.0/totalScale)
        if args.cumulative:
          curPlots[sample] = plotUtils.GetCumulative(curPlots[sample],lowToHighBinsCumulative)
        stack.Add(curPlots[sample])


  # Plot legends, ranges
  legendXOffsetNoSelPrint = 0.18
  legendYOffsetNoSelPrint = 0.1
  if args.data:
    legend = ROOT.TLegend(0.5,0.7,0.91,0.91)
    if args.dataOnly:
      legend = ROOT.TLegend(0.5,0.8,0.89,0.89)
  else:
    legend = ROOT.TLegend(0.5,0.7,0.89,0.89)
  legend.SetLineColor(0)
  legend.SetLineWidth(0)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)
  legend.SetNColumns(2)
  
  for sample in curPlots.keys():
    # Signal
    if "NMSSM_XToYHTo2G2B" in sample:
      model = sample.split("_M")[0]
      mass = sample.split("B_")[1]
      massX = mass.split("_")[1]
      massY = mass.split("_")[3]
      legend.AddEntry(curPlots[sample],sampleLegend[model]+" ("+massX+"/"+massY+") GeV","L")
    # Data
    elif sample=="Data": 
      if plotData:
        legend.AddEntry(curPlots[sample],sampleLegend[sample],"EPL")
    # Bkg
    else:
      legend.AddEntry(curPlots[sample], sampleLegend[sample],"F")
    

  # Define canvas
  canvas = ROOT.TCanvas("canvas","canvas",800,800)

  MCplot = copy.deepcopy(totalSM)
  g_unc = ROOT.TGraphAsymmErrors()
  g_data = ROOT.TGraphAsymmErrors()
  g_data_clone = ROOT.TGraphAsymmErrors()
  g_ratio = ROOT.TGraphAsymmErrors()
  g_ratio_unc = ROOT.TGraphAsymmErrors()
  g_ratio_signal = ROOT.TMultiGraph()

  h_axis = ROOT.TH1F()
  h_axis_ratio = ROOT.TH1F()
  h_axis = ROOT.TH1F("h_axis","", MCplot.GetNbinsX(), MCplot.GetXaxis().GetBinLowEdge(1), MCplot.GetXaxis().GetBinUpEdge(MCplot.GetNbinsX()))
  h_axis_ratio = ROOT.TH1F("h_axis_ratio","", MCplot.GetNbinsX(), MCplot.GetXaxis().GetBinLowEdge(1), MCplot.GetXaxis().GetBinUpEdge(MCplot.GetNbinsX()))
  if logX and MCplot.GetXaxis().GetBinLowEdge(1) < epsilon:
    h_axis.GetXaxis().SetRangeUser(MCplot.GetXaxis().GetBinCenter(1)-0.25*MCplot.GetXaxis().GetBinWidth(1), MCplot.GetXaxis().GetBinUpEdge(MCplot.GetNbinsX()))
    h_axis_ratio.GetXaxis().SetRangeUser(MCplot.GetXaxis().GetBinCenter(1)-0.25*MCplot.GetXaxis().GetBinWidth(1), MCplot.GetXaxis().GetBinUpEdge(MCplot.GetNbinsX()))

  doRatio=False
  if plotData:
    if not args.dataOnly:
      doRatio=True

    #plotUtils.ConvertToPoissonGraph(curPlots["Data"], g_data, drawZeros=True, drawXerr=False)
    plotUtils.ConvertToPoissonGraph(curPlots["Data"], g_data, drawZeros=False, drawXerr=False)
    g_data.SetMarkerStyle(20)
    g_data.SetMarkerSize(1.2)
    g_data.SetLineWidth(1)
    # draw with zero marker size so error bars drawn all the way to x axis in the case of 0 content
    g_data_clone = g_data.Clone()
    g_data_clone.SetMarkerSize(0.0)

    #plotUtils.GetPoissonRatioGraph(MCplot, curPlots["Data"], g_ratio, drawZeros=True, drawXerr=False, useMCErr=False)
    plotUtils.GetPoissonRatioGraph(MCplot, curPlots["Data"], g_ratio, drawZeros=False, drawXerr=False, useMCErr=False)
    g_ratio.SetMarkerStyle(20)
    g_ratio.SetMarkerSize(1.2)
    g_ratio.SetLineWidth(1)

  if not plotData and args.shape and doSignalMCRatio:
    doRatio=True

    for i,sample in enumerate(plotDict.keys()):
      # Signal
      if "NMSSM_XToYHTo2G2B" in sample:
        model = sample.split("_M")[0] 
        mass = sample.split("B_")[1]
        g_signal_temp = ROOT.TGraphAsymmErrors()
        plotUtils.ConvertToPoissonGraph(curPlots[sample], g_signal_temp, drawZeros=False, drawXerr=False, drawYerr=False)
        g_signal_temp.SetMarkerStyle(20)
        g_signal_temp.SetMarkerSize(1.2)
        g_signal_temp.SetLineWidth(1)

        # draw with zero marker size so error bars drawn all the way to x axis in the case of 0 content
        g_signal_temp_clone = g_signal_temp.Clone()
        g_signal_temp_clone.SetMarkerSize(0.0)

        g_ratio_signal_temp = ROOT.TGraphAsymmErrors()
        plotUtils.GetPoissonRatioGraph(MCplot, curPlots[sample], g_ratio_signal_temp, drawZeros=False, drawXerr=False, drawYerr=False, useMCErr=False)
        g_ratio_signal_temp.SetMarkerStyle(20)
        g_ratio_signal_temp.SetMarkerSize(1.2)
        g_ratio_signal_temp.SetMarkerColor(sampleLineColor[model]+i%len(args.signalMass))
        g_ratio_signal_temp.SetLineWidth(1)
        g_ratio_signal.Add(copy.deepcopy(g_ratio_signal_temp))

  for b in range(1,MCplot.GetNbinsX()+1):
    thisPoint = g_ratio_unc.GetN()
    yerror = MCplot.GetBinError(b)
    g_unc.SetPoint(thisPoint, MCplot.GetBinCenter(b), MCplot.GetBinContent(b))
    g_unc.SetPointError(thisPoint, 0.5*MCplot.GetBinWidth(b), 0.5*MCplot.GetBinWidth(b), yerror, yerror)
    if MCplot.GetBinContent(b)>0.0:
      yerror = yerror/MCplot.GetBinContent(b)
    else:
      yerror = 0.0
    g_ratio_unc.SetPoint(thisPoint, MCplot.GetBinCenter(b), 1.0)
    g_ratio_unc.SetPointError(thisPoint, 0.5*MCplot.GetBinWidth(b), 0.5*MCplot.GetBinWidth(b), yerror, yerror)
  g_unc.SetFillStyle(3244)
  g_unc.SetFillColor(ROOT.kGray+3)
  g_ratio_unc.SetFillStyle(1001)
  g_ratio_unc.SetFillColor(ROOT.kGray)

  pads = []
  if doRatio==True:
    minR=0.0
    maxR=2.0
    ty = numpy.array([])
    tmax=maxR
    if args.data:
      ty = g_ratio.GetY()
    else:
      ty = g_ratio_signal.GetY()
    if len(ty)>0:
      tmax = numpy.amax(ty)
    if tmax>maxR:
      maxR=tmax*1.05
    if maxR>5.0:
      minR=0.1
    h_axis_ratio.GetYaxis().SetRangeUser(minR,maxR)
    h_axis_ratio.SetMinimum(minR)
    h_axis_ratio.SetMaximum(maxR)
    h_axis_ratio.SetTitle(";;Data / MC")
    h_axis_ratio.GetYaxis().SetTitleSize(0.16)
    h_axis_ratio.GetYaxis().SetTitleOffset(0.25)
    if logY:
      h_axis_ratio.GetYaxis().SetTitleOffset(0.3)
    h_axis_ratio.GetYaxis().SetLabelSize(0.12)
    h_axis_ratio.GetYaxis().CenterTitle()
    h_axis_ratio.GetYaxis().SetTickLength(0.02)
    h_axis_ratio.GetXaxis().SetLabelSize(0)
    h_axis_ratio.GetXaxis().SetTitle("")
    h_axis_ratio.GetXaxis().SetTickSize(0.06)

    line = ROOT.TLine(h_axis.GetXaxis().GetBinLowEdge(1), 1.0, h_axis.GetXaxis().GetBinUpEdge(h_axis.GetNbinsX()), 1.0)

    pads.append(ROOT.TPad("1","1",0.0,0.18,1.0,1.0))
    pads.append(ROOT.TPad("2","2",0.0,0.0,1.0,0.19))
    pads[0].SetTopMargin(0.08)
    pads[0].SetBottomMargin(0.13)
    pads[0].SetRightMargin(0.05)
    pads[0].SetLeftMargin(0.10)
    pads[1].SetRightMargin(0.05)
    pads[1].SetLeftMargin(0.10)
    pads[0].Draw()
    pads[1].Draw()
    pads[1].cd()
    if maxR>5.0:
      pads[1].SetLogy()
    pads[1].SetTickx()
    if logX:
      h_axis_ratio.GetXaxis().SetMoreLogLabels()
      pads[1].SetLogx()
    if plotData:
      h_axis_ratio.Draw("")
      g_ratio_unc.Draw("SAME,2")
      g_ratio.Draw("SAME,P0")
    else:
      g_ratio_signal.Draw("SAME,P0")
      g_ratio_signal.GetXaxis().SetLimits(h_axis.GetXaxis().GetBinLowEdge(1),h_axis.GetXaxis().GetBinUpEdge(h_axis.GetNbinsX()));
      g_ratio_signal.GetHistogram().GetXaxis().SetRangeUser(h_axis.GetXaxis().GetBinLowEdge(1),h_axis.GetXaxis().GetBinUpEdge(h_axis.GetNbinsX()));
      g_ratio_signal.GetHistogram().GetYaxis().SetRangeUser(0.,2.0);

      g_ratio_signal.GetHistogram().SetTitle(";;Signal / MC")
      g_ratio_signal.GetHistogram().GetYaxis().SetTitleSize(0.16)
      g_ratio_signal.GetHistogram().GetYaxis().SetTitleOffset(0.25)
      g_ratio_signal.GetHistogram().GetYaxis().SetLabelSize(0.12)
      g_ratio_signal.GetHistogram().GetYaxis().CenterTitle()
      g_ratio_signal.GetHistogram().GetYaxis().SetTickLength(0.02)

      g_ratio_signal.GetHistogram().GetXaxis().SetLabelSize(0)
      g_ratio_signal.GetHistogram().GetXaxis().SetTitle("")
      g_ratio_signal.GetHistogram().GetXaxis().SetTickSize(0.06)
      if logX:
        if MCplot.GetXaxis().GetBinLowEdge(1) < epsilon:
          g_ratio_signal.GetHistogram().GetXaxis().SetRangeUser(MCplot.GetXaxis().GetBinCenter(1)-0.25*MCplot.GetXaxis().GetBinWidth(1), MCplot.GetXaxis().GetBinUpEdge(MCplot.GetNbinsX()))
        g_ratio_signal.GetHistogram().GetXaxis().SetMoreLogLabels()
        pads[1].SetLogx()

    #
    line.SetLineStyle(2)
    line.SetLineColor(sampleLineColor["Data"])
    line.SetLineWidth(1)
    line.Draw("SAME")
    #
    #pads[1].RedrawAxis()
    pads[1].Modified();
    pads[1].Update();

  else:
    pads.append(ROOT.TPad("1","1",0,0,1,1))
    pads[0].Draw()

  pads[0].cd()
  if logY:
    pads[0].SetLogy()
  if logX:
    h_axis.GetXaxis().SetMoreLogLabels()
    pads[0].SetLogx()


  #plot data, stack, signal, data  
  h_axis.GetYaxis().SetTitleSize(0.04)
  if args.dataOnly:
    h_axis.GetYaxis().SetTitleOffset(1.35)
  h_axis.GetXaxis().SetTitleSize(0.04)
  h_axis.GetXaxis().SetTitleOffset(1.25)
  h_axis.GetXaxis().SetTitle(MCplot.GetXaxis().GetTitle())
  if args.shape:
    h_axis.GetYaxis().SetTitle("A.U.")
  else:
    h_axis.GetYaxis().SetTitle(MCplot.GetYaxis().GetTitle())
  h_axis.GetYaxis().SetLabelSize(0.03)
  if not args.shape:
    h_axis.GetYaxis().SetMaxDigits(3)
  h_axis.Draw("")
  if not args.dataOnly:
    stack.Draw("HIST,SAME")
    g_unc.Draw("SAME,2")
  histMax = 0.0
  if plotData:
    if histMax < curPlots["Data"].GetMaximum():
      histMax = curPlots["Data"].GetMaximum()
    g_data.Draw("P,SAME")
    g_data_clone.Draw("P,SAME")
  for sample in curPlots.keys():
    if "NMSSM_XToYHTo2G2B" in sample:
      if histMax < curPlots[sample].GetMaximum(): 
        histMax = curPlots[sample].GetMaximum()
      curPlots[sample].Draw("HIST,SAME")

  if histMax < MCplot.GetMaximum(): 
    histMax = MCplot.GetMaximum()
  if logY:
    histMax = histMax*1e3
    h_axis.SetMinimum(1e-3)
  h_axis.SetMaximum(1.1*histMax)

  legend.Draw()
  pads[0].Update()
  pads[0].RedrawAxis()


  # Draw CMS headers
  expoffset=0.03
  if logY or 1.1*histMax<1000.0:
    expoffset=0
  if doRatio:
    latex.DrawLatex(0.95, 0.93+expoffset, yearenergy);
    latexCMS.DrawLatex(0.13,0.88+expoffset,"CMS");
    latexCMSExtra.DrawLatex(0.13,0.835+expoffset, cmsExtra);
  else:
    latex.DrawLatex(0.90, 0.91+expoffset, yearenergy);
    latexCMS.DrawLatex(0.13,0.86+expoffset,"CMS");
    latexCMSExtra.DrawLatex(0.13,0.815+expoffset, cmsExtra);


  # Print and save
  extension = "_"+year
  if plotData:
    if args.dataOnly:
      extension = extension+"_data"
    else:
      extension = extension+"_mcdata"
  else:
    extension = extension+"_mc"
  if logX:
    extension = extension+"_logX"
  if logY:
    extension = extension+"_logY"
  if args.shape:
    extension = extension+"_areaNormalized"
  if args.cumulative:
    extension = extension+"_cumulative"
  
  canvas.SaveAs(args.outDir + plotname.replace("/","Over") + extension + (".pdf" if args.pdf else ".png"))

if __name__=="__main__":
  ROOT.gStyle.SetOptStat(0)
  ROOT.gROOT.SetBatch(1)

  user = os.environ.get("USER")
  today = date.today().strftime("%b-%d-%Y")

  parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("--inDir", default="./cpp/temp_data/", help="Choose input directory. Default: './cpp/temp_data/'")
  parser.add_argument("--outDir", default="/home/users/"+user+"/public_html/XToYHToggbb/plots_"+today, help="Choose output directory. Default: '/home/users/"+user+"/public_html/XToYHToggbb/plots_"+today+"'")
  parser.add_argument("--data", default=False, action="store_true", help="Plot data")
  parser.add_argument("--dataOnly", default=False, action="store_true", help="Plot only data, no MC bkg")
  parser.add_argument("--noSignal", default=False, action="store_true", help="Do not plot signals")
  parser.add_argument("--signalMass", default=[], nargs="+", help="Signal mass points to plot. 'all' plots/prints all of the mass points. Default: 'MX_700_MY_100'")
  parser.add_argument("--shape", default=False, action="store_true", help="Shape normalization")
  parser.add_argument("--cumulative", default=False, action="store_true", help="Cumulative distributions")
  parser.add_argument("--years", default=[], nargs="+", help="List of years to be plotted. Default: all years")
  parser.add_argument("--pdf", default=False, action="store_true", help="Output format: .pdf. Default: .png")
  parser.add_argument("--yields", default=False, action="store_true", help="Print yields instead of plotting")
  args = parser.parse_args()

  args.inDir = args.inDir.rstrip("/")+"/"
  args.outDir = args.outDir.rstrip("/")+"/"

  if not os.path.exists(args.outDir):
    os.makedirs(args.outDir)
  os.system('cp '+os.environ.get("PWD")+'/utils/index.php '+args.outDir)

  if len(args.signalMass)==0: 
    args.signalMass = ["MX_700_MY_100"]
  if args.signalMass==["all"]: 
    args.signalMass = []
    fileNames = glob.glob(args.inDir+"*NMSSM*")
    for fileName in fileNames:
      fileName = fileName.split("/")[-1].split(".")[0].split("B_")[1].split("_201")[0]
      args.signalMass.append(fileName)

  if len(args.years)==0:
    args.years = ["all"]


  # Samples
  samples=[]
  if args.data:
    samples.append("Data")
  # SM MC
  if True: #if not args.dataOnly:
    samples.append("GJets")
    #samples.append("DiPhotonLow")
    samples.append("DiPhoton")
    samples.append("TTGG")
    samples.append("TTGJets")
    samples.append("TTJets")
    samples.append("DY")
    samples.append("WG")
    samples.append("ZG")
    samples.append("WW")
    samples.append("WZ")
    samples.append("ZZ")
    samples.append("VBFH_M125")
    samples.append("VH_M125")
    samples.append("ggHToDiPhoM125")
    samples.append("ttH_M125")
    samples.append("HHbbgg")
  # Signal MC
  if not args.noSignal:
    samples.append("NMSSM_XToYHTo2G2B")


  for year in args.years:
    lumi=0.0 #fb^-1
    if year == "2018":
      lumi = 54.5
    elif year == "2017":
      lumi = 41.5
    elif year == "2016APV":
      lumi = 19.5
    elif year == "2016nonAPV":
      lumi = 16.8
    elif year == "2016":
      lumi = 19.5+16.8
    elif year == "all":
      lumi = 59.83 + 41.48 + 19.5 + 16.8


    # Cuts
    weight = "weight_central"
    cut = "1" # Default value if no cut is to be applied

    # Separate BB. EB, EE photons
    #cut = "fabs(LeadPhoton_eta) < 1.442 && fabs(SubleadPhoton_eta) < 1.442"
    #cut = "( fabs(LeadPhoton_eta) < 1.442 && fabs(SubleadPhoton_eta) > 1.566 ) || ( fabs(LeadPhoton_eta) > 1.566 && fabs(SubleadPhoton_eta) < 1.442 )"
    #cut = "fabs(LeadPhoton_eta) > 1.566 && fabs(SubleadPhoton_eta) > 1.566"

    # Proper mvaPhoID-RunIIFall17-v2-wp90 cut
    #cut = "(fabs(LeadPhoton_eta)<1.442 ? LeadPhoton_mvaID>-0.02 : LeadPhoton_mvaID>-0.26) && (fabs(SubleadPhoton_eta)<1.442 ? SubleadPhoton_mvaID>-0.02 : SubleadPhoton_mvaID>-0.26)"
    # Proper mvaPhoID-RunIIFall17-v2-wp80 cut
    #cut = "(fabs(LeadPhoton_eta)<1.442 ? LeadPhoton_mvaID>0.42 : LeadPhoton_mvaID>0.14) && (fabs(SubleadPhoton_eta)<1.442 ? SubleadPhoton_mvaID>0.42 : SubleadPhoton_mvaID>0.14)"

    # Proper DeepFlavor Loose WP cut
    #cut = "dijet_lead_btagDeepFlavB>0.0490 && dijet_sublead_btagDeepFlavB>0.0490"
    # Proper DeepFlavor Medium WP cut
    #cut = "dijet_lead_btagDeepFlavB>0.2783 && dijet_sublead_btagDeepFlavB>0.2783"

    # Diphoton mass cut
    #cut = "Diphoton_mass > 95"

    # pT/ Mgg cut
    #cut = "LeadPhoton_pt / Diphoton_mass > 0.33 && SubleadPhoton_pt / Diphoton_mass > 0.25"

    cut = weight + "*(" + cut + ")"


    # List of plots
    plotNames = []; plotBins = dict(); plotXTitles = dict()
    plotNames.append("xcand_pt"); plotBins["xcand_pt"] = [50,0,500]; plotXTitles["xcand_pt"] = "p_{T}(X)"
    plotNames.append("xcand_eta"); plotBins["xcand_eta"] = [50,-3,3]; plotXTitles["xcand_eta"] = "#eta(X)"
    plotNames.append("xcand_phi"); plotBins["xcand_phi"] = [50,3.2,3.2]; plotXTitles["xcand_phi"] = "#phi(X)"
    plotNames.append("xcand_mass"); plotBins["xcand_mass"] = [50,200,1000]; plotXTitles["xcand_mass"] = "M(X)"

    plotNames.append("LeadPhoton_pt"); plotBins["LeadPhoton_pt"] = [50,0,500]; plotXTitles["LeadPhoton_pt"] = "p_{T}(#gamma_{1})"
    plotNames.append("LeadPhoton_eta"); plotBins["LeadPhoton_eta"] = [50,-3,3]; plotXTitles["LeadPhoton_eta"] = "#eta(#gamma_{1})"
    plotNames.append("LeadPhoton_phi"); plotBins["LeadPhoton_phi"] = [50,-3.2,3.2]; plotXTitles["LeadPhoton_phi"] = "#phi(#gamma_{1})"
    plotNames.append("LeadPhoton_pixelSeed"); plotBins["LeadPhoton_pixelSeed"] = [2,0,2]; plotXTitles["LeadPhoton_pixelSeed"] = "hasPixelSeed(#gamma_{1})"
    plotNames.append("LeadPhoton_r9"); plotBins["LeadPhoton_r9"] = [50,0,2]; plotXTitles["LeadPhoton_r9"] = "R_{9}(#gamma_{1})"
    plotNames.append("LeadPhoton_sieie"); plotBins["LeadPhoton_sieie"] = [50,0,0.05]; plotXTitles["LeadPhoton_sieie"] = "#sigma_{ieie}(#gamma_{1})"
    plotNames.append("LeadPhoton_pfPhoIso03"); plotBins["LeadPhoton_pfPhoIso03"] = [50,0,20]; plotXTitles["LeadPhoton_pfPhoIso03"] = "PF Iso_{abs}^{#gamma}(#gamma_{1})"
    plotNames.append("LeadPhoton_chargedHadronIso"); plotBins["LeadPhoton_chargedHadronIso"] = [50,0,20]; plotXTitles["LeadPhoton_chargedHadronIso"] = "PF Iso_{abs}^{ch}(#gamma_{1})"
    plotNames.append("LeadPhoton_mvaID"); plotBins["LeadPhoton_mvaID"] = [20,-1,1]; plotXTitles["LeadPhoton_mvaID"] = "MVA ID(#gamma_{1})"
    plotNames.append("LeadPhoton_pt/Diphoton_mass"); plotBins["LeadPhoton_pt/Diphoton_mass"] = [50,0,2]; plotXTitles["LeadPhoton_pt/Diphoton_mass"] = "p_{T}(#gamma_{1}) / M(#gamma#gamma)"

    plotNames.append("SubleadPhoton_pt"); plotBins["SubleadPhoton_pt"] = [50,0,500]; plotXTitles["SubleadPhoton_pt"] = "p_{T}(#gamma_{2})"
    plotNames.append("SubleadPhoton_eta"); plotBins["SubleadPhoton_eta"] = [50,-3,3]; plotXTitles["SubleadPhoton_eta"] = "#eta(#gamma_{2})"
    plotNames.append("SubleadPhoton_phi"); plotBins["SubleadPhoton_phi"] = [50,-3.2,3.2]; plotXTitles["SubleadPhoton_phi"] = "#phi(#gamma_{2})"
    plotNames.append("SubleadPhoton_pixelSeed"); plotBins["SubleadPhoton_pixelSeed"] = [2,0,2]; plotXTitles["SubleadPhoton_pixelSeed"] = "hasPixelSeed(#gamma_{2})"
    plotNames.append("SubleadPhoton_r9"); plotBins["SubleadPhoton_r9"] = [50,0,2]; plotXTitles["SubleadPhoton_r9"] = "R_{9}(#gamma_{2})"
    plotNames.append("SubleadPhoton_sieie"); plotBins["SubleadPhoton_sieie"] = [50,0,0.05]; plotXTitles["SubleadPhoton_sieie"] = "#sigma_{ieie}(#gamma_{2})"
    plotNames.append("SubleadPhoton_pfPhoIso03"); plotBins["SubleadPhoton_pfPhoIso03"] = [50,0,20]; plotXTitles["SubleadPhoton_pfPhoIso03"] = "PF Iso_{abs}^{#gamma}(#gamma_{2})"
    plotNames.append("SubleadPhoton_chargedHadronIso"); plotBins["SubleadPhoton_chargedHadronIso"] = [50,0,20]; plotXTitles["SubleadPhoton_chargedHadronIso"] = "PF Iso_{abs}^{ch}(#gamma_{2})"
    plotNames.append("SubleadPhoton_mvaID"); plotBins["SubleadPhoton_mvaID"] = [20,-1,1]; plotXTitles["SubleadPhoton_mvaID"] = "MVA ID(#gamma_{2})"
    plotNames.append("SubleadPhoton_pt/Diphoton_mass"); plotBins["SubleadPhoton_pt/Diphoton_mass"] = [50,0,2]; plotXTitles["SubleadPhoton_pt/Diphoton_mass"] = "p_{T}(#gamma_{2}) / M(#gamma#gamma)"

    plotNames.append("Diphoton_pt"); plotBins["Diphoton_pt"] = [50,0,500]; plotXTitles["Diphoton_pt"] = "p_{T}(#gamma#gamma)"
    plotNames.append("Diphoton_eta"); plotBins["Diphoton_eta"] = [50,-3,3]; plotXTitles["Diphoton_eta"] = "#eta(#gamma#gamma)"
    plotNames.append("Diphoton_phi"); plotBins["Diphoton_phi"] = [50,-3.2,3.2]; plotXTitles["Diphoton_phi"] = "#phi(#gamma#gamma)"
    plotNames.append("Diphoton_mass"); plotBins["Diphoton_mass"] = [50,60,1000]; plotXTitles["Diphoton_mass"] = "M(#gamma#gamma)"
    plotNames.append("Diphoton_pt_mgg"); plotBins["Diphoton_pt_mgg"] = [20,0,3]; plotXTitles["Diphoton_pt_mgg"] = "p_{T}(#gamma#gamma)/M(#gamma#gamma)"
    plotNames.append("Diphoton_dR"); plotBins["Diphoton_dR"] = [50,0,6]; plotXTitles["Diphoton_dR"] = "#DeltaR(#gamma#gamma)"

    plotNames.append("n_jets"); plotBins["n_jets"] = [5,0,5]; plotXTitles["n_jets"] = "N_{jets}"

    plotNames.append("dijet_lead_pt"); plotBins["dijet_lead_pt"] = [50,0,500]; plotXTitles["dijet_lead_pt"] = "p_{T}(j_{1})"
    plotNames.append("dijet_lead_eta"); plotBins["dijet_lead_eta"] = [50,-3,3]; plotXTitles["dijet_lead_eta"] = "#eta(j_{1})"
    plotNames.append("dijet_lead_phi"); plotBins["dijet_lead_phi"] = [50,-3.2,3.2]; plotXTitles["dijet_lead_phi"] = "#phi(j_{1})"
    plotNames.append("dijet_lead_mass"); plotBins["dijet_lead_mass"] = [50,0,100]; plotXTitles["dijet_lead_mass"] = "M(j_{1})"
    plotNames.append("dijet_lead_btagDeepFlavB"); plotBins["dijet_lead_btagDeepFlavB"] = [50,0,1]; plotXTitles["dijet_lead_btagDeepFlavB"] = "btagDeepFlavB(j_{1})"

    plotNames.append("dijet_sublead_pt"); plotBins["dijet_sublead_pt"] = [50,0,500]; plotXTitles["dijet_sublead_pt"] = "p_{T}(j_{2})"
    plotNames.append("dijet_sublead_eta"); plotBins["dijet_sublead_eta"] = [50,-3,3]; plotXTitles["dijet_sublead_eta"] = "#eta(j_{2})"
    plotNames.append("dijet_sublead_phi"); plotBins["dijet_sublead_phi"] = [50,-3.2,3.2]; plotXTitles["dijet_sublead_phi"] = "#phi(j_{2})"
    plotNames.append("dijet_sublead_mass"); plotBins["dijet_sublead_mass"] = [50,0,100]; plotXTitles["dijet_sublead_mass"] = "M(j_{2})"
    plotNames.append("dijet_sublead_btagDeepFlavB"); plotBins["dijet_sublead_btagDeepFlavB"] = [50,0,1]; plotXTitles["dijet_sublead_btagDeepFlavB"] = "btagDeepFlavB(j_{2})"

    plotNames.append("dijet_pt"); plotBins["dijet_pt"] = [50,0,500]; plotXTitles["dijet_pt"] = "p_{T}(jj)"
    plotNames.append("dijet_eta"); plotBins["dijet_eta"] = [50,-3,3]; plotXTitles["dijet_eta"] = "#eta(jj)"
    plotNames.append("dijet_phi"); plotBins["dijet_phi"] = [50,-3.2,3.2]; plotXTitles["dijet_phi"] = "#phi(jj)"
    plotNames.append("dijet_mass"); plotBins["dijet_mass"] = [50,60,1000]; plotXTitles["dijet_mass"] = "M(jj)"
    plotNames.append("dijet_dR"); plotBins["dijet_dR"] = [50,0,6]; plotXTitles["dijet_dR"] = "#DeltaR(jj)"

    plotNames.append("pfmet_pt"); plotBins["pfmet_pt"] = [50,0,200]; plotXTitles["pfmet_pt"] = "PF MET P_{T}"
    plotNames.append("puppimet_pt"); plotBins["puppimet_pt"] = [50,0,200]; plotXTitles["puppimet_pt"] = "PUPPI MET P_{T}"
    toexclude = []

    if args.yields:
      plotname = "LeadPhoton_pixelSeed"
      plotDict = get_plots(samples, year, plotname, cut, plotBins, plotXTitles)
      get_yields(plotDict, plotname, lumi, year, args.data)
    else:
      for plotname in plotNames:
        if plotname in toexclude:
            continue
        # Open files and get trees and return plots
        plotDict = get_plots(samples, year, plotname, cut, plotBins, plotXTitles)
        draw_plot(plotDict, plotname, lumi, year, True, False, args.data)
