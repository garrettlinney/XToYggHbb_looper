from collections import OrderedDict
import ROOT
import copy
import os
from datetime import date
import plotUtils
from tdrStyle import *
setTDRStyle()

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
#logger.setLevel(logging.INFO)

formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')

file_handler = logging.FileHandler('logs/compare_shape.log')
import sys
stdout_handler = logging.StreamHandler(sys.stdout)

file_handler.setFormatter(formatter)

logger.addHandler(file_handler)
logger.addHandler(stdout_handler)

user = os.environ.get("USER")
today= date.today().strftime("%b-%d-%Y")

inputbase = "../cpp/temp_data/"

xLabels_g = {
        "leadPho_sieie": "lead photon #sigma(i#eta i#eta)", 
        "leadPho_phoIso": "lead photon Iso_{pho}",
        "leadPho_trkIso": "lead photon Iso_{track}", 
        "leadPho_chgIso": "lead photon Iso_{charged}", 
        "subleadPho_sieie": "sublead photon #sigma(i#eta i#eta)",
        "subleadPho_phoIso": "sublead photon Iso_{pho}",
        "subleadPho_trkIso": "sublead photon Iso_{track}",
        "subleadPho_chgIso": "sublead photon Iso_{charged}" 
        }

colors_g = ["#D10E0E", "#3476CB", "#34CB4D", "#CB8034"]

def getHistNames(samples, years):

    fnames = {} 
    for sample in samples:
        for year in years:
            tag = "{}_{}".format(sample, year) 

            fname = "{}output_{}_{}.root".format(inputbase, sample, year)
            logger.debug("[getHistNames]: file to be opened: {}".format(fname))
            fnames[tag] = fname 

    return fnames

def makePlot(fnames, variable, properties, outdir, tag):

    style = properties["style"]
    doRatio = properties["doRatio"]
    doLog = properties["doLog"]

    hists = {}
    xRanges = [] # determine from input hists in root file
    yRanges = [] # determine from content from hists in root file
    xLabel = xLabels_g[variable]
    yLabel = "Events" #getYLabel(xRanges)

    for key, value in fnames.items():
        # open root file and get histogram
        logger.debug("[makePlot]: open file {}".format(value))
        f_ = ROOT.TFile(value)
        logger.debug("[makePlot]: get histogram {} from root file {}".format(variable, value))
        hist = f_.Get(variable)
        hist.Scale(1/hist.Integral())
        hists[key] = copy.deepcopy(hist)
        logger.debug("[makePlot]: histogram with name {} saved with key {}, integral is {}".format(hists[key].GetName(), key, hists[key].Integral()))
        xRanges.append(getXRange(hist) )
        yRanges.append(getYRange(hist) )

    xmin, xmax = getLimits(xRanges)
    ymin, ymax = getLimits(yRanges)

    # make a dummy canvas
    c = ROOT.TCanvas('c', '', 800, 800)

    #c.SetGrid()
    dummy = ROOT.TH1D("dummy","dummy",1,xmin,xmax)
    dummy.SetMinimum(ymin)
    dummy.SetMaximum(ymax*1.3)
    dummy.SetLineColor(0)
    dummy.SetMarkerColor(0)
    dummy.SetLineWidth(0)
    dummy.SetMarkerSize(0)
    dummy.GetYaxis().SetTitle(yLabel)
    dummy.GetXaxis().SetTitle(xLabel)
    dummy.GetXaxis().SetLabelOffset(0.015)
    dummy.Draw()

    # legend
    legend = ROOT.TLegend(.5,.75,.9,.9)
    legend.SetTextColor(1)

    # put hists onto canvas
    counter_color = 0
    for key, value in hists.items():
        logger.debug("[makePlot]: draw histogram from {} with name {}, integral is {}".format(key, value.GetName(), value.Integral()))
        value.Draw("{} same".format(style))

        color = colors_g[counter_color]
        if style == "HIST":
            value.SetLineColor(ROOT.TColor.GetColor(color))
            value.SetLineWidth(2)
            legend.AddEntry(value, key, "l")

        counter_color += 1

    legend.SetShadowColor(0);
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.SetTextSize(0.03)
    legend.Draw("same")

    # save .png/.pdf
    c.SaveAs("{}{}_{}.png".format(outdir, variable, tag))
    c.SaveAs("{}{}_{}.pdf".format(outdir, variable, tag))

def getYLabel():
    pass

def getXRange(hist):
    xmin = hist.GetXaxis().GetBinLowEdge(1) 
    xmax = hist.GetXaxis().GetBinUpEdge(hist.GetNbinsX()) 
    return xmin, xmax

def getYRange(hist):
    return hist.GetMinimum(), hist.GetMaximum()

def getLimits(ranges):
    min = 9999
    max = -9999
    for range in ranges:
        if range[0] < min:
            min = range[0]
        if range[1] > max:
            max = range[1]

    return min, max

def compare_shape(config):

    samples = config["samples"]
    years = config["years"]
    fnames = getHistNames(samples, years)

    variables = config["variables"]
    properties = config["properties"]

    outdir = config["outdirbase"] + config["outdirtag"]
    for variable in variables:
        # make one variable plot at a time
        makePlot(fnames, variable, properties, outdir, "")
        # save .json to the same place as plot

'''
Internal content holder
- [ ] a list of hists names, saved in ordered dict, keys contains: year, sample indentifier
- [ ] a list of variables to be plotted
- [ ] a list of histogram properties, style (marker/hist), width, color, doRatio, log scale...
    + colors should be a global dict: color1 = colors[sample]
- [ ] global settings: plot tag
'''

'''
User operation:
In a config:
    - sample, and year combination (two lists)
    - variables to be plotted
    - properties: style, doRatio, log or not (a subdict) 
'''

test_config = {
        "samples": ["ggHToDiPhoM125", "diPhoton", "HHggtautau"],
        "years": ["2018"],
        "variables": [
            "leadPho_sieie", "leadPho_phoIso", "leadPho_trkIso", "leadPho_chgIso", 
            "subleadPho_sieie", "subleadPho_phoIso", "subleadPho_trkIso", "subleadPho_chgIso"
            ],
        "outdirbase": "/home/users/hmei/public_html/Hgg/triggerStudy/checkVariables/",
        "outdirtag": "",
        "properties": {
            "style": "HIST",
            "doRatio": False,
            "doLog": False
            }
        }

# make it configurable, either with commandline or shell script
# when use it, just compare_shape.py -c configname 
compare_shape(test_config)
