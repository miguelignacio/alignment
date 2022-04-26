import ROOT, numpy as np, pandas as pd
import sys
n = 6

indir = "~/alignment/clas12/scripts/generic/cache/data_noLC"
if len(sys.argv)>1:
    indir=sys.argv[1]
if len(sys.argv)>2:
    n = int(sys.argv[2])

ROOT.gSystem.Load("~/alignment/event/libAlignEvent.so")
#ROOT.gStyle.SetTitleSize(0.05,\"XYZT\");
#ROOT.gStyle.SetPadLeftMargin(.15);\n",
#ROOT.gStyle.SetPadBottomMargin(.13);\n",
#ROOT.gStyle.SetOptStat(0);\n",
ROOT.TCanvas()
means = []
hists = []
#n = 6
for i in range(0,n):
    inputFileName=f"{indir}/plots_pass_{i+1}/prealign.root"
    inputFile = ROOT.TFile(inputFileName)
    AlignTree =inputFile.Get("AlignTree")
    AlignTree.Draw(f"AlignEvent.fChi2/AlignEvent.fNdof>>h{i+1}(100, 0, 20)", "AlignEvent.fChi2/AlignEvent.fNdof<20", "N")
    means.append(AlignTree.GetHistogram().GetMean())
    hists.append(AlignTree.GetHistogram().Clone())
print(means)
h = ROOT.TH1D("h", "#chi^{2} progress;iteration;#chi^{2}", n, 0.5, n+.5)
for i in range(0,n):
    h.SetBinContent(i+1, means[i])
h.SetMarkerStyle(21)
h.SetMinimum(0)
h.Draw("HIST")
h2 = h.Clone()
h2.Draw("SAME P")
ROOT.gStyle.SetOptStat(0)

ROOT.gPad.Draw()
ROOT.gPad.SaveAs("chi2 progress.pdf")
