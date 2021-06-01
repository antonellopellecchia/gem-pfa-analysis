#! /usr/bin/env python

import os, sys, re
import argparse
import numpy as np
import pandas as pd

import ROOT as rt
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptFit(1)
rt.gStyle.SetOptStat(0)

import hit

def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('--input')
    ap.add_argument('--drphi', type=float)
    ap.add_argument('--output')
    options = ap.parse_args(sys.argv[1:])

    for d in ['', 'residuals', 'efficiency']:
        try: os.makedirs('%s/%s'%(options.output, d))
        except OSError: pass
  
    rootFile = rt.TFile(options.input)
    histogramNames = rootFile.GetListOfKeys() # unused for now

    residualTwoGausFit = rt.TF1('residualTwoGausFit', 'gaus(0)+gaus(3)+expo(6)', -0.4, 0.4)
    residualTwoGausFit.SetParameters(30, 0, 0.05, 20, 0, 0.04)
    residualTwoGausFit.SetParNames('k_{1}', '#mu_{1}', '#sigma_{1}', 'k_{2}', '#mu_{2}', '#sigma_{2}', 'k_{3}', '#mu_{3}')

    ietaRange = np.arange(1,9,1)
    spaceResolutions = list()
    for ieta in ietaRange:
        residualCanvas = rt.TCanvas('ResidualCanvas%d'%(ieta), '', 800, 600)
        residualCanvas.cd()

        residualHistogram = rootFile.Get('residuals%d_%1.2f'%(ieta, options.drphi))
        residualHistogram.Draw('e1')

        residualHistogram.Fit(residualTwoGausFit)
        gaus1 = rt.TF1('g1', 'gaus(0)', -0.4, 0.4)
        gaus1.SetParameters(residualTwoGausFit.GetParameter(0),residualTwoGausFit.GetParameter(1),residualTwoGausFit.GetParameter(2))
        gaus1.SetLineColor(rt.kBlue)
        gaus1.SetLineStyle(2)
        gaus2 = rt.TF1('g2', 'gaus(0)', -0.4, 0.4)
        gaus2.SetParameters(residualTwoGausFit.GetParameter(3),residualTwoGausFit.GetParameter(4),residualTwoGausFit.GetParameter(5))
        gaus2.SetLineColor(rt.kGreen)
        gaus2.SetLineStyle(2)
        expo = rt.TF1('e', 'expo(0)', -0.4, 0.4)
        expo.SetParameters(residualTwoGausFit.GetParameter(6),residualTwoGausFit.GetParameter(7))
        expo.SetLineColor(rt.kBlack)
        expo.SetLineStyle(2)
        gaus1.Draw('same')
        gaus2.Draw('same')
        expo.Draw('same')

        a1, a2, sigma1, sigma2 = gaus1.GetParameter(0), gaus2.GetParameter(0), gaus1.GetParameter(2), gaus2.GetParameter(2)
        sigma = np.sqrt( ((a1*sigma1)**2+(a2*sigma2)**2)/(a1**2+a2**2) )
        spaceResolutions.append(sigma*1e4)
        
        residualCanvas.Update()
        residualCanvas.Draw()
        fitBox = residualHistogram.FindObject('stats')
        fitBox.SetTextSize(0.025)
        fitBox.SetTextColor(rt.kBlack)
        fitBox.SetLineColor(rt.kBlack)
        fitBox.SetFillColor(0)
        fitBox.SetFillStyle(1001)
        fitBox.SetBorderSize(1)
        x1,y1,x2,y2 = 0.6,0.45,0.9,0.9
        fitBox.SetX1NDC(x1)
        fitBox.SetY1NDC(y1)
        fitBox.SetX2NDC(x2)
        fitBox.SetY2NDC(y2)

        latex = rt.TLatex()
        latex.SetTextSize(0.043)
        latex.DrawLatexNDC(.1, .92, 'CMS #it{Simulation}')
        latex.SetTextSize(0.035)
        latex.DrawLatexNDC(.15, .83, '#bf{Z #rightarrow #mu#mu Monte Carlo sample}')
        latex.DrawLatexNDC(.15, .77, '#bf{Chamber #eta partition %d}'%(ieta))
        latex.SetTextAlign(31)
        latex.DrawLatexNDC(.9, .92, '14 TeV, PU=0')

        residualCanvas.SaveAs('%s/residuals/residuals_ieta%d.eps'%(options.output, ieta))

        print 'Space resolution %d um'%(spaceResolutions[-1])


        efficiencyCanvas = rt.TCanvas('EfficiencyCanvas%d'%(ieta), '', 800, 600)
        efficiencyCanvas.cd()
        efficiencyLegend = rt.TLegend(0.12, 0.32, 0.5, 0.13)
        efficiencyLegend.SetHeader('Cut on R#Delta#phi', 'c')
        efficiencyLegend.SetNColumns(2)

        efficiencyHistograms = dict()
        for key in histogramNames:
            name = key.GetName()
            nameMatch = re.match('efficiency_(\d)_(\d.\d\d)', name)
            if not nameMatch: continue
            eta, drphi = nameMatch.groups()
            eta, drphi = int(eta), float(drphi)
            if eta==ieta:
                efficiencyHistograms[(eta,drphi)] = rootFile.Get('efficiency_%d_%1.2f'%(eta, drphi))
                if len(efficiencyHistograms)==1:
                    efficiencyHistograms[(eta,drphi)].Draw('e1')
                    efficiencyHistograms[(eta,drphi)].GetYaxis().SetRangeUser(0, 1.2)
                else: efficiencyHistograms[(eta,drphi)].Draw('e1same')
                efficiencyLegend.AddEntry(efficiencyHistograms[(eta,drphi)], '%1.2f cm'%(drphi), 'p')

        latex = rt.TLatex()
        latex.SetTextSize(0.043)
        latex.DrawLatexNDC(.1, .92, 'CMS #it{Simulation}')
        latex.SetTextSize(0.035)
        latex.DrawLatexNDC(.15, .83, '#bf{Z #rightarrow #mu#mu Monte Carlo sample}')
        latex.DrawLatexNDC(.15, .77, '#bf{Chamber #eta partition %d}'%(ieta))
        latex.SetTextAlign(31)
        latex.DrawLatexNDC(.9, .92, '14 TeV, PU=0')

        efficiencyLegend.Draw()
        efficiencyCanvas.SaveAs('%s/efficiency/efficiency_ieta%d.eps'%(options.output, ieta))

    spaceResolutions = np.array(spaceResolutions, dtype=float)

    spaceResolutionCanvas = rt.TCanvas('SpaceResolutionCanvas', '', 800, 600)
    spaceResolutionGraph = rt.TGraph(len(ietaRange), np.array(ietaRange, dtype=float), spaceResolutions)
    spaceResolutionGraph.SetTitle(';GE1/1 eta partition;Space resolution (#mum)')
    spaceResolutionGraph.SetMarkerStyle(0)
    spaceResolutionGraph.Draw('AP')

    latex = rt.TLatex()
    latex.SetTextSize(0.043)
    latex.DrawLatexNDC(.1, .92, 'CMS #it{Simulation}')
    latex.SetTextSize(0.035)
    latex.SetTextAlign(31)
    latex.DrawLatexNDC(.85, .83, '#bf{Z #rightarrow #mu#mu Monte Carlo sample}')
    latex.DrawLatexNDC(.85, .77, '#bf{Chamber #eta partition %d}'%(ieta))
    latex.SetTextAlign(31)
    latex.DrawLatexNDC(.9, .92, '14 TeV, PU=0')

    spaceResolutionCanvas.SaveAs('%s/SpaceResolution.eps'%(options.output))
    
if __name__=='__main__': main()
