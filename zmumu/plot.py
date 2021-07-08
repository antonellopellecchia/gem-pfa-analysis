#! /usr/bin/env python3

import os, sys, re
import argparse
import numpy as np

import ROOT as rt
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptFit(1)
rt.gStyle.SetOptStat(0)

import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import hit

def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('--input')
    ap.add_argument('--output')
    ap.add_argument('--borderR', nargs='+', type=float)
    ap.add_argument('--borderPhi', nargs='+', type=float)
    options = ap.parse_args(sys.argv[1:])

    for d in ['', 'efficiency/borderCut', 'residuals/borderCut']:
        try: os.makedirs('%s/%s'%(options.output, d))
        except OSError: pass

    ietaRange = np.arange(1,9,1)
    drphi = 1 # matching drphi
    #drphi = 0.15

    if options.borderR:
        rmin, rmax, rstep = options.borderR
        fiducialCutRRange = np.array(np.arange(rmin, rmax, rstep), dtype=float) # cm from eta partition borders
    else: fiducialCutRRange = np.array(np.arange(0.0, 1.1, 0.1), dtype=float) # cm from eta partition borders
    if options.borderPhi:
        phimin, phimax, phistep = options.borderPhi
        fiducialCutPhiRange = np.array(np.arange(phimin, phimax, phistep), dtype=float) # rad from chamber borders
    else: fiducialCutPhiRange = np.array(np.arange(0e-3, 10e-3, 1e-3), dtype=float) # rad from chamber borders

    #fiducialCutPhiRange = np.arange(5e-3, 90e-3, 10e-3) # mrad from chamber borders
    #fiducialCutRRange = np.arange(0.0, 1.1, 0.1) # cm from eta partition borders
  
    rootFile = rt.TFile(options.input)

    residualTwoGausFit = rt.TF1('residualTwoGausFit', 'gaus(0)+gaus(3)+expo(6)', -0.4, 0.4)
    residualTwoGausFit.SetParameters(30, 0, 0.05, 20, 0, 0.04)
    residualTwoGausFit.SetParNames('k_{1}', '#mu_{1}', '#sigma_{1}', 'k_{2}', '#mu_{2}', '#sigma_{2}', 'k_{3}', '#mu_{3}')

    spaceResolutions = list()


    for ieta in ietaRange:
        # plot efficiency vs fiducial cut in r, phi
        efficiencyVsPtCanvas = rt.TCanvas(f'EfficiencyVsPtCanvas%d{ieta}', '', 800, 600)
        efficiencyLegend = rt.TLegend(0.14, 0.30, 0.38, 0.15)
        efficiencyLegend.SetHeader('Border cut in R (cm)', 'c')
        efficiencyLegend.SetNColumns(2)

        efficiencyVsPhiCutPlot = rt.TGraphErrors()
        efficiencyVsPhiCutPlot.SetName(f'efficiencyvsphicut_{ieta}')
        efficiencyVsPhiCutPlot.SetTitle(';Border cut in R (cm);Average sector efficiency')

        efficiencyNumbers = list()

        # plot efficiencies and residuals for various fiducial cuts in phi
        iPoint = -1
        for iR,fiducialCutR in enumerate(fiducialCutRRange):
            for iPhi,fiducialCutPhi in enumerate(fiducialCutPhiRange):
                iPoint += 1
                key = 'eta%d_fiducialR%1.2e_fiducialPhi%1.2e'%(ieta, fiducialCutR, fiducialCutPhi)

                # plot residuals
                residualCanvas = rt.TCanvas(f'ResidualCanvas{key}', '', 800, 600)
                residualHistogram = rootFile.Get(f'residuals/residuals_{key}')
                print(key)
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
                latex.DrawLatexNDC(.15, .71, '#bf{Space resolution %1.1f #mum}'%(spaceResolutions[-1]))
                latex.SetTextAlign(31)
                latex.DrawLatexNDC(.9, .92, '14 TeV, PU=0')

                residualCanvas.SaveAs(f'{options.output}/residuals/borderCut/residual_{key}.eps')

                print('Space resolution %d um'%(spaceResolutions[-1]))
        
                # plot efficiency
                efficiencyVsPtCanvas = rt.TCanvas(f'EfficiencyVsPtCanvas_{key}', '', 800, 600)

                efficiencyMatchedHistogram = rootFile.Get(f'efficiency/matched/ptMatched_{key}')
                efficiencyTotalHistogram = rootFile.Get(f'efficiency/total/ptPropagated_{key}')
                efficiencyMatchedHistogram.Sumw2()
                efficiencyTotalHistogram.Sumw2()
                efficiencyHistogram = efficiencyMatchedHistogram.Clone()
                efficiencyHistogram.Divide(efficiencyMatchedHistogram, efficiencyTotalHistogram, 1, 1, 'B')
                
                # rebin to group high pt values
                nbins = efficiencyHistogram.GetNbinsX()
                oldBins = np.array([
                    efficiencyHistogram.GetXaxis().GetBinLowEdge(i+1)
                    for i in range(nbins)
                ], dtype=float)
                newBins = oldBins[oldBins<=80]
                newBins = np.append(newBins, efficiencyHistogram.GetXaxis().GetBinLowEdge(nbins+1))
                efficiencyHistogram = efficiencyHistogram.Rebin(len(newBins)-1, 'ptEfficiency', newBins)

                efficiencyHistogram.SetMarkerStyle(20)
                efficiencyHistogram.SetMarkerColor(iPhi+2)
                efficiencyHistogram.SetLineColor(iPhi+2)

                efficiencyHistogram.Draw('e1')
                efficiencyHistogram.GetYaxis().SetRangeUser(0.6, 1.1)

                # save efficiency canvas vs pt
                latex = rt.TLatex()
                latex.SetTextSize(0.043)
                latex.DrawLatexNDC(.1, .92, 'CMS #it{Simulation}')
                latex.SetTextSize(0.035)
                latex.DrawLatexNDC(.15, .83, '#bf{Z #rightarrow #mu#mu Monte Carlo sample}')
                latex.DrawLatexNDC(.15, .77, '#bf{Chamber #eta partition %d}'%(ieta))
                latex.SetTextAlign(31)
                latex.DrawLatexNDC(.9, .92, '14 TeV, PU=0')

                efficiencyVsPtCanvas.SaveAs(f'{options.output}/efficiency/borderCut/efficiency_{key}.eps')
                
                '''if iPhi==0:
                    efficiencyHistogram.Draw('e1')
                    efficiencyHistogram.GetYaxis().SetRangeUser(0.6, 1.1)
                else: efficiencyHistogram.Draw('e1same')
                efficiencyLegend.AddEntry(efficiencyHistogram, '%d mm'%(fiducialCutR*10), 'l')'''

                # calculate average efficiency over pt by rebinning to 1-bin histo
                efficiencyMatched = efficiencyMatchedHistogram.Rebin(efficiencyMatchedHistogram.GetNbinsX(), 'efficiencyMatched_{key}')
                efficiencyTotal = efficiencyTotalHistogram.Rebin(efficiencyTotalHistogram.GetNbinsX(), 'efficiencyTotal_{key}')
                efficiency = efficiencyMatched.Clone()
                efficiency.Divide(efficiencyMatched, efficiencyTotal)
                efficiencyNumber = efficiency.Integral()
                efficiencyNumbers.append(efficiencyNumber)
                #efficiencyVsBorderCutPlot.SetPoint(iPoint, fiducialCutR, fiducialCutPhi*1e3, efficiencyNumber)

        #efficiencyLegend.Draw()        
        # save efficiency canvas vs phi fiducial cut
        efficiencyNumbers = np.array(efficiencyNumbers)
        print(fiducialCutRRange, len(fiducialCutRRange))
        print(fiducialCutPhiRange, len(fiducialCutPhiRange))
        print(efficiencyNumbers, len(efficiencyNumbers))

        #efficiencyVsBorderCutPlot = rt.TGraph2D(len(efficiencyNumbers), fiducialCutRRange, fiducialCutPhiRange, efficiencyNumbers)
        efficiencyVsBorderCutPlot = rt.TGraph2D()
        efficiencyVsBorderCutPlot.SetName(f'efficiencyvsBorderCut_{ieta}')
        efficiencyVsBorderCutPlot.SetTitle(';Border cut in r (cm);Border cut in #phi (mrad);Average sector efficiency')

        iPoint = -1
        for iR,fiducialCutR in enumerate(fiducialCutRRange):
            for iPhi,fiducialCutPhi in enumerate(fiducialCutPhiRange):
                iPoint += 1
                efficiencyVsBorderCutPlot.SetPoint(iPoint, fiducialCutR, fiducialCutPhi, efficiencyNumbers[iPoint])

        efficiencyVsBorderCutPlot.SaveAs(f'{options.output}/efficiency/efficiency_{ieta}.root')

        efficiencyVsBorderCutCanvas = rt.TCanvas('EfficiencyVsFiducialCutRCanvas%d'%(ieta), '', 800, 600)
        efficiencyVsBorderCutPlot.SetMarkerStyle(0)
        efficiencyVsBorderCutPlot.Draw('colz')
        #efficiencyVsBorderCutPlot.GetYaxis().SetRangeUser(0.7, 1.1)
        latex = rt.TLatex()
        latex.SetTextSize(0.043)
        latex.DrawLatexNDC(.1, .92, 'CMS #it{Simulation}')
        latex.SetTextSize(0.035)
        latex.DrawLatexNDC(.15, .83, '#bf{Z #rightarrow #mu#mu Monte Carlo sample}')
        latex.DrawLatexNDC(.15, .77, '#bf{Chamber #eta partition %d}'%(ieta))
        latex.SetTextAlign(31)
        latex.DrawLatexNDC(.9, .92, '14 TeV, PU=0')

        efficiencyVsBorderCutCanvas.SaveAs(f'{options.output}/efficiency/efficiency_{ieta}.eps')

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
