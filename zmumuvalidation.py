#! /usr/bin/env python

import os, sys
import argparse
import numpy as np

import ROOT as rt
rt.gROOT.SetBatch(True)

import pandas as pd

import hit

def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('--sample')
    options = ap.parse_args(sys.argv[1:])
  
    ietaRange = np.arange(1,9,1)
    drphiRange = np.arange(0.02,0.22,0.04)
    residualHistograms = dict()
    efficiencyHistograms = dict()
    for ieta in ietaRange:
        for drphi in drphiRange:
            residualHistograms[(ieta,drphi)] = rt.TH1F('residuals%d_%1.2f'%(ieta,drphi), ';#DeltaR#phi (cm);', 40, -0.25, 0.25)
            efficiencyHistograms[(ieta,drphi)] = rt.TEfficiency('efficiency_%d_%1.2f'%(ieta,drphi), ';p_{t} (GeV);', 30, 0, 100)

    rtFiles = os.listdir(options.sample)
    try:
        for f in rtFiles:
            print 'Processing file%s...'%(f)
            rtSample = rt.TFile(options.sample+'/'+f)
            sampleTree = rtSample.Get('muNtupleProducer').Get('MuDPGTree')
            sampleTree = rtSample.Get('muNtupleProducer/MuDPGTree')
            for evtId,event in enumerate(sampleTree):
                if (evtId%10000==0): print 'Processed %d/%d events...'%(evtId,sampleTree.GetEntries())
            
                # match propagated hits with rechits:
                nRechits = event.gemRecHit_nRecHits
                nProphits = event.mu_propagated_chamber.size()
                
                # chack all rechits and prophits for matches:
                for iRechit in range(nRechits):
                    for iProphit in range(nProphits):
                        recHit = hit.RecHit(event, iRechit)
                        propHit = hit.PropHit(event, iProphit)

                        for drphiMax in drphiRange: # check match for all possible drphi cuts
                            if recHit.matches(propHit, drphiMax=drphiMax):
                                # if prophit and rechit match for chosen drphi, fill histogram for hit ieta
                                residual = recHit.residual(propHit)
                                ieta = recHit.ieta
                                residualHistograms[(ieta,drphiMax)].Fill(residual)

    except KeyboardInterrupt: print 'Saving output...'

    # save residual and efficiency plots:
    outFile = rt.TFile('results/zmumu_mc/results.root', 'RECREATE')

    # fit all residual plots with two gauss and save canvases per ieta:
    residualTwoGausFit = rt.TF1('g2', 'gaus(0)+gaus(3)', -0.4, 0.4)
    residualTwoGausFit.SetParameters(30, 0, 0.05, 20, 0, 0.04)
    for ieta in ietaRange:
        residualCanvas = rt.TCanvas('ResidualCanvas', '', 800, 600)
        residualStack = rt.THStack('ResidualStack%d'%(ieta), ';#DeltaR#phi (cm);')

        residualLegend = rt.TLegend(0.12, 0.83, 0.5, 0.73)
        residualLegend.SetHeader('Cut on R#Delta#phi', 'c')
        residualLegend.SetNColumns(2)

        for i,drphi in enumerate(drphiRange):
            h = residualHistograms[(ieta,drphi)]
            h.Fit(residualTwoGausFit)
            h.SetMarkerStyle(20)
            h.SetMarkerColor(i+2)
            h.Write()
            residualStack.Add(h)
            residualLegend.AddEntry(h, '%1.2f cm'%(drphi), 'f')
        residualStack.Draw('e1nostack')
        residualLegend.Draw()
        residualCanvas.SaveAs('results/zmumu_mc/residuals/residuals_ieta%d.eps'%(ieta))

    outFile.Write()
    outFile.Close()
    
if __name__=='__main__': main()
