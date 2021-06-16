#! /usr/bin/env python

import os, sys
import argparse
import numpy as np
import pandas as pd

import ROOT as rt
rt.gROOT.SetBatch(True)

import hit

def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('--sample')
    ap.add_argument('--output')
    ap.add_argument('--nevents', type=int)
    options = ap.parse_args(sys.argv[1:])

    for d in ['', 'residuals', 'efficiency']:
        try: os.makedirs('%s/%s'%(options.output, d))
        except OSError: pass
  
    ietaRange = np.arange(1,9,1)
    drphi = 0.15
    fiducialCutPhiRange = np.arange(5e-3, 100e-3, 5e-3) # mrad from chamber borders
    #fiducialCutRRange = np.arange(0.1, 0.10, 0.1) # mm from eta partition borders
    fiducialCutR = 0.5
    residualHistograms = dict()
    efficiencyMatchedHistograms = dict()
    efficiencyTotalHistograms = dict()
    efficiencyHistograms = dict()
    for ieta in ietaRange:
        for fiducialCutPhi in fiducialCutPhiRange:
            key = 'eta%d_fiducialPhi%1.2e'%(ieta, fiducialCutPhi)
            residualHistograms[key] = rt.TH1F('residuals_%s'%(key), ';#DeltaR#phi (cm);', 100, -0.25, 0.25)
            efficiencyMatchedHistograms[key] = rt.TH1F('ptMatched_%s'%(key), ';p_{t} (GeV);Events', 15, 0, 100)
            efficiencyTotalHistograms[key] = rt.TH1F('ptPropagated_%s'%(key), ';p_{t} (GeV);Events', 15, 0, 100)
            efficiencyHistograms[key] = rt.TH1F('ptEfficiency_%s'%(key), ';p_{t} (GeV);Efficiency', 15, 0, 100)

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
                if nRechits==0 or nProphits==0: continue
                
                '''print 'mu', event.mu_isGEM.size(), event.mu_isME11.size(), event.mu_isStandalone.size(), event.mu_isTight.size(), event.mu_nMuons
                print 'prop', event.mu_propagated_chamber.size(), event.mu_propagated_isME11.size(), event.mu_propagatedGlb_phi.size()
                print 'rec', event.gemRecHit_nRecHits, event.gemRecHit_g_phi.size()
                print 'standalone', event.mu_isStandalone
                print 'isME11', event.mu_isME11
                print '''''

                # check all rechits and prophits for matches:
                for iProphit in range(nProphits):
                    propHit = hit.PropHit(event, iProphit)
                    # if not propHit.isGEM: continue # isGEM not implemented yet
                    ieta = propHit.ieta
                    for fiducialCutPhi in fiducialCutPhiRange:
                        key = 'eta%d_fiducialPhi%1.2e'%(ieta, fiducialCutPhi)
                        if not propHit.checksFiducialCuts(fiducialCutR, fiducialCutPhi): continue

                        for iRechit in range(nRechits):
                            recHit = hit.RecHit(event, iRechit)
                            if recHit.matches(propHit, drphiMax=drphi):
                                residual = recHit.residual(propHit)
                                residualHistograms[key].Fill(residual)
                                # fill list of rechits matching with present prophit:
                                #residual = recHit.residual(propHit)
                                #matchedRechits.append((recHit, residual))
                                
                        if propHit.isME11: # include in efficiency counts only prophits from ME11
                            if propHit.matchFound:
                                efficiencyMatchedHistograms[key].Fill(propHit.pt)
                            efficiencyTotalHistograms[key].Fill(propHit.pt)

    except KeyboardInterrupt: print 'Saving output...' # save smaller sample

    # save residual and efficiency plots:
    outFile = rt.TFile('%s/results.root'%(options.output), 'RECREATE')
    outFile.mkdir('residuals')
    outFile.mkdir('efficiency')
    outFile.mkdir('efficiency/matched')
    outFile.mkdir('efficiency/total')

    # save and draw all residual and efficiency plots per ieta:
    '''residualTwoGausFit = rt.TF1('g2', 'gaus(0)+gaus(3)', -0.4, 0.4)
    residualTwoGausFit.SetParameters(30, 0, 0.05, 20, 0, 0.04)'''
    rt.gStyle.SetOptStat(0)
    for ieta in ietaRange:
        for fiducialCutPhi in fiducialCutPhiRange:
            key = 'eta%d_fiducialPhi%1.2e'%(ieta, fiducialCutPhi)
            outFile.cd('residuals')
            residualHistograms[key].Write()
            outFile.cd('efficiency/matched')
            efficiencyMatchedHistograms[key].Write()
            outFile.cd('efficiency/total')
            efficiencyTotalHistograms[key].Write()

        '''residualStack = rt.THStack('ResidualStack%d'%(ieta), ';#DeltaR#phi (cm);')

        residualLegend = rt.TLegend(0.12, 0.9, 0.5, 0.75)
        residualLegend.SetHeader('Cut on R#Delta#phi', 'c')
        residualLegend.SetNColumns(2)

        efficiencyLegend = rt.TLegend(0.12, 0.28, 0.5, 0.13)
        efficiencyLegend.SetHeader('Cut on R#Delta#phi', 'c')
        efficiencyLegend.SetNColumns(2)

        for i,drphi in enumerate(reversed(drphiRange)):
            hResidual = residualHistograms[(ieta,drphi)]
            h.Fit(residualTwoGausFit)
            hResidual.SetMarkerStyle(20)
            hResidual.SetFillColor(i+2)
            hResidual.Write()
            residualStack.Add(hResidual)
            residualLegend.AddEntry(hResidual, '%1.2f cm'%(drphi), 'f')
            
            efficiencyMatchedHistograms[(ieta,drphi)].Write()
            efficiencyTotalHistograms[(ieta,drphi)].Write()
            #efficiencyMatchedHistograms[(ieta,drphi)].Divide(efficiencyTotalHistograms[(ieta,drphi)])
            #hEfficiency = efficiencyMatchedHistograms[(ieta,drphi)]
            efficiencyMatchedHistograms[(ieta,drphi)].Sumw2()
            efficiencyTotalHistograms[(ieta,drphi)].Sumw2()
            efficiencyHistograms[(ieta,drphi)].Divide(efficiencyMatchedHistograms[(ieta,drphi)],efficiencyTotalHistograms[(ieta,drphi)], 1, 1, 'B')
            hEfficiency = efficiencyHistograms[(ieta,drphi)]
            hEfficiency.SetMarkerStyle(20)
            hEfficiency.SetMarkerColor(i+2)
            hEfficiency.Write()
            efficiencyLegend.AddEntry(hEfficiency, '%1.2f cm'%(drphi), 'p')

        residualCanvas = rt.TCanvas('ResidualCanvas%d'%(ieta), '', 800, 600)
        residualCanvas.cd()
        residualStack.Draw('nostack')
        residualLegend.Draw()
        residualCanvas.SaveAs('%s/residuals/residuals_ieta%d.eps'%(options.output, ieta))
        
        efficiencyCanvas = rt.TCanvas('EfficiencyCanvas%d'%(ieta), '', 800, 600)
        efficiencyCanvas.cd()
        for i,drphi in enumerate(drphiRange):
            if i==0: efficiencyHistograms[(ieta,drphi)].Draw('e1')
            else: efficiencyHistograms[(ieta,drphi)].Draw('e1same')
        efficiencyLegend.Draw()
        efficiencyCanvas.SaveAs('%s/efficiency/efficiency_ieta%d.eps'%(options.output, ieta))'''

    outFile.Write()
    outFile.Close()
    
if __name__=='__main__': main()
