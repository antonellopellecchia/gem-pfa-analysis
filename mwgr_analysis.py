#! /usr/bin/env python

import os, sys
import argparse
import numpy as np
import pandas as pd

import ROOT as rt
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptStat(0)

import hit

def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('--sample')
    ap.add_argument('--output')
    options = ap.parse_args(sys.argv[1:])

    for d in ['', 'residuals', 'efficiency']:
        try: os.makedirs('%s/%s'%(options.output, d))
        except OSError: pass
  
    ietaRange = np.arange(1,9,1)
    #drphiRange = np.arange(0.03,0.3,0.03)
    #drphiRange = np.arange(1, 6, 1)
    drphiRange = np.arange(0.005, 0.03, 0.005)
    residualHistograms = dict()
    efficiencyMatchedHistograms = dict()
    efficiencyTotalHistograms = dict()
    efficiencyHistograms = dict()
    for ieta in ietaRange:
        for drphi in drphiRange:
            residualHistograms[(ieta,drphi)] = rt.TH1F('residuals_%d_%1.2f'%(ieta,drphi), ';#DeltaR#phi (cm);', 100, -7., 7.)
            efficiencyMatchedHistograms[(ieta,drphi)] = rt.TH1F('efficiencyMatched_%d_%1.2f'%(ieta,drphi), ';p_{t} (GeV);Events', 15, 0, 100)
            efficiencyTotalHistograms[(ieta,drphi)] = rt.TH1F('efficiencyTotal_%d_%1.2f'%(ieta,drphi), ';p_{t} (GeV);Events', 15, 0, 100)
            efficiencyHistograms[(ieta,drphi)] = rt.TH1F('efficiency_%d_%1.2f'%(ieta,drphi), ';p_{t} (GeV);Efficiency', 15, 0, 100)

    prophitMap = rt.TH2F('ProphitHistogram', ';Propagated x position(cm);Propagated y position (cm);Entries', 100, -300, 300, 100, -300, 300)    
    rechitMap = rt.TH2F('RechitHistogram', ';Rechit x position(cm);Rechit y position (cm);Entries', 100, -300, 300, 100, -300, 300)

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

                 # check match for all possible drphi cuts:
                for drphiMax in drphiRange:
                    # check all rechits and prophits for matches:
                    for iProphit in range(nProphits):
                        propHit = hit.PropHit(event, iProphit)
                        prophitMap.Fill(propHit.globalX, propHit.globalY)
                        # if not propHit.isGEM: continue # isGEM not implemented yet
                        if not propHit.checksFiducialCuts: continue

                        matchedRechits = list() # all rechits matching with present prophit
                        for iRechit in range(nRechits):
                            recHit = hit.RecHit(event, iRechit)
                            rechitMap.Fill(recHit.globalX, recHit.globalY)
                            if recHit.matches(propHit, drphiMax=drphiMax):
                                # fill list of rechits matching with present prophit:
                                residual = recHit.residual(propHit)
                                matchedRechits.append((recHit, residual))
                                # if prophit and rechit match for chosen drphi, fill histogram for hit ieta
                                # break # stop loop on rechits and continue to another
                        if len(matchedRechits)>0:
                            matchedRechit, matchedResidual = min(matchedRechits, key=lambda h:abs(h[1]))
                            ieta = matchedRechit.ieta
                            residualHistograms[(ieta,drphiMax)].Fill(matchedResidual)
                        # after checking matches with all rechits,
                        # add prophit to efficiency in numerator or denominator:
                        #print bool(propHit.isME11)
                        if propHit.isME11: # include in efficiency counts only prophits from ME11
                            #print (ieta,drphiMax), propHit.pt
                            if propHit.matchFound:
                                efficiencyMatchedHistograms[(ieta,drphiMax)].Fill(propHit.pt)
                            efficiencyTotalHistograms[(ieta,drphiMax)].Fill(propHit.pt)

    except KeyboardInterrupt: print 'Saving output...' # save smaller sample

    # save residual and efficiency plots:
    outFile = rt.TFile('%s/results.root'%(options.output), 'RECREATE')

    # save rechit and prophit maps
    prophitCanvas = rt.TCanvas('ProphitCanvas', '', 800, 600)
    prophitMap.Write()
    prophitMap.Draw('colz')
    prophitCanvas.SaveAs('%s/prophits.eps'%(options.output))

    rechitCanvas = rt.TCanvas('RechitCanvas', '', 800, 600)
    rechitMap.Write()
    rechitMap.Draw('colz')
    rechitCanvas.SaveAs('%s/rechits.eps'%(options.output))

    # save and draw all residual and efficiency plots per ieta:
    '''residualTwoGausFit = rt.TF1('g2', 'gaus(0)+gaus(3)', -0.4, 0.4)
    residualTwoGausFit.SetParameters(30, 0, 0.05, 20, 0, 0.04)'''
    rt.gStyle.SetOptStat(0)
    for ieta in ietaRange:
        residualStack = rt.THStack('ResidualStack%d'%(ieta), ';#DeltaR#phi (cm);')

        residualLegend = rt.TLegend(0.12, 0.9, 0.5, 0.75)
        residualLegend.SetHeader('Cut on R#Delta#phi', 'c')
        residualLegend.SetNColumns(2)

        efficiencyLegend = rt.TLegend(0.12, 0.28, 0.5, 0.13)
        efficiencyLegend.SetHeader('Cut on R#Delta#phi', 'c')
        efficiencyLegend.SetNColumns(2)

        for i,drphi in enumerate(reversed(drphiRange)):
            hResidual = residualHistograms[(ieta,drphi)]
            '''h.Fit(residualTwoGausFit)'''
            hResidual.SetMarkerStyle(20)
            hResidual.SetFillColor(i+2)
            hResidual.Write()
            residualStack.Add(hResidual)
            residualLegend.AddEntry(hResidual, '%1.3f cm'%(drphi), 'f')
            
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
        efficiencyCanvas.SaveAs('%s/efficiency/efficiency_ieta%d.eps'%(options.output, ieta))

    outFile.Write()
    outFile.Close()
    
if __name__=='__main__': main()
