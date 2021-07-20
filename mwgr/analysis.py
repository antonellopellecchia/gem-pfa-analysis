#! /usr/bin/env python3

import os, sys
import argparse
import numpy as np
import multiprocessing
import subprocess
import time

import ROOT as rt
rt.gROOT.SetBatch(True)

import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import hit

ietaRange = np.arange(1,9,1)
regionRange = [-1, 1]
layerRange = [1, 2]

def processRootFile(processNumber, rootFileIn, rootFileOut, runParameters):
    # unpack analysis parameters:
    drphi, fiducialCutR, fiducialCutPhi = runParameters

    residualHistograms = dict()
    rechitMaps = dict()
    prophitMaps = dict()
    matchedRechitMaps = dict()
    for layer in layerRange:
        for region in regionRange:
            key = f'layer_{layer}_region_{region}'
            residualHistograms[key] = rt.TH1F(f'residuals_{key}', ';#DeltaR#phi (cm);', 100, 5, 5)
            rechitMaps[key] = rt.TH2F(f'rechits_{key}', ';Chamber;Eta partition;Events', 36, 0, 36, 8, 1, 9)
            prophitMaps[key] = rt.TH2F(f'prophits_{key}', ';Chamber;Eta partition;Events', 36, 0, 36, 8, 1, 9)
            matchedRechitMaps[key] = rt.TH2F(f'matchedRechits_{key}', ';Chamber;Eta partition;Events', 36, 0, 36, 8, 1, 9)

    rtSample = rt.TFile(rootFileIn)
    #sampleTree = rtSample.Get('muNtupleProducer').Get('MuDPGTree')
    sampleTree = rtSample.Get('muNtupleProducer/MuDPGTree')
    entries = sampleTree.GetEntries()
    try:
        for evtId,event in enumerate(sampleTree):
            if (evtId%10000==0): print(f'Process {processNumber} ran on {evtId}/{entries} events...')
        
            # match propagated hits with rechits:
            nRechits = event.gemRecHit_nRecHits
            nProphits = event.mu_propagated_chamber.size()
            
            '''print 'mu', event.mu_isGEM.size(), event.mu_isME11.size(), event.mu_isStandalone.size(), event.mu_isTight.size(), event.mu_nMuons
            print 'prop', event.mu_propagated_chamber.size(), event.mu_propagated_isME11.size(), event.mu_propagatedGlb_phi.size()
            print 'rec', event.gemRecHit_nRecHits, event.gemRecHit_g_phi.size()
            print 'standalone', event.mu_isStandalone
            print 'isME11', event.mu_isME11
            print '''''

            # fill rechit map:
            for iRechit in range(nRechits):
                recHit = hit.RecHit(event, iRechit)
                ieta = recHit.ieta
                rechitKey = f'layer_{recHit.layer}_region_{recHit.region}'
                rechitMaps[rechitKey].Fill(recHit.chamber, recHit.ieta)

            # check all rechits and prophits for matches:
            for iProphit in range(nProphits):
                propHit = hit.PropHit(event, iProphit)
                # if not propHit.isGEM: continue # isGEM not implemented yet
                ieta = propHit.ieta
                prophitKey = f'layer_{propHit.layer}_region_{propHit.region}'
                prophitMaps[prophitKey].Fill(propHit.chamber, propHit.ieta)
                if not propHit.checksFiducialCuts(fiducialCutR, fiducialCutPhi): continue

                for iRechit in range(nRechits):
                    recHit = hit.RecHit(event, iRechit)
                    if recHit.matches(propHit, drphiMax=drphi):
                        #print('matches')
                        residual = recHit.residual(propHit)
                        key = f'layer_{recHit.layer}_region_{recHit.region}'
                        residualHistograms[key].Fill(residual)
                        matchedRechitMaps[key].Fill(recHit.chamber, recHit.ieta)
                        # fill list of rechits matching with present prophit:
                        #residual = recHit.residual(propHit)
                        #matchedRechits.append((recHit, residual))
                        
                '''if propHit.isME11: # include in efficiency counts only prophits from ME11
                    if propHit.matchFound:
                        efficiencyMatchedHistograms[key].Fill(propHit.pt)
                    efficiencyTotalHistograms[key].Fill(propHit.pt)'''
    except KeyboardInterrupt:
        print(f'Stopping process {processNumber} and saving output...') # save smaller sample

    # save residual and efficiency plots:
    outFile = rt.TFile(rootFileOut, 'RECREATE')
    outFile.mkdir('residuals')
    outFile.mkdir('rechits')
    outFile.mkdir('prophits')
    outFile.mkdir('rechits/matched')

    # save and draw all residual and efficiency plots per ieta:
    '''residualTwoGausFit = rt.TF1('g2', 'gaus(0)+gaus(3)', -0.4, 0.4)
    residualTwoGausFit.SetParameters(30, 0, 0.05, 20, 0, 0.04)'''
    rt.gStyle.SetOptStat(0)

    for layer in layerRange:
        for region in regionRange:
            key = f'layer_{layer}_region_{region}'
            outFile.cd('residuals')
            residualHistograms[key].Write()
            outFile.cd('rechits')
            rechitMaps[key].Write()
            outFile.cd('prophits')
            prophitMaps[key].Write()
            outFile.cd('rechits/matched')
            matchedRechitMaps[key].Write()

    outFile.Write()
    outFile.Close()

    '''for ieta in ietaRange:
        for fiducialCutR in fiducialCutRRange:
            for fiducialCutPhi in fiducialCutPhiRange:
                key = 'eta%d_fiducialR%1.2e_fiducialPhi%1.2e'%(ieta, fiducialCutR, fiducialCutPhi)
                outFile.cd('residuals')
                residualHistograms[key].Write()
                outFile.cd('efficiency/matched')
                efficiencyMatchedHistograms[key].Sumw2()
                efficiencyMatchedHistograms[key].Write()
                outFile.cd('efficiency/total')
                efficiencyTotalHistograms[key].Sumw2()
                efficiencyTotalHistograms[key].Write()'''

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

def mergeRootFiles(inputFilePaths, outputFilePath):
    print('Merging temporary output files...')
    subprocess.run(['hadd', '-f', outputFilePath]+inputFilePaths)

    return
    # save residual and efficiency plots:
    outFile = rt.TFile(outputFilePath, 'RECREATE')
    outFile.mkdir('residuals')
    outFile.mkdir('rechits')
    outFile.mkdir('prophits')
    outFile.mkdir('rechits/matched')

    tempRootFiles = list()
    for inputPath in inputFilePaths:
        tempRootFiles.append(rt.TFile(outputFilePath))
    
    for layer in layerRange:
        for region in regionRange:
            continue
            key = f'layer_{layer}_region_{region}'
            outFile.cd('residuals')
            residualHistograms[key].Write()
            outFile.cd('rechits')
            rechitMaps[key].Write()
            outFile.cd('prophits')
            prophitMaps[key].Write()
            outFile.cd('rechits/matched')
            matchedRechitMaps[key].Write()


def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument('--sample')
    ap.add_argument('--output')
    ap.add_argument('--drphi', type=float)
    #ap.add_argument('--borderR', nargs='+', type=float)
    #ap.add_argument('--borderPhi', nargs='+', type=float)
    ap.add_argument('--nevents', type=int)
    options = ap.parse_args(sys.argv[1:])

    for d in ['', 'tmp']:
        try: os.makedirs(f'{options.output}/{d}')
        except OSError: pass
  
    # analysis workflow parameters:
    if options.drphi: drphi = options.drphi # matching drphi in cm
    else: drphi = 4 # default 4 cm
    fiducialCutR = 1 # cm
    fiducialCutPhi = 5e-3 # mrad

    runParameters = (drphi, fiducialCutR, fiducialCutPhi)

    rtFiles = os.listdir(options.sample)
    rtFiles = [ f for f in rtFiles if f[-5:]=='.root' ]
    tempOutPaths = list()
    processes = list()
    for nfile, f in enumerate(rtFiles):
        # start a process related to the present file
        print(f'Starting process {nfile} for file {f}...')
        inPath = f'{options.sample}/{f}'
        tempOutPath = f'{options.output}/tmp/{f}'
        process = multiprocessing.Process(target=processRootFile, args=(nfile, inPath, tempOutPath, runParameters))
        tempOutPaths.append(tempOutPath)
        processes.append(process)
        process.start()

    try:
        for process in processes: process.join()
    except KeyboardInterrupt:
        print('All processes stopped.')
    time.sleep(1)

    ''' after all the processes have finished,
        merge all the temporary output files '''
    mergeRootFiles(tempOutPaths, f'{options.output}/results.root')

    '''residualHistograms = dict()
    efficiencyMatchedHistograms = dict()
    efficiencyPropagatedHistograms = dict()
    efficiencyHistograms = dict()
    recHitHistograms = dict()
    for ieta in ietaRange:
        for 
                key = 'eta%d_fiducialR%1.2e_fiducialPhi%1.2e'%(ieta, fiducialCutR, fiducialCutPhi)
                residualHistograms[key] = rt.TH1F('residuals_%s'%(key), ';#DeltaR#phi (cm);', 100, -0.25, 0.25)
                efficiencyMatchedHistograms[key] = rt.TH1F('ptMatched_%s'%(key), ';p_{t} (GeV);Events', 15, 0, 100)
                efficiencyTotalHistograms[key] = rt.TH1F('ptPropagated_%s'%(key), ';p_{t} (GeV);Events', 15, 0, 100)
                efficiencyHistograms[key] = rt.TH1F('ptEfficiency_%s'%(key), ';p_{t} (GeV);Efficiency', 15, 0, 100)'''
    
if __name__=='__main__': main()
