class Hit:
    # generic hit, can be reco or propagated
    def __init__(self, event, index):
        self.event = event # event in the sample tree
        self.index = index # index of the hit in the reco or prop branches

class RecHit(Hit):
    @property
    def region(self): return self.event.gemRecHit_region[self.index]

    @property
    def chamber(self): return self.event.gemRecHit_chamber[self.index]

    @property
    def ieta(self): return self.event.gemRecHit_etaPartition[self.index]

    @property
    def layer(self): return self.event.gemRecHit_layer[self.index]

    @property
    def phi(self): return self.event.gemRecHit_g_phi[self.index]

    @property
    def r(self): return self.event.gemRecHit_g_r[self.index]

    @property
    def globalX(self): return self.event.gemRecHit_g_x[self.index]

    @property
    def globalY(self): return self.event.gemRecHit_g_y[self.index]

    '''def matches(self, propHit, drphiMax):
        matchesRegion = self.region==propHit.region
        matchesLayer = self.layer==propHit.layer
        matchesChamber = self.chamber==propHit.chamber
        #matchesEta = self.ieta==propHit.ieta
        matchesEta = abs(self.ieta-propHit.ieta) < 2
        matchesDRPhi = abs(self.phi-propHit.phi)*propHit.r<=drphiMax
        matches = matchesRegion and matchesChamber and matchesEta and matchesLayer and matchesDRPhi
        if matches: propHit.matchFound = matches
        return matches'''
    def matches(self, propHit, drphiMax):
        # check region, layer, chamber, eta, drphi
        if not self.region==propHit.region: return False
        if not self.layer==propHit.layer: return False
        if not self.chamber==propHit.chamber: return False
        if not self.ieta==propHit.ieta: return False
        # alternatively to allow neighbour eta partitions: if not abs(self.ieta-propHit.ieta) < 2: return False
        if not abs(self.residual(propHit))<=drphiMax: return False
        # if you arrive here then everything matches:
        propHit.matchFound = True
        return True
    
    def residual(self, propHit):
        return (self.phi-propHit.phi)*propHit.r

class PropHit(Hit):
    def __init__(self, event, index):
        Hit.__init__(self, event, index)
        self.matchFound = False

    @property
    def matchFound(self): return self._matchFound

    @matchFound.setter
    def matchFound(self, m): self._matchFound = m

    @property
    def checksFiducialCuts(self):
        return True # to be implemented

    @property
    def isGEM(self): return self.event.m_propagated_isGEM[self.index]

    @property
    def pt(self): return self.event.mu_propagated_pt[self.index]

    @property
    def region(self): return self.event.mu_propagated_region[self.index]

    @property
    def chamber(self): return self.event.mu_propagated_chamber[self.index]

    @property
    def ieta(self): return self.event.mu_propagated_etaP[self.index]

    @property
    def layer(self): return self.event.mu_propagated_layer[self.index]

    @property
    def phi(self): return self.event.mu_propagatedGlb_phi[self.index]

    @property
    def r(self): return self.event.mu_propagatedGlb_r[self.index]

    @property
    def globalX(self): return self.event.mu_propagatedGlb_x[self.index]

    @property
    def globalY(self): return self.event.mu_propagatedGlb_y[self.index]

    @property
    def isME11(self):
        return self.event.mu_propagated_isME11[self.index]
