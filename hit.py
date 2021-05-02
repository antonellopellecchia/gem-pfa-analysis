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

    def matches(self, propHit, drphiMax):
        matchesRegion = self.region==propHit.region
        matchesChamber = self.chamber==propHit.chamber
        matchesEta = self.ieta==propHit.ieta
        matchesLayer = self.layer==propHit.layer
        matchesDRPhi = abs(self.phi-propHit.phi)*self.r<=drphiMax
        return matchesRegion and matchesChamber and matchesEta and matchesLayer and matchesDRPhi
    
    def residual(self, propHit):
        return (self.phi-propHit.phi)*self.r

class PropHit(Hit):
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