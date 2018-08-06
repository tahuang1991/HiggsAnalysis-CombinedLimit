from HiggsAnalysis.CombinedLimit.PhysicsModel import PhysicsModel

class DYDataDrivenEstimationModel(PhysicsModel):
    def __init__(self):
        PhysicsModel.__init__(self)
	#self.signal_xsec = 10.0

    def doParametersOfInterest(self):
        PhysicsModel.doParametersOfInterest(self)

        self.modelBuilder.doVar("MinusOne[-1]")
        #self.modelBuilder.doVar("Signal_xsec[%f]"%self.signal_xsec)

    def getYieldScale(self, bin, process):
        if self.DC.isSignal[process]:
	    print "signal process ",process
            return "r"

        if "nobtag_to_btagM" in process or "untagged" in process:
	    print "nobtag_to_btagM, process ",process
            if 'data' in process:
                return 1
            else:
                return "MinusOne"

        return 1

DYDataDrivenEstimationModelInstance = DYDataDrivenEstimationModel()
