from HiggsAnalysis.CombinedLimit.PhysicsModel import PhysicsModel

class DYDataDrivenEstimationModel(PhysicsModel):
    def __init__(self):
        PhysicsModel.__init__(self)

    def doParametersOfInterest(self):
        PhysicsModel.doParametersOfInterest(self)

        self.modelBuilder.doVar("MinusOne[-1]")

    def getYieldScale(self, bin, process):
        if self.DC.isSignal[process]:
            return "r"

        if "nobtag_to_btagM" in process or "untagged" in process:
            if 'data' in process:
                return 1
            else:
                return "MinusOne"

        return 1

DYDataDrivenEstimationModelInstance = DYDataDrivenEstimationModel()
