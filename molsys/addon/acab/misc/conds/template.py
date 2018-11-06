"""
constraint template for pyscipopt
"""

class NewConshdlr(Conshdlr):

    def conscheck(self, *args):
        #if self.findSubtours(checkonly = True, sol = solution):
        if not self.addcut(checkonly = True, sol = solution):
            return {"result": SCIP_RESULT.INFEASIBLE}
        else:
            return {"result": SCIP_RESULT.FEASIBLE}

    def consenfolp(self, *args):
        #if self.findSubtours(checkonly = False, sol = None):
        if self.addCuts(checkonly = False):
            return {"result": SCIP_RESULT.CONSADDED}
        else:
            return {"result": SCIP_RESULT.FEASIBLE}


    def conslock(self, *args):
        pass

if __name__ == "__main__":
    model.init()
    conshdlr = Conshdlr_type()
    model.includeConshdlr(conshdlr, "VRP", "VRP constraint handler",
        sepapriority = 0, # default
        enfopriority = 1,
        chckpriority = 1,
        sepafreq = -1, # default
        propfreq = -1, # default
        eagerfreq = -1,
        maxprerounds = 0,
        delaysepa = False, # default
        delayprop = False, # default
        needscons = False,
        presoltiming = SCIP_PRESOLTIMING.FAST,
        proptiming = SCIP_PROPTIMING.BEFORELP) # default
    model.includeConshdlr(conshdlr, "TSP", "TSP subtour eliminator",
        sepapriority = -1,
        enfopriority = -1,
        chckpriority = -1,
        sepafreq = -1, # default
        propfreq = -1, # default
        eagerfreq = -1,
        maxprerounds = 0,
        delaysepa = False, # default
        delayprop = False, # default
        needscons = False,
        presoltiming = SCIP_PRESOLTIMING.FAST,
        proptiming = SCIP_PROPTIMING.BEFORELP) # default
    model.optimize()
