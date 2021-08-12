#########################################################################
#
# @file
#
# Created:            09-08-2021
#
# Author:             Daniel Oelert (add mail)
#
# Short Description:  Timer class to measure timings.
#
# Last modified:
#
# Copyright:
#
# CMC group, Ruhr-University Bochum 
# The code may only be used and/or copied with the written permission
# of the author or in accordance with the terms and conditions under
# which the program was supplied.  The code is provided "as is"
# without any expressed or implied warranty.
#
#########################################################################

import sys, time

NOT_STARTED = 0
RUNNING = 1
DONE = 2

class Timer(object):
    """Timer object.

    Use like this::

        timer = Timer("Name")
        timer.start()
        # do something
        timer.stop()

    or::

        with timer as mt:
            # do something

    To get a summary call::

        timer.report()

    """

    def __init__(self, desc : str):
        self.desc = desc
        self.status = NOT_STARTED
        self.children = []
        self.t1 = 0
        self.t2 = 0

    @property
    def elapsed(self):
        if self.status == RUNNING:
            return time.time() - self.t1
        return self.t2-self.t1

    def start(self):
        self.status = RUNNING
        self.t1 = time.time()
    
    def stop(self):
        self.t2 = time.time()
        self.status = DONE

    def __enter__(self):
        self.start()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.stop()
        if exc_type:
            raise
        

    def report(self, indent="  ", out=sys.stdout):

        rep = self._report()

        def rep2str(rep,elapsed,level=0):
            
            replist = []

            repstr = "| " + level*indent
            repstr += ("{:<"+str(39-level*len(indent))+"}").format(rep["desc"])
            repstr += " |{:<10}|".format("-"*int(rep["elapsed"]*10/self.elapsed))
            repstr += " {:>7}s  {:5.1f}%  {:5.1f}% |".format(
                "{:.6g}".format(rep["elapsed"])[:7],rep["elapsed"]/elapsed*100,rep["elapsed"]/self.elapsed*100
                )
            # HACK: There is a small error happening when truncating the elapsed time to a width of 7 characters without rounding here

            # HACK: No indication yet for running timers
            # if rep["status"] == RUNNING:
            #     repstr += " running ..."

            replist.append(repstr)

            for r in rep["children"]:
                replist += rep2str(r,rep["elapsed"],level=level+1)

            return replist

        out.write("| Timer report"+39*" "+"| elapsed   rel.%   tot.%  |\n|"+"-"*79+"|\n")
        for i in rep2str(rep,self.elapsed):
            out.write(i+"\n")
        out.write("|"+"-"*79+"|\n")
        return
    
    def _report(self):
        reps = []

        if self.children:
            for ch in self.children:
                reps.append(ch._report())
        return {"status":self.status,"elapsed":self.elapsed,"desc":self.desc,"children":reps}

    def fork(self, desc : str):
        if self.status != RUNNING:
            raise RuntimeError("Unable to fork timer that is not running.")
        new_timer = Timer(desc)
        self.children.append(new_timer)
        return new_timer

