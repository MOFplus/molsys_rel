import time

from molsys.util.timer import Timer

maintimer = Timer("Test timer")
with maintimer as mt:
    time.sleep(6)
    with mt.fork("sub timer 1") as subtimer1:
        time.sleep(0.001)
    with mt.fork("sub timer 2") as subtimer2:
        time.sleep(2)
        with subtimer2.fork("subsub timer") as subsubtimer:
            time.sleep(1)
    time.sleep(1)
    with mt.fork("sub timer 3") as subtimer3:
        time.sleep(1)

maintimer.report()

