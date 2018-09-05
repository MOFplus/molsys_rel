import os
import glob

def _makedirs(folder):
    """
    Make leaf directory and intermediates if it does not exists. Else:do nothing
    It is the standard os.makedirs with error handler in case directory aldready
    exists
    """
    try:
        os.makedirs(folder)
    except OSError:
        pass

def _checkrundir(folder, basename):
    """Check if numbered basename directory in folder exists; create a new one.
    The function searches for directories in the folder with a starting
    indices (e.g. `0_run` where `0` is the index and `run` the basename).
    If no directory is found, a new `0_$basename` directory is created.
    Else: the maximum index in the found directories is increased and a new
    directory with the increased index is created.
    The function is intended for output directories after multiple run in the
    same folder."""
    # find directories
    globs = glob.glob("%s%s[0-9]*_%s%s" % (folder, os.sep, basename, os.sep))
    # split removing basename
    splits = [''.join(g.split("_%s%s" % (basename,os.sep))[:-1]) for g in globs]
    # split removing folder
    splits = [''.join(s.split("%s%s" % (folder, os.sep))[-1]) for s in splits]
    ints = []
    # only integers are stored
    for split in splits:
        try:
            i = int(split)
            ints.append(i)
        except ValueError:
            pass
    if ints: # if at least one directory with an index is found, increase max
        inew = max(ints) + 1
    else: # take 0
        inew = 0
    # create new directory
    # N.B.: now os.makedirs cannot raise OSError since the directory is new
    os.makedirs("%s%s%s_%s%s" % (folder, os.sep, inew, basename, os.sep))
    # return directory name
    return "%s%s%s_%s%s" % (folder, os.sep, inew, basename, os.sep)
