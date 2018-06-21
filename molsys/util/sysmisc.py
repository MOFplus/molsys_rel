import os
import glob

def _makedirs(folder):
    """
    standard os.makedirs with error handler in case directory
    aldready exists
    """
    try:
        os.makedirs(folder)
    except OSError:
        pass

def _checkrundir(folder, basename):
    globs = glob.glob("%s%s[0-9]*_%s%s" % (folder, os.sep, basename, os.sep))
    splits = [''.join(g.split("_%s%s" % (basename,os.sep))[:-1]) for g in globs]
    splits = [''.join(s.split("%s%s" % (folder, os.sep))[-1]) for s in splits]
    ints = []
    for split in splits:
        try:
            i = int(split)
            ints.append(i)
        except ValueError:
            pass
    if ints:
        inew = max(ints) + 1
    else:
        inew = 0
    os.makedirs("%s%s%s_%s%s" % (folder, os.sep, inew, basename, os.sep))
    return "%s%s%s_%s%s" % (folder, os.sep, inew, basename, os.sep)
