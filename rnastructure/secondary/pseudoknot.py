from subprocess import Popen

from rnastructure.secondary.connect import Parser as Connect

class Remove(object):
    def __init__(self, filename=None):
        self.name = filename
        self.program = "RemovePseudoknots"

    def remove(self, structure, options={}):
        args = [self.program]
        args.append(in_ct)
        args.append(out_ct)
        for key, value in options:
            flag = "--%s"
            if len(key) == 1:
                flag = "-%s"
            replace = value or ''
            flag = flag % replace
            args.extend([flag, value])
        Popen(args)
        out_file = open(out_ct)
        removed = Connect(out_file)
        out_file.close()
        return removed
