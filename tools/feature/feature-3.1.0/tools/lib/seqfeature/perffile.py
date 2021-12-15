#! /usr/bin/env python

class PerformanceFile:
    class Entry:
        def __init__(self,line):
            fields = line.split()
            self.score = float(fields[0])
            self.fpr = float(fields[1])
            self.tpr = float(fields[2])
            self.posCount = int(fields[3])
            self.negCount = int(fields[4])

    def __init__(self,filename):
        self.Load(filename)

    def Load(self,filename):
        infh = file(filename)
        line = infh.next()

        fields = line.split()
        self.auc = float(fields[1])
        self.numPos = int(fields[2])
        self.numNeg = int(fields[3])

        self.entries = []
        for line in infh:
            entry = self.Entry(line)
            self.entries.append(entry)

    def GetAuc(self):
        return self.auc
