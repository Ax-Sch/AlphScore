#! /usr/bin/env python
import agrepy

class AgrepPattern:
    def __init__(self,pattern,edit=1):
        self.pattern = pattern
        self.pattern_obj = agrepy.compile(pattern,len(pattern),edit)
    def search(self,line,gotoend=0):
        return agrepy.agrepy(self.pattern,len(self.pattern),line,len(line),gotoend,self.pattern_obj)

def agrep(line,pattern,edit=1,gotoend=0):
    return AgrepPattern(pattern,edit).search(line,gotoend)
