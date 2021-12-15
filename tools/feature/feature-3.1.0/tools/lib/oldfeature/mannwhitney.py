#! /bin/env python

# Implements the Mann-Whitney Non-Parametric test of difference in means
#
# Reference:
#   Mathematical Statistics and Data Analysis / John A. Rice
#   2nd ed. 1995

from math import sqrt,exp,pi
from mystat import normal_cdf


# ====================================================================
def sum(aList):
    return reduce(lambda a,b:a+b, aList)


# ====================================================================
def mean(aList):
    return sum(aList)/float(len(aList))


# ====================================================================
def ranker(small, large):
    # Combine the lists and sort
    rankList = small+large
    rankList.sort()

    # Assign ranks for each value
    rankMapList = {}
    for index in range(len(rankList)):
        rankMapList.setdefault(rankList[index],[]).append(index+1)

    # Break ties
    rankMap = {}
    for (value, ranks) in rankMapList.items():
        rankMap[value] = mean(ranks)

    # Return ranks for each list
    smallRank = [rankMap[value] for value in small]
    largeRank = [rankMap[value] for value in large]
    
    return (smallRank, largeRank)


# ====================================================================
zmap = { 0.20: 1.28,
         0.10: 1.64,
         0.05: 1.96,
         0.01: 2.58 }


# ====================================================================
def _simpson(f,a,b,n):
    """ Simpson's Rule """
    m = int(n/2)
    h = float(b-a)/(m*2)
    Sn = f(a) + 4.0*f(b-h) + f(b)
    for k in range(1,m):
        Sn += 2.0*f(a+2*k*h)
        Sn += 4.0*f(a+2*k*h-h)
    Sn *= h/3.0
    return Sn


# ====================================================================
def _normal_pdf(x):
    return exp(-x**2/2.0)/sqrt(2.0*pi)


# ====================================================================
def _normal_cdf(z):
    n = 100
    if z < 0:
        return 0.5-simpson(normal_pdf,0,-z,n)
    return 0.5+simpson(normal_pdf,0,z,n)


# ====================================================================
def ranksumtest(listA, listB, alpha):
    # make sure alpha is supported
    if alpha not in zmap.keys():
        return None
    z = zmap[alpha]
    
    # determine larger/smaller lists
    if len(listA) > len(listB):
        (small,large) = (listB,listA)
    else:
        (small,large) = (listA,listB)

    # calculate ranks
    (rankSmall, rankLarge) = ranker(small, large)
    n = float(len(small))
    m = float(len(large))

    # calculate test statistic
    R = sum(rankSmall)
    Rprime = n*(m+n+1) - R
    Rstar = min(R,Rprime)

    # calculate critical value
    critval = n*(n+m+1)/2.0 - z*sqrt(n*m*(n+m+1)/12.0)

    # if rstar is less than critical value, reject NULL hypothesis
    if Rstar < critval:
        return 1
    return 0


# ====================================================================
def ranksumtestAppx(listA, listB, alpha, optret=None):
    if not listA or not listB:
        return None
    
    # determine larger/smaller lists
    if len(listA) > len(listB):
        (small,large) = (listB,listA)
    else:
        (small,large) = (listA,listB)
    (small,large) = (listA, listB)

    # Calculate Ranks
    (rankSmall, rankLarge) = ranker(small,large)
    n = float(len(small))
    m = float(len(large))

    # Calculate Test Statistic
    Ty = sum(rankSmall)
    ETy = n*(n+m+1)/2.0
    VarTy = m*n*(m+n+1)/12.0
    zT = (Ty-ETy)/sqrt(VarTy)

    # calculate p-value
    p = 2*(1-normal_cdf(abs(zT)))    # two-tailed test

    # Return optional parameters
    if optret != None:
        optret.extend((zT, p, Ty, ETy, VarTy))

    if p < alpha:
        return 1
    return 0


# ====================================================================
if __name__ == "__main__":
    import sys
    args = sys.argv[1:]
    if len(args) != 3:
        print "Usage: %s <list1> <list2> <alpha>" % sys.argv[0]
        sys.exit(2)

    list1 = [float(line) for line in open(args[0]).readlines()]
    list2 = [float(line) for line in open(args[1]).readlines()]
    alpha = float(args[2])

    if ranksumtestAppx(list1,list2,alpha):
        print "yes"
    else:
        print "no"
