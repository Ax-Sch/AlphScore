import dlmio
import time

help(dlmio.dlmread)
filename = '/home/mliang/Projects/Research/DataSetGenerator/astral_40.serog.features'
#filename = '/home/mliang/Projects/Research/Selector/trypsin_ser_40_og_site.features'
infh = file(filename)
starttime = time.time()
x = dlmio.dlmread(filename)
print x.shape
print time.time()-starttime
starttime = time.time()
y = dlmio.dlmread(infh,delim='\t',row=100,col=100)
print y.shape
print time.time()-starttime
