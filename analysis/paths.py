import os
import socket

if 'nmpost' in socket.gethostname():
    root = '/lustre/aginsbur/orion/2016.1.00165.S/'
else:
    root = os.path.expanduser('~/work/orion/alma_lb/')

datapath = os.path.join(root, 'FITS/')
regpath = os.path.join(root, 'regions/')
figurepath = os.path.join(root, 'figures/')
reductionpath = os.path.join(root, 'reduction/')
paperpath = os.path.join(root, 'paper_sourceImass')
paperfigpath = os.path.join(root, 'paper_sourceImass/figures/')



def dpath(x, datapath=datapath):
    return os.path.join(datapath, x)

def rpath(x, regpath=regpath):
    return os.path.join(regpath, x)

def fpath(x, figurepath=figurepath):
    return os.path.join(figurepath, x)

def redpath(x, redpath=reductionpath):
    return os.path.join(redpath, x)

def texpath(x, paperpath=paperpath):
    return os.path.join(paperpath, x)

def texfpath(x, paperfigpath=paperfigpath):
    return os.path.join(paperfigpath, x)
