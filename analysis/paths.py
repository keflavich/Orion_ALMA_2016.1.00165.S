import os

root = os.path.expanduser('~/work/orion/alma_lb/')

datapath = os.path.join(root, 'FITS/')
regpath = os.path.join(root, 'regions/')
figurepath = os.path.join(root, 'figures/')
redpath = os.path.join(root, 'reduction/')

def dpath(x, datapath=datapath):
    return os.path.join(datapath, x)

def rpath(x, regpath=regpath):
    return os.path.join(regpath, x)

def fpath(x, figurepath=figurepath):
    return os.path.join(figurepath, x)

def redpath(x, redpath=redpath):
    return os.path.join(redpath, x)
