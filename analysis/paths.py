import os
import socket

if 'nmpost' in socket.gethostname():
    root = '/lustre/aginsbur/orion/2016.1.00165.S/'
    imagingpath = os.path.join(root, 'imaging/')
elif 'rng9000' in socket.gethostname():
    root = '/home/rng90003/orion/2016.1.00165.S'
    imagingpath = os.path.join(root, 'imaging/')
elif 'ufhpc' in socket.gethostname():
    root = '/orange/adamginsburg/orion/2016.1.00165.S'
    imagingpath = os.path.join(root, 'imaging/')    
else:
    root = os.path.expanduser('~/work/orion/alma_lb/')
    imagingpath = fullcubepath = '/Volumes/external/orion/'

datapath = os.path.join(root, 'FITS/')
regpath = os.path.join(root, 'regions/')
figurepath = os.path.join(root, 'figures/')
reductionpath = os.path.join(root, 'reduction/')
paperpath = os.path.join(root, 'paper_sourceImass')
paperfigpath = os.path.join(root, 'paper_sourceImass/figures/')
paperpath2 = os.path.join(root, 'paper_salts')
tablepath = os.path.join(root, 'tables')
saltpath = os.path.join(root, 'salt_data')



def dpath(x, datapath=datapath):
    return os.path.join(datapath, x)

def rpath(x, regpath=regpath):
    return os.path.join(regpath, x)

def fpath(x, figurepath=figurepath):
    return os.path.join(figurepath, x)

def tpath(x, figurepath=tablepath):
    return os.path.join(tablepath, x)

def redpath(x, redpath=reductionpath):
    return os.path.join(redpath, x)

def texpath(x, paperpath=paperpath):
    return os.path.join(paperpath, x)

def texpath2(x, paperpath=paperpath2):
    return os.path.join(paperpath, x)

def texfpath(x, paperfigpath=paperfigpath):
    return os.path.join(paperfigpath, x)

def imagingpath(x, imagingpath=imagingpath):
    return os.path.join(imagingpath, x)

def fcp(x):
    # "full cube path"
    return imagingpath(x)

def salty(x, path=saltpath):
    return os.path.join(path, x)
