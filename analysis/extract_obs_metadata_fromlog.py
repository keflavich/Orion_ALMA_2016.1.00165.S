import paths
import datetime
import glob
import os
from astropy.table import Table

def get_ant_pos(fn):
    ants = {}
    with open(fn, 'r') as fh:
        start = False
        for line in fh.readlines():
            if start and 'East' not in line and 'Station' not in line:
                ls = line.split()
                ID, east, north = ls[0],ls[7],ls[8]
                size = ls[3]
                ants[ID] = (east,north)
            if 'Antennas:' in line:
                start = True

    return ants, size

def ants_to_baselines(ants):
    bls = {}
    for ant1 in ants:
        for ant2 in ants:
            if ant1==ant2:
                continue
            else:
                x1,y1 = float(ants[ant1][0]), float(ants[ant1][1])
                x2,y2 = float(ants[ant2][0]), float(ants[ant2][1])
                d = ((x1-x2)**2 + (y1-y2)**2)**0.5
                bls[ant1+"-"+ant2] = d
    return bls


def longest_baseline(fn):
    return max(ants_to_baselines(get_ant_pos(fn)[0]).values())
def shortest_baseline(fn):
    return min(ants_to_baselines(get_ant_pos(fn)[0]).values())

def get_freqtable(fn):
    with open(fn, 'r') as fh:
        lines = fh.readlines()
    return get_freqs(lines)

def get_freqs(lines):
    start = False
    data = []
    for line in lines:
        if 'Spectral Windows:' in line:
            start = True
            continue
        if 'Sources:' in line:
            start = False
            break
        if start:
            if 'SpwID' in line:
                header = line.split()
            else:
                data.append(line.split())

    if header[-2] != 'Num':
        raise ValueError
    header[-3] = 'BBC Num'
    header[-2] = 'Corrs'
    header.pop(-1)

    for row in data:
        # XX YY -> XXYY
        if row[-1] == 'YY' and row[-2] == 'XX':
            row[-2] = row[-2]+row[-1]
            row.pop(-1)

    data_T = list(map(list, zip(*data)))

    tbl = Table(names=header, data=data_T,
                dtype=[int, str, int, str, float, float, float, float, int, str])

    return tbl

def get_date(fn):
    with open(fn, 'r') as fh:
        lines = fh.readlines()
        for line in lines:
            if 'Observed from' in line:
                dateline = line
                break
        start_date = dateline.split()[2].split("/")[0]
        end_date = dateline.split()[4].split("/")[0]

    return start_date, end_date

def get_total_time(fn):
    with open(fn, 'r') as fh:
        lines = fh.readlines()
        for line in lines:
            if 'Total elapsed time' in line:
                timeline = line
                break
        time = timeline.split()[-2]
    return time

def get_metadata_line(fn, DEBUG=False):
    ants,size = get_ant_pos(fn)
    start_date, end_date = get_date(fn)
    baselines = ants_to_baselines(ants)
    longestbl = max(baselines.values())
    shortestbl = min(baselines.values())
    integrationtime = get_total_time(fn)

    frqtbl = get_freqtable(fn)
    # http://www.eso.org/public/usa/teles-instr/alma/receiver-bands/
    avfrq = frqtbl['Ch0(MHz)'].mean() / 1e3
    if (avfrq > 84) & (avfrq < 116):
        band = 3
    elif (avfrq > 125) & (avfrq < 163):
        band = 4
    elif (avfrq > 163) & (avfrq < 211):
        band = 5
    elif (avfrq > 211) & (avfrq < 275):
        band = 6
    elif (avfrq > 275) & (avfrq < 373):
        band = 7
    elif (avfrq > 385) & (avfrq < 500):
        band = 8
    elif (avfrq > 602) & (avfrq < 720):
        band = 9
    elif (avfrq > 787) & (avfrq < 950):
        band = 10

    if start_date != end_date:
        # skip because we don't want to include the merged, only the original, MS metadata
        return None
    
    line = ("{0} & {6} & {1}m & {2} & {3}-{4} & {5}\\\\"
            .format(start_date, int(round(float(size))),
                    int(round(float(integrationtime))),
                    int(round(shortestbl)),
                    int(round(longestbl)),
                    len(ants),
                    band,
                   )
           )
    if DEBUG:
        print(line)
    return line, (start_date, size, integrationtime, shortestbl, longestbl, len(ants))

def make_meta_tables(listobspath=os.path.join(paths.reductionpath, 'listobs'), DEBUG=False):
    listfiles=glob.glob(os.path.join(listobspath, '*.listobs'))

    lines = [get_metadata_line(fn, DEBUG=DEBUG) for fn in listfiles]
    lines = [l for l in lines if l is not None]
    print()
    print()

    dct = {}
    dct1 = {}
    for prtline,(date, size, integrationtime, shortestbl, longestbl, nants) in lines:
        if date in dct:
            if DEBUG:
                print("old integration time for {0} = {1}".format(date, integrationtime))
                print("will add {0}".format(dct1[date][2]))
                print("new sum is {0}".format(float(dct1[date][2])+float(integrationtime)))
            integrationtime = float(dct1[date][2]) + float(integrationtime)
            if DEBUG:
                print("new integration time for {0} = {1}".format(date, integrationtime))
            dct[date] = "{0} & {1}m & {2} & {3}-{4} & {5}\\\\".format(date,
                                                                      int(round(float(size))),
                                                                      int(round(float(integrationtime))),
                                                                      int(round(shortestbl)),
                                                                      int(round(longestbl)),
                                                                      nants)
            dct1[date] = [date,size,integrationtime,shortestbl,longestbl]
        else:
            dct[date] = prtline
            dct1[date] = [date,size,integrationtime,shortestbl,longestbl]


    for k,v in sorted(dct.items(), key=lambda x: datetime.datetime.strptime(x[0][:11],'%d-%b-%Y')):
        print(v)

    return lines

if __name__ == "__main__":
    formatted, data = zip(*make_meta_tables(DEBUG=False))
    dates = [datetime.datetime.strptime(x[0], '%d-%b-%Y') for x in data]

    lines = "\n".join([x for _,x in sorted(set(zip(dates,formatted)))])
    print()
    print(lines)

    basetable = r"""
\begin{table*}[htp]
\centering
\caption{Observation Summary}
\begin{tabular}{llllll}
\label{tab:observations}
Date & Band & Array & Observation Duration &  Baseline Length Range  & \# of antennae\\
     &      &       & seconds              & meters                    & \\
\hline
DATA
\hline
\end{tabular}
\end{table*}
"""

    with open(paths.texpath('obs_metadata.tex'), 'w') as fh:
        fh.write(basetable.replace("DATA", lines))
