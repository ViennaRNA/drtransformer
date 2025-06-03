#!/usr/bin/env python
#
# drtransformer.plotting
# 
# All plotting utilities should be available through
# the commandline script DrPlotter --help
#

import logging
drlog = logging.getLogger(__name__)

import numpy as np
import matplotlib.pyplot as plt

from . import __version__
from RNA import ptable, loopidx_from_ptable

def plot_nxy(stream, basename, formats,
             title = '',
             plim = 1e-2,
             lines = [],
             xscale = 'log',
             ylim = None,
             xlim = None):
    """ Plot a list of trajectories.

    Args:
      stream: Usually stdin from a nxy file.
      basename:
      formats:
      title (str, optional): Name of the title for the plot.
      plim (float, optional): Minimal occupancy to plot a trajectory. Defaults to 0.01
      lines ([int,..], optional): Selected list of lines to be plotted.
      xscale (str, optional): *lin* or *log*. Default: *log*.
      xlim ((float,float), optional): matplotlib xlim.
      ylim ((float,float), optional): matplotlib ylim.

    Returns:
      [str]: Name of the output file.
    """

    # Prepare the plotting
    fig = plt.figure()
    fig.set_size_inches(7, 3)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xscale(xscale)
    ax.set_ylim([-0.02, 1.02])
    if ylim: ax.set_ylim(ylim)
    if xlim: ax.set_xlim(xlim)

    lines = set(lines)
    nxy = [list(map(float, line.strip().split())) for line in stream if line[0] != '#']

    for e, traject in enumerate(zip(*nxy)):
        if e == 0:
            time = traject
            continue
        if lines and e not in lines:
            continue
        if plim and max(traject) < plim:
            continue
        p, = ax.plot(time, traject, '-', lw = 1.5)
        p.set_label("ID {:d}".format(e))

    ax.set_ylabel('occupancy', fontsize = 11)
    ax.set_xlabel('time [s]', ha = 'center', va = 'center', fontsize = 11)
    plt.legend(facecolor='white', framealpha=1)

    # Add ticks for 1 minute, 1 hour, 1 day, 1 year
    #ax.axvline(x = 60, linewidth = 1, color = 'black', linestyle = '--') # 1 minute
    #ax.axvline(x = 3600, linewidth = 1, color = 'black', linestyle = '--') # 1 hour
    #ax.axvline(x = 86400, linewidth = 1, color = 'black', linestyle = '--') # 1 day
    #ax.axvline(x = 31536000, linewidth = 1, color = 'black', linestyle = '--') # 1 year

    for ending in formats:
        pfile = basename + '.' + ending
        plt.savefig(pfile, bbox_inches = 'tight')
    return 

def plot_simulation(trajectories, basename, formats, 
                    lin_time, log_time, plim = 0.05, tlen = None, 
                    motifs = None, extlen = False, title = ''):
    """DrTransformer standard plotting function.
    """

    fig, axs = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1.8, 1]})
    #fig.set_size_inches(4, 1.95) # paper dims
    fig.set_size_inches(7, 3)
    plt.subplots_adjust(wspace=0)

    [ax1, ax2] = axs
    ax1.set_xscale('linear')
    ax1.spines['right'].set_visible(False)
    ax1.set_xlim((0, lin_time))
    if not extlen:
        ax1.set_ylim([-0.02, 1.02])
    # Tick business
    xticks = [round(x, 2) for x in np.linspace(0, lin_time, 6)]
    ax1.set_xticks(xticks)

    ax2.set_xscale('log')
    ax2.spines['left'].set_visible(False)
    ax2.set_xlim((lin_time, log_time))
    if not extlen:
        ax2.set_ylim([-0.02, 1.02])
    # Tick business
    ax2.yaxis.tick_right()
    ax2.set_yticklabels([])

    ax1t = ax1.twiny()
    ax1t.spines['right'].set_visible(False)
    ax1t.set_xlim(ax1.get_xlim())
    # Tick business
    att, atl = zip(*tlen)
    tl = [int(x) for x in np.linspace(atl[0], atl[-1], 10)]
    tt = [t for e, (t, l) in enumerate(tlen) if l in tl]
    ax1t.set_xticks(tt)
    ax1t.set_xticklabels(tl)
    ax1t.set_xticks(att, minor = True)

    # Mark the end of transcription.
    ax1.axvline(x = lin_time, linewidth = 2, color = 'black', linestyle = '-') 
    ax2.axvline(x = lin_time, linewidth = 2, color = 'black', linestyle = '-') 

    # Set the grid lines.
    ax1.grid(axis = 'y', which = 'major', alpha = 0.7,
             color='gray', linestyle='--', linewidth=0.5)
    ax2.grid(axis = 'y', which = 'major', alpha = 0.7, 
             color='gray', linestyle='--', linewidth=0.5)
    ax1t.grid(which = 'major', alpha = 0.7,
             color='gray', linestyle='--', linewidth=0.5)

    # Plot the data.
    if motifs:
        for i, (name, _) in enumerate(motifs):
            line = trajectories[f'{name}']
            t, o = list(zip(*line))
            (lw, dh, zo) = (2, '-', 1) if max(o) >= plim else (1, '--', 3)
            p, = ax1.plot(t, o, dh, lw=lw, color = f'C{i}', zorder = zo)
            l, = ax2.plot(t, o, dh, lw=lw, color = f'C{i}', zorder = zo)
            l.set_label(f'{name}')

    else:
        for ni in sorted(trajectories):
            line = trajectories[ni]
            t, o = list(zip(*line))
            # Determine which lines are part of the legend:
            # like this, it is only those that are populated
            # at the end of transcription and if they reach
            # an occupancy of 10% or higher
            if t[-1] > lin_time:
                p, = ax1.plot(t, o, '-', lw=1.5)
                l, = ax2.plot(t, o, '-', lw=1.5)
                if max(o) >= 0.1:
                    l.set_label(f"ID {ni}")
            elif max(o) >= 0.1:
                p, = ax1.plot(t, o, '-', lw=1)
                l, = ax2.plot(t, o, '-', lw=1)

    # Legends and labels.
    if extlen:
        ax1.set_ylabel('exterior loop length')
    else:
        ax1.set_ylabel('occupancy')
    ax1.set_xlabel('time (seconds)')
    ax1.xaxis.set_label_coords(0.7, -0.15)
    #ax1.xaxis.set_label_coords(0.7, -0.20) # paper dims
    ax1t.set_xlabel(f'transcript length')
    ax1t.xaxis.set_label_coords(0.7, 1.15)
    #ax1t.xaxis.set_label_coords(0.7, 1.25) # paper dims
    if not extlen:
        ax2.legend(ncol=1, loc = "center right", frameon = True, facecolor = 'white', framealpha = 0.8)
    plt.suptitle(f'{title}', fontsize = 16, y=1.15)
    #plt.suptitle(f'{title}', fontsize = 16, y=1.25) # paper dims


    # Save a file.
    for ending in formats:
        pfile = basename + '.' + ending
        plt.savefig(pfile, bbox_inches = 'tight')
    return

def plot_xmgrace(trajectories, filename):
    """Plot trajectories into a file in grace format.
    """
    head = """
@with line
@line on
@line loctype world
@line g0
@line linewidth 2
@line linestyle 1
@line color 7
@line arrow 0
@line arrow type 0
@line arrow length 1.000000
@line arrow layout 1.000000, 1.000000
@line def
"""
    with open(filename, 'w') as gfh:
        gfh.write(head)
        for nid in sorted(trajectories):
            course = trajectories[nid]
            t, o = list(zip(*course))
            for i in range(len(t)):
                gfh.write("{:f} {:f}\n".format(t[i], o[i]))
            gfh.write("&\n")
    return

def parse_drf(stream):
    xydata = {}
    time_len = []
    llen, lint, ltime, lltime = 0, 0, '0', '0'
    for e, line in enumerate(stream):
        if e == 0:
            continue
        [idx, stime, occu, ss, en] = line.split()
        time = float(stime)
        occu = float(occu)
        xydata[idx] = (xydata.get(idx, [])) + [(time, occu)]
        if len(ss) > llen:
            llen = len(ss)
            lint = float(lltime)
            time_len.append((float(lltime), len(ss)))
        if ltime != stime:
            ltime = stime
            lltime = ltime
    logt = time if time > lint + 60 else lint + 60
    return xydata, lint, logt, time_len

def parse_drf_motifs(stream, motifs):
    xydata = {}
    for (name, m) in motifs:
        xydata[f'{name}'] = [[0, 0]] 
    time_len = []
    llen, lint, ltime, lltime = 0, 0, '0', '0'
    for e, line in enumerate(stream):
        if e == 0:
            continue
        [idx, stime, occu, ss, en] = line.split()
        time = float(stime)
        occu = float(occu)
        for (name, m) in motifs:
            if stime != ltime: # Initialzie all possibilities to 0
                xydata[f'{name}'].append([time, 0])
            if motif_finder(ss, m):
                xydata[f'{name}'][-1][1] += occu
        if len(ss) > llen:
            llen = len(ss)
            lint = float(lltime)
            time_len.append((float(lltime), len(ss)))
        if ltime != stime:
            ltime = stime
            lltime = ltime
    logt = time if time > lint + 60 else lint + 60
    return xydata, lint, logt, time_len

def parse_drf_extlen(stream, maxlen = False):
    xydata = {'exterior': [[0, 0]]}
    time_len = []
    llen, lint, ltime, lltime = 0, 0, '0', '0'
    for e, line in enumerate(stream):
        if e == 0:
            continue
        [idx, stime, occu, ss, en] = line.split()
        time = float(stime)
        occu = float(occu)
        ltable = loopidx_from_ptable(ptable(ss))
        elen = ltable[1:].count(0)
        welen = elen * occu
        if stime != ltime: # Initialzie all possibilities to 0
            xydata['exterior'].append([time, 0])
        if maxlen:
            xydata['exterior'][-1][1] = max(xydata['exterior'][-1][1], elen)
        else:
            xydata['exterior'][-1][1] += welen

        if len(ss) > llen:
            llen = len(ss)
            lint = float(lltime)
            time_len.append((float(lltime), len(ss)))
        if ltime != stime:
            ltime = stime
            lltime = ltime
    logt = time if time > lint + 60 else lint + 60
    return xydata, lint, logt, time_len

def get_motifs(mfile):
    motifs = dict()
    def dbr_to_motif(dbr):
        pt = ptable(dbr)
        motP = [(i, j) for i, j in enumerate(pt) if i and i < j]
        motA = [i for i, j in enumerate(dbr, 1) if j == 'x']
        return (motP, motA)
    with open(mfile, 'r') as mf:
        for line in mf:
            if line[0] == '#':
                continue
            if not line.strip():
                continue
            dbr, name = line.split()[0:2]
            if name in motifs:
                raise SystemExit(f"Motif '{name}' is specified multiple times.")
            motif = dbr_to_motif(dbr) 
            motifs[name] = motif
    return motifs

def motif_finder(ss, motif):
    pt = ptable(ss)
    try:
        return all(pt[i] == j for (i, j) in motif[0]) and all(pt[i] == 0 for i in motif[1])
    except IndexError:
        return False

def get_uprobs(stream):
    uprobs = []
    llen, ltime, lltime = 0, '0', '0'
    up = []
    for i, line in enumerate(stream):
        if i == 0:
            assert line == "id time occupancy structure energy\n"
            continue
        [_, stime, occ, ss, en] = line.split()
        occu = float(occ)
        if ltime == stime and len(ss) == llen:
            for j, b in enumerate(ss):
                 if b == '.': 
                     up[j] += occu
        else:
            #if len(ss) > llen:
            uprobs.append(up)
            up = [0 for _ in range(len(ss))]
            for j, b in enumerate(ss):
                 if b == '.': 
                     up[j] += occu
        llen = len(ss)
        ltime = stime
    uprobs.append(up)
    return uprobs

def main():
    """ DrPlotter -- visulization of cotranscriptional folding simulations.
    """
    import sys, argparse
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'DrPlotter: visualization of the DrForna file format.')

    parser.add_argument('--version', action = 'version', 
            version = '%(prog)s ' + __version__)
    parser.add_argument("-v", "--verbose", action = 'count', default = 0,
            help = """Track process using verbose output.""")
    parser.add_argument("-n", "--name", default = 'DrPlotter', metavar = '<str>',
            help = """Name your output files, name the header of your plots, etc.
            this option overwrites the fasta-header.""")
    parser.add_argument("--molecule", default = '', metavar = '<str>',
            help = """Include molecule name in the plot. """)
    parser.add_argument("--exterior-length", action = "store_true", 
            help = """Plot the mean length of the exterior loop. """)
    parser.add_argument("-f", "--formats", nargs = '+', default = ['pdf'],
            choices = ('pdf', 'svg', 'png', 'gr', 'eps'),
            help = """Plot the simulation using matplotlib (pdf, svg, png)
            and/or write an input file for xmgrace (gr). The legend uses 
            the identities of structures as specifed in the logfile. """)
    parser.add_argument("--motiffile", default = None, metavar = '<str>',
            help = """Specify motifs in a file using dot-bracket notation.""")
    parser.add_argument("-m", "--motifs", nargs = '+', default = None, metavar = '<str>',
            help = """Select specific motif names for plotting. The motifs must be defined in the motiffile.""")
    parser.add_argument("--nxy", action = "store_true", 
            help = "Use nxy format instead of *.drf. This option ignores all other arguments.")
    parser.add_argument("--plim", type = float, default = 0.05,
            help = "Define occupancy limit for legend.")
    parser.add_argument("--accessibilities", action = "store_true", 
            help = """Use heatmaps to show accessibility profiles during transcription (preliminary feature; requires seaborn). """)

    args = parser.parse_args()

    if args.nxy:
        # A hack to plot a simple nxy format,
        mplf = [f for f in args.formats if f != 'gr']
        plot_nxy(sys.stdin, args.name, mplf, title = args.molecule)
        return

    if args.accessibilities:
        import seaborn as sns
        si, sl = 0, 30
        uprobs = get_uprobs(sys.stdin) 
        d1 = len(uprobs)
        d2 = len(uprobs[-1])
        mx = np.zeros((d1, d2))
        mask = np.zeros_like(mx)
        ylabs = []
        for i, l in enumerate(uprobs):
            mx[i][:len(l)] = l
            mask[i][len(l):] = True
            if (not si) and (len(l) == sl):
                si = i
            if si:
                ylabs.append(len(l))
        ax = sns.heatmap(mx[si:], mask = mask[si:], square = True, cmap="viridis", linewidth=0)
        yticks = np.linspace(0, len(ylabs) - 1, 10, dtype=np.int)
        ax.set_yticks(yticks)
        ylabels = [ylabs[idx] for idx in yticks]
        ax.set_yticklabels(ylabels)
        for ending in args.formats:
            pfile = args.name + '.' + ending
            plt.savefig(pfile, bbox_inches = 'tight')
        return

    if args.motifs:
        motifd = get_motifs(args.motiffile)
        for n in args.motifs:
            if n not in motifd:
                raise SystemExit(f"Motif '{n}' not specified.")
        motifs = [(n, motifd[n]) for n in args.motifs]
        xydata, lint, logt, tlen = parse_drf_motifs(sys.stdin, motifs)
    elif args.exterior_length:
        motifs = None
        xydata, lint, logt, tlen = parse_drf_extlen(sys.stdin)
    else:
        motifs = None
        xydata, lint, logt, tlen = parse_drf(sys.stdin)

    if 'gr' in args.formats:
        plot_xmgrace(xydata, args.name + '.gr')

    if not args.molecule and args.name != 'DrPlotter':
        args.molecule = args.name
    mplf = [f for f in args.formats if f != 'gr']
    if mplf:
        with plt.rc_context({
            'axes.spines.right': True,
            'axes.spines.top': True,
            'axes.spines.left': True,
            'axes.spines.bottom': True,
            'axes.grid': False,
            'xtick.direction': 'out',
            'ytick.direction': 'out',
            'xtick.major.size': 5,
            'ytick.major.size': 5,
            'font.size': 12,
            'axes.titlesize': 14,
            'axes.labelsize': 12}):
            plot_simulation(xydata, args.name, mplf, lint, logt, args.plim, 
                            tlen = tlen, motifs = motifs, 
                            extlen = args.exterior_length,
                            title = args.molecule)


if __name__ == '__main__':
    main()
