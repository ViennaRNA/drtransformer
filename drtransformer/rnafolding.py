#
# drtransformer.rnafolding
# 
# Utilities using the ViennaRNA package and more.
#

import logging
drlog = logging.getLogger(__name__)

import math
import numpy as np
from itertools import chain, combinations

import RNA

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Custom error definitions                                                     #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
class RNAfoldingError(Exception):
    pass

def parse_model_details(parser):
    """ ViennaRNA Model Details Argument Parser.  """
    model = parser.add_argument_group('ViennaRNA model details')

    model.add_argument("-T", "--temp", type = float, default = 37.0, 
        metavar = '<flt>',
        help = 'Rescale energy parameters to a temperature of temp C.')

    model.add_argument("-4", "--noTetra", action = "store_true",
        help = """Do not include special tabulated stabilizing 
        energies for tri-, tetra- and hexaloop hairpins.""")

    model.add_argument("-d", "--dangles", type = int, default = 2, 
        metavar = '<int>',
        help = """How to treat dangling end energies for bases adjacent to
        helices in free ends and multi-loops.""")

    model.add_argument("--noGU", action = "store_true",
        help = 'Do not allow GU/GT pairs.')

    model.add_argument("--noClosingGU", action = "store_true",
        help = 'Do not allow GU/GT pairs at the end of helices.')

    model.add_argument("--noLP", action = "store_true",
        help = ("Only consider structures without lonely pairs. "
                "(This is never enforced when estimating energy barriers.)"))

    model.add_argument("-P", "--paramFile", action = "store", default = None,
        metavar = '<str>',
        help = """Read energy parameters from paramfile, instead of 
        using the default parameter set.""")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Find fraying neighbors (with caching)                                        #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
def find_fraying_neighbors(seq, md, parents, mfree = 6):
    """ Generates all secondary structures with missing fraying helices.

    Args:
        seq (str): RNA sequence.
        md (RNA.md()): ViennaRNA model details.
        parents (list[str]): A list of parent secondary structures.
        mfree (int, optional): Minimum number of freed bases to stop recursion.

    Returns:
        dict: new_nodes[parent]=set(neigbors)
    """
    assert len(parents[0]) == len(seq), "length of sequence and parents doesn't match."
    ecache = dict()
    new_nodes = dict() 
    for pss in parents:
        nbrs = set()
        # NOTE: potential optimization via passing pair tables.
        for oss in open_fraying_helices(seq, pss, mfree):
            nss, ecache = fold_exterior_loop(seq, oss, md, ecache)
            nbrs.add(nss)
        nbrs.discard(pss)
        new_nodes[pss] = nbrs
    return new_nodes

def open_fraying_helices(seq, ss, mfree = 6):
    """ Generate structures with opened fraying helices. 
    
    Fraying helices share a base-pair with an exterior loop region. The
    function returns n structures, where the first n-1 correspond to individual
    fraying helices opened respectively and structure n has all fraying helices
    opened at once.

    Args:
        seq (str): Primary structure.
        ss (str): Secondary structure.
        mfree (int, optional): Minimal number of bases freed by a fraying helix
            move. If less bases are freed, and there exists a nested helix, then
            that helix is opened as well. Defaults to 6.
    """
    # mutable pairtable and secondary structure
    pt = list(RNA.ptable(ss))
    nbr = list(ss)

    def gen_nbrs(ss, mb, pt, myrange, mfree):
        """Recursive helix opening.
    
        Args:
            ss (str): Reference secondary structure.
            mb (list): A mutable version of ss, to open *all* fraying helices.
            pt (list): ViennaRNA pair table.
            myrange (int, int): The currently interesting range (n, m) of pt.
            mfree: The number of bases to be freed.
    
        Yields:
            (str): secondary structures without fraying helices.
        """
        (n, m) = myrange
        skip = 0  # fast forward in case we have deleted stuff
        for i in range(n, m):
            j = pt[i + 1] - 1
            if j == -1 or i < skip:
                continue
    
            nb = list(ss)
            [o, l] = [0, 0]
            [p, q] = [i, j]
    
            while p < q and (l == 0 or o < mfree):
                pp, qq = p+1, q+1
                if pt[pp] != qq or pp != pt[qq]:
                    # it is a multiloop.
                    for nb in gen_nbrs(''.join(nb), mb, pt, (p, q), mfree - o):
                        yield nb
                    break
                # remove the base-pairs
                pt[pp] = pt[qq] = 0
                nb[p] = nb[q] = '.'
                mb[p] = mb[q] = '.'
                o += 2  # one base-pair deleted, two bases freed
    
                l = 0  # reset interior-loop size
                while (p < q and pt[pp + 1] == 0):
                    p, pp, l = p + 1, pp + 1, l + 1
                p += 1
                while (p < q and pt[qq - 1] == 0):
                    q, qq, l = q - 1, qq - 1, l + 1
                q -= 1
                o += l
            else: # it was *not* a multiloop.
                yield ''.join(nb)
            skip = j + 1
        return

    count = 0
    for nb in gen_nbrs(ss, nbr, pt, (0, len(ss)), mfree):
        count += 1
        yield nb
    if count > 1:
        yield ''.join(nbr)
    return

def fold_exterior_loop(seq, con, md, cache = None, spacer = 'NNN'):
    """ Constrained folding of the exterior loop.

    All constrained helices are replaced with the motif:
        NNNNNNN
        ((xxx))
    for example a helix with the closing-stack CG-UG:
        CG ~ UG -> CGNNNUG
        (( ~ )) -> ((xxx))
    This reduces the sequence length (n) and therefore the runtime O(n^3),
    and it enables the identification of independent structures with the same
    exterior loop features.

    Args:
        seq (str):          RNA sequence
        con (str):          RNA structure constraint
        md (RNA.md()):      ViennaRNA model details (temperature, noLP, etc.)
        spacer (str, optional): Specify the sequence of the spacer region.
                            Defaults to 'NNN'.

    Returns:
      (str, str):
    """
    pt = RNA.ptable(con)
    ext_seq = ''
    ext_con = ''

    # shrink the sequences
    skip = 0
    for ii, jj in enumerate(pt[1:], 1):
        i, j = ii-1, jj-1 # pt indices to string indices
        if i < skip:
            continue
        if j == -1:
            ext_seq += seq[i]
            ext_con += '.'
        else:
            ext_seq += seq[i] + seq[i + 1]
            ext_seq += spacer
            ext_seq += seq[j - 1] + seq[j]
            ext_con += '(('
            ext_con += 'x' * len(spacer)
            ext_con += '))'
            skip = j + 1

    # Check cache if we have seen this exterior loop before.
    if cache is None:
        cache = {}
    if ext_seq in cache:
        css = cache[ext_seq]
    else:
        fc_tmp = RNA.fold_compound(ext_seq, md)
        fc_tmp.constraints_add(ext_con, 
                               RNA.CONSTRAINT_DB_DEFAULT | 
                               RNA.CONSTRAINT_DB_ENFORCE_BP)
        css, cfe = fc_tmp.mfe()
        cache[ext_seq] = css
        del fc_tmp

    # replace characters in constraint
    c, skip = 0, 0
    for ii, jj in enumerate(pt[1:], 1):
        i, j = ii-1, jj-1 # pt indices to string indices
        if i < skip:
            continue
        if j == -1:
            con = con[:i] + css[c] + con[i + 1:]
            c += 1
        else:
            c += len(spacer) + 4
            skip = j + 1

    return con, cache

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# More constrained folding utils                                               #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
def mfe_constrained(seq, md, dbcon):
    """Constrained folding (base-pairs are not enforced)."""
    fc = RNA.fold_compound(seq, md)
    fc.constraints_add(dbcon, RNA.CONSTRAINT_DB_DEFAULT)
    mss, mfe = fc.mfe()
    del fc
    return mss, int(round(mfe * 100))

def forbid_all_basepairs(seq, fc):
    """Modify fold compound to disable all base-pairs."""
    for i in range(1, len(seq) + 1):
        for j in range(i + 4, len(seq) + 1):
            fc.hc_add_bp(i, j, RNA.CONSTRAINT_CONTEXT_NONE | 
                               RNA.CONSTRAINT_CONTEXT_NO_REMOVE)
    return fc

def get_basepairs(dotbrackets):
    bps = set()
    for db in dotbrackets:
        pt = RNA.ptable(db)
        bps |= set((i, j) for i, j in enumerate(pt[1:], 1) if j > i)
    return bps

def mfe_intersect(seq, md, bps, fc_empty = None):
    """ Return the MFE structure given allowed base pairs.

    Makes a temporary fold compound where all bases are forbidden and only the
    specified base-pairs are allowed. If fc is specified (for efficiency
    reasons only), then it is assumed that all base-pairs are already forbidden!
    The function will return the fold compound with forbidden base-pairs again.
    """
    fc_tmp = fc_empty is None 
    if fc_empty is None:
        fc_empty = forbid_all_basepairs(seq, RNA.fold_compound(seq, md))

    for bp in bps: # enable base pairs
        fc_empty.hc_add_bp(bp[0], bp[1], RNA.CONSTRAINT_CONTEXT_ALL_LOOPS | 
                                         RNA.CONSTRAINT_CONTEXT_NO_REMOVE)
    mss, mfe = fc_empty.mfe()
    if fc_tmp:
        del fc_empty
    else:
        for bp in bps:
            fc_empty.hc_add_bp(bp[0], bp[1], RNA.CONSTRAINT_CONTEXT_NONE | 
                                             RNA.CONSTRAINT_CONTEXT_NO_REMOVE)
    return mss, int(round(mfe * 100))
 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Guide Landscape construction                                                 #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
BPD_CACHE = {}
def clear_bpd_cache():
    global BPD_CACHE
    BPD_CACHE = {}

def get_bpd_cache(s1, s2):
    if (s1, s2) not in BPD_CACHE:
        assert (s2, s1) not in BPD_CACHE
        dist = RNA.bp_distance(s1, s2)
        BPD_CACHE[(s1, s2)] = dist
        BPD_CACHE[(s2, s1)] = dist
    return BPD_CACHE[(s1, s2)]

BPD_I_CACHE = {}
def get_bpd_i_cache(p, q):
    """An alternative distance from p and q to their intersection (p and q).
    
    Note: Currently not in use.
    """
    def intersect(p, q):
        ptp = RNA.ptable(p)
        ptq = RNA.ptable(q)
        bpp = set((i, j) for i, j in enumerate(ptp[1:], 1) if j > i)
        bpq = set((i, j) for i, j in enumerate(ptq[1:], 1) if j > i)
        bpI = bpp & bpq # intersect has always less basepairs
        dpI = len(bpp) - len(bpI)
        dIq = len(bpq) - len(bpI)
        return dpI, dIq
    if (p, q) not in BPD_I_CACHE:
        assert (q, p) not in BPD_I_CACHE
        dpI, dIq = intersect(p, q)
        BPD_I_CACHE[(p, q)] = dpI
        BPD_I_CACHE[(q, p)] = dIq
    return BPD_I_CACHE[(p, q)]

def shortcut_edge_search(nodes, edges):
    """Add edges d(x,y) if E edges to n s.t. [d(x,n),d(n,y) > d(x,y)]
    """
    guide_nbrs = {n: set() for n in nodes}
    for (p, q) in edges: 
        assert p != q
        guide_nbrs[p].add(q)
    for n in nodes:
        for (p, q) in combinations(guide_nbrs[n], 2):
            assert (p, n) in edges
            assert (n, q) in edges
            if (p, q) in edges:
                continue
            dpq = get_bpd_cache(p, q)
            dpn = get_bpd_cache(p, n) 
            dnq = get_bpd_cache(n, q)
            if dpn + dnq > dpq:
                edges.add((p, q))
                edges.add((q, p))
    return edges

def guiding_edge_search(nodes, edges = None):
    """ Find all edges of the guiding neighborhood.

    Annoying case, cannot add pq at dpq = 8
     - p ..........((((.....((((.((.........)).)))).))))... 
     - q ...((((...)))).....((((.((.........)).))))........
     - i .(((......)))......((((.((.........)).))))........

    Args:
        nodes (iterable): An iterable of secondary structures.
        edges (set, optional): Provide a set of edges. They will not be deleted!

    Returns:
        (set): A set of tuples (p, q) which are neighboring structures.
    """
    if edges is None:
        edges = set()
    for p, q in combinations(nodes, 2):
        if (p, q) in edges:
            assert (q, p) in edges
            continue
        for i in nodes:
            if i == p or i == q:
                continue
            dpq = get_bpd_cache(p, q)
            dpi = get_bpd_cache(p, i)
            diq = get_bpd_cache(i, q)
            if max([dpi, diq]) < dpq:
                break
        else: # there is no i in between p and q.
            edges.add((p, q))
            edges.add((q, p))
    return edges

def guiding_node_search(seq, md, nodes, edges, fc_empty, mind = 5):
    """ For every edge in the graph, find the MFE intersect.

    Args:
        seq: The RNA sequence.
        md: ViennaRNA model details.
        nodes: An set of input structures.
        edges: A set of input reactions (node pairs).
        fc_empty (fold compound): A fold compound where
            all base-pairs are forbidden.
        mind (int, optional): A minimum base-pair distance
            to do constrained folding. Defaults to 5.

    Returns:
        (set): A set of tuples (ss, en) for new nodes that should be added.
    """
    seen = set()
    lmins = set()
    for (s1, s2) in edges:
        if get_bpd_cache(s1, s2) < mind:
            continue
        if (s2, s1) in seen:
            continue
        bps = get_basepairs([s1, s2])
        mss, mfe = mfe_intersect(seq, md, bps, fc_empty)
        if mss not in nodes:
            lmins.add((mss, mfe))
        seen.add((s1, s2))
    return lmins

def get_guide_graph(seq, md, gnodes, gedges = None, tgn_folded = None):
    """ Guide graph constuction.

    Args:
        seq: The RNA sequence (for constrained folding).
        md: ViennaRNA model details (for constrained folding).
        gnodes: A set of input guide nodes.
        gedges (optional): A set/list of input guide edges which will always be returned.
        tgn_folded: A cache to look up previous results of constrained folding for a node.

    Returns:
        (list): A list of *new* nodes.
        (set): A set of all edges.
    """
    fc_empty = forbid_all_basepairs(seq, RNA.fold_compound(seq, md))
    if tgn_folded is None:
        tgn_folded = dict()
    all_gnodes = set(gnodes) # make sure it is a set.
    new_gnodes = []
    while True:
        # Do the constrained folding for all gnodes.
        for con in list(all_gnodes):
            if not tgn_folded.get(con, False):
                css, cen = mfe_constrained(seq, md, con)
                tgn_folded[con] = True
                tgn_folded[css] = True
                if css not in all_gnodes:
                    all_gnodes.add(css)
                    new_gnodes.append((css, cen))
 
        # Find all guide edges in the graph (only force input gedges).
        all_gedges = guiding_edge_search(all_gnodes, gedges)

        # Now add all shortcut edges.
        all_gedges = shortcut_edge_search(all_gnodes, all_gedges)

        # Find new nodes that should be added to the graph.
        add_gnodes = guiding_node_search(seq, md, all_gnodes, all_gedges, fc_empty)
        if not add_gnodes:
            break
        for (ss, en) in add_gnodes:
            assert ss not in all_gnodes
            assert ss not in tgn_folded
            all_gnodes.add(ss) 
            new_gnodes.append((ss, en))
    del fc_empty
    return new_gnodes, all_gedges

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Findpath stuff                                                               #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
def common_exterior_bases(pt1, pt2):
    """ Find common unpaired bases in the exterior loop.

    Args:
        pt1: ViennaRNA pair table.
        pt2: ViennaRNA pair table.

    Yields:
        All common unpaired bases in the exterior loop.
    """
    skip = -1 
    for e in range(1, len(pt1)):
        if e > skip and pt1[e] == 0 and pt2[e] == 0: # a common unpaired base
            yield e
        elif pt1[e] > e or pt2[e] > e: # found a base-pair
            skip = max([skip, pt1[e], pt2[e]])

def common_basepairs(pt1, pt2):
    """Find common base-pairs in two pair tables.

    Args:
        pt1: ViennaRNA pair table.
        pt2: ViennaRNA pair table.

    Yields:
        Iterates over all common base-pairs.
    """
    return iter((i, j1) for i, (j1, j2) in enumerate(zip(pt1[1:], pt2[1:]), 1) 
                                            if j1 == j2 and i < j1)

def split_struct(struct, i, j, spacer = '...'):
    """Split structure into two (given vrna pairtable indices!)."""
    if j is None:
        return struct[:i], struct[i-1:]
    else:
        return struct[:i] + spacer + struct[j-1:], struct[i-1:j]

def merge_struct(outside, inside, i, j, slen = 4):
    """Merge two structures into one (given vrna pairtable indices!)."""
    if j is None:
        return outside[:i-1] + inside
    else:
        return outside[:i-1] + inside[:-1] + outside[j-len(inside)+slen:]

def findpath_merge(outside, inside, i, j):
    """ Merge two composable findpath runs.

    This is a simplified variant which just appends the two runs ...
    """
    starten = outside[0][1] + inside[0][1]
    bh1 = bh2 = 0

    path1 = []
    stepI, enI = inside[0]
    for stepO, enO in outside:
        merged = merge_struct(stepO, stepI, i, j)
        energy = enO + enI - starten
        path1.append((merged, energy))
        bh1 = max(bh1, energy)
    stepO, enO = outside[-1]
    for stepI, enI in inside[1:]:
        merged = merge_struct(stepO, stepI, i, j)
        energy = enO + enI - starten
        path1.append((merged, energy))
        bh1 = max(bh1, energy)

    path2 = []
    stepO, enO = outside[0]
    for stepI, enI in inside:
        merged = merge_struct(stepO, stepI, i, j)
        energy = enO + enI - starten
        path2.append((merged, energy))
        bh2 = max(bh2, energy)
    stepI, enI = inside[-1]
    for stepO, enO in outside[1:]:
        merged = merge_struct(stepO, stepI, i, j)
        energy = enO + enI - starten
        bh2 = max(bh2, energy)
        path2.append((merged, energy))
    return (path1, bh1) if bh1 < bh2 else (path2, bh2)

def findpath_split(seq, ss1, ss2, md, th = 5, fpwm = 4, mxb = float('inf')):
    """ Calculate findpath barriers for smaller components.

    Args:
        seq: RNA sequence.
        ss1: Structure 1.
        ss2: Structure 2.
        md: ViennaRNA model details.
        th: Threshold of how many basepairs must change for an 
            independent findpath run. Defaults to 5.
        w: Findpath width. Defaults to None.

    Returns:
        path, barrier: The folding path and the barrier height. 
            WARNING: If path splitting actually took place, then energy values
            given in the path data are only relative to the starting structure.
    """
    pt1 = RNA.ptable(ss1)
    pt2 = RNA.ptable(ss2)
    mindiff = None
    recurse = None
    for ij in chain(common_exterior_bases(pt1, pt2),
                    common_basepairs(pt1, pt2)):
        (i, j) = ij if isinstance(ij, tuple) else (ij, None)
        st1O, st1I = split_struct(ss1, i, j, spacer = '...')
        st2O, st2I = split_struct(ss2, i, j, spacer = '...')
        do = RNA.bp_distance(st1O, st2O)
        if do < th: continue
        di = RNA.bp_distance(st1I, st2I)
        if di < th: continue
        diff = abs(di-do)
        if mindiff is None or diff < mindiff:
            mindiff = diff
            seqO, seqI = split_struct(seq, i, j, spacer = 'NNN')
            recurse = ((i, j),
                       (seqO, st1O, st2O), 
                       (seqI, st1I, st2I))
        elif mindiff is not None and diff > mindiff:
            # No need to check the rest if we are getting worse.
            break

    if mindiff is not None:
        pathO, _ = findpath_split(*recurse[1], md, th, fpwm, mxb)
        pathI, _ = findpath_split(*recurse[2], md, th, fpwm, mxb)
        return findpath_merge(pathO, pathI, *recurse[0])
    else:
        w = fpwm * RNA.bp_distance(ss1, ss2)
        return call_findpath(seq, ss1, ss2, md, fpw = w, mxb = mxb)

def call_findpath(seq, ss1, ss2, md, fpw, mxb = float('inf')):
    """ Call ViennaRNA findpath.

    TODO: Find minimal motif and look it up in a cache.
    """
    fc = RNA.fold_compound(seq, md)
    if mxb == float('inf'):
        path = fc.path_findpath(ss1, ss2, width = fpw)
    else:
        e1 = round(fc.eval_structure(ss1), 2)
        dcal_bound = int(round((mxb + e1) * 100))
        path = fc.path_findpath(ss1, ss2, maxE = dcal_bound, width = fpw)
    del fc

    if len(path):
        mypath = []
        barrier = None
        for step in path:
            struct = step.s
            energy = int(round(step.en*100))
            mypath.append((struct, energy))
            if barrier is None or barrier < energy:
                barrier = energy
        barrier -= int(round(path[0].en*100))
        del step, path # potentially a good idea.
        return mypath, barrier
    return None, None

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Path flooding                                                                #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
def path_flooding(path, minh, maxlm = None):
    """ Use flooding algorithm to determine local minima on a folding path.
    
    Identifies the lowest energy of the local minimum, one representative
    stucture and "other structures" associated with that minimum. Beware that
    the "other structures" can contain saddle components when the function
    is applied to paths with degenerate saddle components.

    Args:
        path (list[(str, int)]): A list of tuples (structure, energy). 
            The energy must be an integer with unit dcal/mol.
        minh (int): Minimal height of an energy barrier separatig two
            basins [dcal/mol].
        maxlm (optional, int): The highest possible energy of a local minimum.

    Returns:
        [dict, dict]: properties of the local minima.
    """
    assert isinstance(minh, int) # dcal/mol
    assert isinstance(path[0][1], int) # dcal/mol
    indexpath = [[ss, en, step] for step, [ss, en] in enumerate(path)]
    sconf = sorted(indexpath, key = lambda x: x[1]) # by energy

    lmins = dict() # lmins[ss] = [en, set(steps)]
    ssmap = dict() # ssmap[step] = step
    for e, (ss, en, step) in enumerate(sconf):
        if step == 0 or step == len(path) - 1: # edge cases
            nbr = step + 1 if step == 0 else step - 1
            if nbr in ssmap:
                nbr = ssmap[nbr]
                lms, lme = path[nbr]
                assert en >= lme
                assert lme == lmins[lms][0]
                lmins[lms][1].add(step)
                ssmap[step] = nbr
            else:
                lmins[ss] = [en, set([step])]
                ssmap[step] = step
        else: 
            assert 0 < step < len(path)-1
            left, right = step-1, step+1
            if left in ssmap and right in ssmap: # let's connect!
                left, right = ssmap[left], ssmap[right]
                assert left != right # yay, saddle!
                lmsl, lmel = path[left]
                lmsr, lmer = path[right]
                [lower, higher] = [left, right] if lmel < lmer else [right, left]
                [lowlm, highlm] = [lmsl, lmsr] if lmel < lmer else [lmsr, lmsl]
                [lowen, highen] = [lmel, lmer] if lmel < lmer else [lmer, lmel]
                # now merge higher minimum into lower one ...
                if en - highen < minh or (maxlm and highen > maxlm):
                    lmins[lowlm][1] |= lmins[highlm][1] | set([step]) 
                    for hlm in lmins[highlm][1]:
                        ssmap[hlm] = lower
                    ssmap[step] = lower
                    del lmins[highlm]
                else: # ... or keep them separated.
                    ssmap[step] = [higher, lower]
            elif left in ssmap:
                lms, lme = path[ssmap[left]]
                assert en >= lme
                lmins[lms][1].add(step)
                ssmap[step] = ssmap[left]
            elif right in ssmap:
                lms, lme = path[ssmap[right]]
                assert en >= lme
                lmins[lms][1].add(step)
                ssmap[step] = ssmap[right]
            else:
                lmins[ss] = [en, set([step])]
                ssmap[step] = step
    return ssmap

def edge_flooding(fp, s1, s2, e1, e2, minh = None):
    """ Connect two arbitrary secondary structures.
        Args:
            minh (int): Minimal height of an energy barrier separatig 
                        two basins in dcal/mol.
    """
    if isinstance(fp, tuple):
        (seq, md, fpwm) = fp
        path, barrier = findpath_split(seq, s1, s2, md, fpwm = fpwm)
        if path[0][1] == 0:
            for i in range(len(path)):
                path[i] = (path[i][0], path[i][1] + e1)
        sen = e1 + barrier
    else:
        raise NotImplementedError('This is just a placeholder for new findpath object.')
        #sen, path = findpath_max(fp, s1, s2)

    assert path[0][1] == e1

    if minh is not None:
        # NOTE: We only allow to add minima with energy between start/stop structure.
        maxlme = max(path[0][1], path[-1][1])
        ssmap = path_flooding(path, minh, maxlm = maxlme)
        for e, si in enumerate(sorted(ssmap)):
            assert e == si # actually, no need for e, right?
            lm = ssmap[si]
            if isinstance(lm, list):
                [si1, si2] = sorted(lm)
                (ssB, enB) = path[si]
                (ss1, en1) = path[si1]
                (ss2, en2) = path[si2]
                yield ss1, en1, ssB, enB, ss2, en2
            elif e == 0 and si != lm:
                #drlog.warning(f'Starting structure {s1} is not a mimimum!')
                (ss2, en2) = path[lm]
                ssB, enB = s1, e1
                for bb in range(si+1, lm):
                    if path[bb][1] > enB:
                        (ssB, enB) = path[bb]
                yield s1, e1, ssB, enB, ss2, en2
            elif e == len(ssmap) - 1 and si != lm:
                #drlog.warning(f'Final structure {s2} is not a mimimum!')
                (ss1, en1) = path[lm]
                ssB, enB = s2, e2
                for bb in range(lm+1, si):
                    if path[bb][1] > enB:
                        (ssB, enB) = path[bb]
                yield ss1, en1, ssB, enB, s2, e2
    else:
        yield s1, e1, None, sen, s2, e2

def neighborhood_flooding(fp, ndata, gedges, tedges = None, minh = None):
    """ Calculate flooded paths for all edges.

    Note: Modifies arguments ndata and tedges.

    Starting with ndata and tedges (which may be None), 
    we add all the guide-edges to the graph.

    Args:
        ndata:
        gedges:
        tedges: 
    """
    if tedges is None:
        tedges = dict()

    # guide nbrs for each node that are not in known transition edges.
    guide_nbrs = {k: set() for k in ndata}
    for (p, q) in gedges: 
        assert p != q
        if (p, q) not in tedges:
            guide_nbrs[p].add(q)

    # all neighbors from transition edges
    tstep_nbrs = {k: set() for k in ndata}
    for (p, q) in tedges: 
        assert p != q
        assert tedges[(p, q)]['saddle_energy'] is not None
        tstep_nbrs[p].add(q)

    while not all(len(v) == 0 for v in guide_nbrs.values()):
        seen = set()
        new_gedges = []
        for s2 in sorted(ndata, key = lambda x: (ndata[x]['energy'], x), 
                         reverse = True):
            for s1 in sorted(guide_nbrs[s2], key = lambda x: (ndata[x]['energy'], x), 
                             reverse = True):
                assert s1 != s2
                if (s1, s2) in seen:
                    continue

                assert ndata[s1]['energy'] <= ndata[s2]['energy']
                e2 = ndata[s2]['energy']
                e1 = ndata[s1]['energy']
                
                for (ss2, en2, ssB, enB, ss1, en1) in edge_flooding(fp, s2, s1, 
                                                                    e2, e1, minh):
                    assert ss1 != ss2
                    if ss2 == s2 and ss1 == s1: 
                        # adding the direct connection.
                        assert ndata[ss2]['energy'] == en2
                        assert ndata[ss1]['energy'] == en1
                        tedges[(ss2, ss1)] = {'saddle_energy': enB}
                        tedges[(ss1, ss2)] = {'saddle_energy': enB}
                        tstep_nbrs[ss2].add(ss1)
                        tstep_nbrs[ss1].add(ss2)
                        guide_nbrs[ss2].remove(ss1)
                        guide_nbrs[ss1].remove(ss2)
                    elif ss2 in tstep_nbrs.get(ss1, set()): 
                        # we might have found a better transition energy.
                        assert ss1 in tstep_nbrs[ss2]
                        # TODO: is this ncecessary?
                        #assert ss2 not in guide_nbrs[ss1]
                        #assert ss1 not in guide_nbrs[ss2]
                        assert ndata[ss2]['energy'] == en2
                        assert ndata[ss1]['energy'] == en1
                        if tedges[(ss2, ss1)]['saddle_energy'] is not None:
                            enB = min(enB, tedges[(ss2, ss1)]['saddle_energy'])
                        tedges[(ss2, ss1)] = {'saddle_energy': enB}
                        tedges[(ss1, ss2)] = {'saddle_energy': enB}
                    else: 
                        # postpone evaluation for next time.
                        assert ss2 not in ndata or ndata[ss2]['energy'] == en2
                        assert ss1 not in ndata or ndata[ss1]['energy'] == en1
                        ndata[ss2] = ndata.get(ss2, {'energy': en2})
                        ndata[ss1] = ndata.get(ss1, {'energy': en1})
                        guide_nbrs[ss2] = guide_nbrs.get(ss2, set())
                        guide_nbrs[ss1] = guide_nbrs.get(ss1, set())
                        tstep_nbrs[ss2] = tstep_nbrs.get(ss2, set())
                        tstep_nbrs[ss1] = tstep_nbrs.get(ss1, set())
                        new_gedges.append((ss2, ss1))
                if s2 not in tstep_nbrs[s1]: 
                    assert s1 != s2
                    assert s1 not in tstep_nbrs[s2] 
                    guide_nbrs[s2].remove(s1)
                    guide_nbrs[s1].remove(s2)
                    tedges[(s1, s2)] = {'saddle_energy': None}
                    tedges[(s2, s1)] = {'saddle_energy': None}
                seen.add((s2, s1))
        for (p, q) in new_gedges:
            assert p != q
            guide_nbrs[p].add(q)
            guide_nbrs[q].add(p)
    return ndata, tedges

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Top-Down coarse graining                                                     #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
def top_down_coarse_graining(ndata, edata, minh = 1):
    """Top down coarse graining procedure.

    Returns:
        cg_nodes: All local minima.
        cg_edges: All edges connecting local minima.
        cg_basins: A mapping from representatives to their basins.
    """
    # minh in dcal/mol
    assert minh is not None and minh > 0
    cg_ndata = {k: v for (k, v) in ndata.items()}
    cg_edata = {k: v for (k, v) in edata.items() if v['saddle_energy'] is not None}
    cg_basin = {}

    successors = {k: set() for k in ndata}
    for (x, y) in cg_edata:
        assert (y, x) in cg_edata
        successors[x].add(y)

    for node in sorted(ndata, key = lambda x: (ndata[x]['energy'], x), reverse = True):
        en = cg_ndata[node]['energy']
        nbrs = successors[node]
        for nb in nbrs: # If there is a lower/equal-energy neighbor to merge with.
            if en < cg_ndata[nb]['energy']:
                # Assert that the higher energy structure is a delta-minimum
                assert (cg_edata[(node, nb)]['saddle_energy'] - 
                        cg_ndata[nb]['energy']) >= minh
                continue
            barrier = cg_edata[(node, nb)]['saddle_energy'] - en
            if barrier < minh:
                break
        else: # It is a local minimum!
            cg_basin[node] = cg_basin.get(node, set())
            continue

        # It is a transition structure!
        fastn = sorted([n for n in nbrs if 
                            cg_edata[(node, n)]['saddle_energy'] - en < minh], 
                            key = lambda x: (cg_ndata[x]['energy'], x))

        # starting maximum barrier is just a tick lower than minh
        seen = set()
        reps, min_se = set(), en + minh - 1
        for nbr1 in fastn:
            assert ndata[nbr1]['energy'] <= en
            se1 = cg_edata[(node, nbr1)]['saddle_energy']
            if se1 < min_se: 
                # a neighbor with so far lowest saddle energy.
                reps, min_se = set([nbr1]), se1
            elif se1 == min_se: 
                # a neighbor with equally best saddle energy.
                reps.add(nbr1)
            for nbr2 in nbrs:
                if nbr1 == nbr2 or (nbr2, nbr1) in seen:
                    continue
                se2 = cg_edata[(nbr2, node)]['saddle_energy']
                se = max(se1, se2)
                if (nbr1, nbr2) in cg_edata:
                    se = min(se, cg_edata[(nbr1, nbr2)]['saddle_energy'])
                assert isinstance(se, int)
                cg_edata[(nbr1, nbr2)] = {'saddle_energy': se}
                cg_edata[(nbr2, nbr1)] = {'saddle_energy': se}
                successors[nbr1].add(nbr2)
                successors[nbr2].add(nbr1)
                seen.add((nbr1, nbr2))

        # Update the basins given the current representatives.
        if node in cg_basin:
            basin = cg_basin[node]
            del cg_basin[node]
        else:
            basin = set()
        for rep in reps: # Assume reps are lmins (for now).
            cg_basin[rep] = cg_basin.get(rep, set())
            cg_basin[rep].add(node)
            cg_basin[rep] |= basin

        # Remove the node.
        for nb in nbrs:
            del cg_edata[(nb, node)]
            del cg_edata[(node, nb)]
            successors[nb].remove(node)
        del cg_ndata[node]
        del successors[node]
    return cg_ndata, cg_edata, cg_basin

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Main stuff                                                                   #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
def costruct(seq, cut = None):
    """ Translate between 'input' and 'internal' RNAcofold sequences """
    delim = '&'
    if delim in seq:
        cut = seq.index(delim) + 1
        table = str.maketrans(dict.fromkeys('&'))
        seq = seq.translate(table)
    elif cut and cut != -1:
        seq = seq[:cut-1] + delim + seq[cut-1:]
        cut = None
    return seq, cut

def as_barfile(seq, cgn, cge, cgm):
    out = '  ID {}  Energy  Entropy  {}\n'.format(seq,
            " ".join(map("{:7d}".format, range(1, len(cgn)+1))))
    for e, node in enumerate(sorted(cgn, key = lambda x: cgn[x]['energy']), 1):
        ne = cgn[node]['energy']
        # Calculate barrier heights to all other basins.
        barstr = ''
        for other in sorted(cgn, key = lambda x: cgn[x]['energy']):
            oe = cgn[other]['energy']
            sE = cge[(node, other)]['saddle_energy'] if (node, other) in cge else None
            if sE is not None:
                barstr += ' {:7.2f}'.format((sE - ne)/100)
            else:
                barstr += ' {:7.2f}'.format(float('nan'))
        out += "{:>4d} {}  {:>6.2f} {:>8d} {} \n".format(e, 
                node, cgn[node]['energy']/100, len(cgm[node]), barstr)
    return out

def call_findpath_exe():
    """ A wrapper for ViennaRNA findpath functions. """
    import sys
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action='count', default = 0,
            help = "Verbose output, e.g. the folding pathway. (-vv increases verbosity.)")
    parser.add_argument("-w","--search-width-multiplier", type = int, default = 4, 
            help = "Adjust upper bound for findpath search.")
    parser.add_argument("-m","--max-barrier", type = float, default = None,
            help = """Specify upper bound for barrier energy relative to 
                      first stucture (kcal/mol).""")
    parser.add_argument("--split", action = "store_true",
            help = "Split findpath into subproblems, if possible.")
    parse_model_details(parser)
    args = parser.parse_args()

    # Set model details.
    md = RNA.md()
    md.temperature = args.temp
    md.dangles = args.dangles
    md.special_hp = not args.noTetra
    md.logML = 0
    md.noGU = args.noGU
    md.noGUclosure = args.noClosingGU

    for e, line in enumerate(sys.stdin):
        if e == 0 :
            seq = line.strip()
        elif e == 1 :
            ss1, cut = costruct(line.strip())
        elif e == 2 :
            ss2, cut_ = costruct(line.strip())
    if not (cut is cut_ is None):
        raise NotImplementedError(f'No cut-point support atm {cut=} {cut_=}.')

    mxb = float('inf') if args.max_barrier is None else args.max_barrier
    if args.split:
        fpwm = args.search_width_multiplier
        path, bar = findpath_split(seq, ss1, ss2, md, th = 5, fpwm = fpwm, mxb = mxb)
    else:
        fpw = args.search_width_multiplier * RNA.bp_distance(ss1, ss2)
        path, bar = call_findpath(seq, ss1, ss2, md, fpw = fpw, mxb = mxb)

    if bar is not None:
        e1 = path[0][1]
        se = e1 + bar
        print(f"Saddle: {se/100:6.2f} kcal/mol | Barrier: {bar/100:6.2f} kcal/mol | Search width: {fpw}")
    else:
        print("No valid saddle found below energy: {:6.2f} kcal/mol | Search width parameter: {}".format(args.max_barrier, fpw))

    if args.verbose:
        print(f"    {seq}")
        for e, (ss, en) in enumerate(path):
            print(f'{e:>3d} {ss} {en/100:>6.2f}')

def top_down_coarse_graining_exe():
    import sys
    import argparse
    """ A wrapper for ViennaRNA findpath functions. """
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action='count', default = 0,
            help = """Track process using verbose output.""")
    parser.add_argument("--minh", type = float, default = None, metavar = '<float>',
            help = "Set a minimum barrier height for path flooding (kcal/mol).")
    parser.add_argument("-w","--search-width-multiplier", type = int, default = 4, 
            help = "Adjust upper bound for findpath search.")
    parser.add_argument("-e", "--elementary-moves", action = "store_true",
            help = "Elementary base-pair moves only.")
    parse_model_details(parser)
    args = parser.parse_args()

    # Set model details.
    md = RNA.md()
    md.temperature = args.temp
    md.dangles = args.dangles
    md.special_hp = not args.noTetra
    md.logML = 0
    md.noGU = args.noGU
    md.noGUclosure = args.noClosingGU
    md.noLP = args.noLP #NOTE: findpath/flooding cannot do noLP!
    fpwm = args.search_width_multiplier
    minh = None if not args.minh else int(round(args.minh*100))

    # Parsing RNAsubopt-like file.
    seq = None
    ndata = dict()
    for e, line in enumerate(sys.stdin):
        if e == 0:
            seq = line.split()[0]
            continue
        l = line.split()
        if l[0]:
            en = int(round(float(l[1])*100)) if len(l) > 1 else None
            ndata[l[0]] = {'energy': en}
            if len(l) > 2:
                ndata[l[0]]['occupancy'] = float(l[2]) 
        else:
            print(f'Troubles with input line {e+1} {line=}')

    fc = RNA.fold_compound(seq, md)
    for ss in ndata:
        if ndata[ss]['energy'] is None:
            ndata[ss]['energy'] = int(round(fc.eval_structure(ss) * 100))

    # Get the guide graph for all inputs:
    print(f'# Guide graph construction with {len(ndata)} structures.')
    gnodes, gedges = get_guide_graph(seq, md, ndata)

    if len(gnodes):
        print('# Some important structures have been included automatically:')
        for gn in gnodes:
            print(gn[0], gn[1])
            ndata[gn[0]] = {'energy': gn[1]}

    if args.elementary_moves: 
        # Do you want to have a base-pair neighborhood only?
        # You probably could have gotten that much cheaper!
        print('# Using elementary moves only.')
        edata = dict()
        for (x, y) in gedges:
            if get_bpd_cache(x, y) == 1:
                edata[(x,y)] = {'saddle_energy': max(ndata[x]['energy'], 
                                                     ndata[y]['energy'])}
            else:
                edata[(x,y)] = {'saddle_energy': None}
    else:
        print('# Neighborhood flooding ...')
        ndata, edata = neighborhood_flooding((seq, md, fpwm), ndata, gedges, minh = minh)

    for x in ndata:
        ndata[x]['occupancy'] = ndata[x].get('occupancy', 0)

    if args.minh:
        print(f'# Top down coarse graining with {args.minh=} ({len(ndata)=}, {len(edata)=}).')
        edata = {k: v for k, v in edata.items() if v['saddle_energy'] is not None}
        cgn, cge, cgm = top_down_coarse_graining(ndata, edata, minh = minh)
    else:
        print(f'# No coarse graining.')
        cgn, cge, cgm = ndata, edata, {n:[] for n in ndata}
    print('# Total hidden nodes:', sum(len(cgm[n]) for n in cgn))
    print(as_barfile(seq, cgn, cge, cgm))

    mapping = {n: set() for n in ndata}
    for lmin, hidden in cgm.items():
        assert lmin in ndata
        for hn in hidden:
            assert hn in ndata
            assert 'occupancy' in ndata[hn]
            mapping[hn].add(lmin)
            
    p0 = np.zeros(len(cgn))
    for e, x in enumerate(sorted(cgn, key = lambda x:ndata[x]['energy'])):
        p0[e] = ndata[x]['occupancy'] + sum(ndata[h]['occupancy']/len(mapping[h]) for h in cgm[x])
    print(p0, sum(p0))

