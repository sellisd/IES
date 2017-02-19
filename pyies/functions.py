from __future__ import division
from __future__ import print_function
from ete3 import NodeStyle, TextFace, Tree
from collections import Counter
import numpy as np
import re
import os.path
import matplotlib
import matplotlib.cm
from userOptions import basePath

colorFile = 'colors.hex'
colorFile9 = 'colors9.hex'

def makeNodeTable(phyldogTreeF, revBayesTreeF):
    """Create node dictioary to link two trees."""
    rbDict = revBayesTree2key(revBayesTreeF)
    phDict = nhx2key(phyldogTreeF)
    L = []
    if not sorted(rbDict.keys()) == sorted(phDict.keys()):
        print("Warning: keys not identical")
    else:
        for k in rbDict:
            L.append((rbDict[k], phDict[k]))
    return L

def revBayesTree2key(file):
    """Parse a revBayes node index tree file and create a key for each node."""
    t=Tree(file, format = 1)
    keyD = {}
    for node in t.traverse():
        k = "|".join(sorted([re.sub('\[&index=\d+\]','', n) for n in node.get_leaf_names()]))
        keyD[k] = nodeIndexFromString(node.name)
    return(keyD)

def nodeIndexFromString(nodeName):
    """Parse index number from string."""
    m = re.search('\[&index=(\d+)\]', nodeName)
    if m:
        nodeId = m.group(1)
        return(nodeId)
    else:
        print(node.name)
        quit(1)

def nhx2key(nhxtree):
    """Parse a PHYLDOG nhx file or string and create key for each node."""
    t = Tree(nhxtree)
    keyD = {}
    for node in t.traverse():
        k = "|".join(sorted([n for n in node.get_leaf_names()]))
        keyD[k] = node.ND
    return(keyD)

def speciesAbr2bn(abr, opt = "full", sep = "_"):
    """ Species abbreviation to binomial."""
    # opt can be full|dot
    binomial = ""
    d = {
        'ppr':('Paramecium', 'primaurelia'),
        'pbi':('Paramecium', 'biaurelia'),
        'pte':('Paramecium', 'tetraurelia'),
        'ppe':('Paramecium', 'pentaurelia'),
        'pse':('Paramecium', 'sexaurelia'),
        'poc':('Paramecium', 'octaurelia'),
        'ptr':('Paramecium', 'tredecaurelia'),
        'pso':('Paramecium', 'sonneborni'),
        'pca':('Paramecium', 'caudatum'),
        'tth':('Tetrahymena', 'thermophila')
    }
    (genus, specificEpithet) = d[abr]
    if opt == "full":
        binomial = genus + sep + specificEpithet
    elif opt == "dot":
        binomial = genus[0] + "." + sep + specificEpithet
    else:
        sys.exit(1)
    return binomial

def scaleCol(L):
    """ Create color scale from list."""
    m = matplotlib.cm.ScalarMappable(cmap = "coolwarm")
    R = {}
    for (i, rgba) in enumerate(m.to_rgba(L)):
        rgb = rgba[:3]
        R[L[i]] = matplotlib.colors.rgb2hex(rgb)
    return(R)

def readScafLength(f, abr, d):
    """ Read file with scaffold lengths into dictionary."""
    f.readline() #header
    for line in f:
        line = line.rstrip()
        (scaffold, length, gc) = line.split()
        d[abr + '.' + scaffold] = int(length)

def loadScafLengths(basePath, scafL):
    """ Load scaffold lengths from all species."""
    # read scaffold lenghts
    pprs = open(os.path.join(basePath, 'analysis/filtscaf/ppr.scaf'), 'r')
    pbis = open(os.path.join(basePath, 'analysis/filtscaf/pbi.scaf'), 'r')
    ptes = open(os.path.join(basePath, 'analysis/filtscaf/pte.scaf'), 'r')
    ppes = open(os.path.join(basePath, 'analysis/filtscaf/ppe.scaf'), 'r')
    pses = open(os.path.join(basePath, 'analysis/filtscaf/pse.scaf'), 'r')
    pocs = open(os.path.join(basePath, 'analysis/filtscaf/poc.scaf'), 'r')
    ptrs = open(os.path.join(basePath, 'analysis/filtscaf/ptr.scaf'), 'r')
    psos = open(os.path.join(basePath, 'analysis/filtscaf/pso.scaf'), 'r')
    pcas = open(os.path.join(basePath, 'analysis/filtscaf/pca.scaf'), 'r')

    readScafLength(pprs, "ppr", scafL)
    readScafLength(pbis, "pbi", scafL)
    readScafLength(ptes, "pte", scafL)
    readScafLength(ppes, "ppe", scafL)
    readScafLength(pses, "pse", scafL)
    readScafLength(pocs, "poc", scafL)
    readScafLength(ptrs, "ptr", scafL)
    readScafLength(psos, "pso", scafL)
    readScafLength(pcas, "pca", scafL)

    pprs.close()
    pbis.close()
    ptes.close()
    ppes.close()
    pses.close()
    pocs.close()
    ptrs.close()
    psos.close()
    pcas.close()

def numbered2name(string):
    """Extract species name from PHYLDOG numbered leaf."""
    return(re.sub(r'(.+_.+)_\d+', r'\1', string))

def phyldogSpeciesTree(phyldogTreeFile, brlenTreeFile, outgroupName):
    """Add branch lengths to PHYLDOG tree from a topologically equivalent tree."""
    b = Tree(brlenTreeFile)
    b.set_outgroup(b&outgroupName)

    brlenD = {}
    for node in b.traverse():
        leaveNames = [x.name for x in node.get_leaves()]
        leaveNames.sort()
        brlenD[tuple(leaveNames)] = node.dist

    t = Tree(phyldogTreeFile)

    for node in t.traverse():
        PHYLDOGid = ''
        if node.is_leaf():
            PHYLDOGid = (re.sub(r'.+_.+_(\d+)', r'\1', node.name))
            node.name = numbered2name(node.name)
        elif node.is_root():
            PHYLDOGid = '0'
        else:
            PHYLDOGid = str(int(node.support))
        node.add_feature("PHYLDOGid", PHYLDOGid)
        leaveNames = [numbered2name(x.name) for x in node.get_leaves()]
        leaveNames.sort()
        node.dist = brlenD[tuple(leaveNames)]

    return(t)

def parseClientOut(files, loglk, run):
    """Parse PHYLDOG output files and extract gene tree likelihoods."""
    for f in files:
        for line in f:
            line = line.rstrip()
            m = re.search("^Gene Family:.*?(\d+)\.opt total logLk:\s+(.*)\s+;.*$", line)
            if m:
                loglk[(m.group(1), run)] = m.group(2)
    return(loglk)

def up2sp(node):
    """Go upwards in a tree until a speciation node is found or return root."""
    while node.up:
        node = node.up
        if node.Ev == 'S':
            break
    return node

def isDescendant(nodeA, nodeB):
    """find if nodeB is descendant of nodeA"""
    if nodeB in nodeA.get_descendants():
        return 1
    else:
        return 0

def choseAncList(nodeL):
    """return ancestor node from list"""
    l = len(nodeL)
    rL = np.array([ ['-1']*l ]*l)
    for i in range(l):
        for j in range(l):
            rL[i][j] = isDescendant(nodeL[i], nodeL[j])
    for i in range(l):
        notSelf = [ns for ns in range(l) if ns != i]
        row = rL[i,notSelf]
        col = rL[notSelf,i]
        if all(v == '1' for v in row) and all(v == '0' for v in col):
            return nodeL[i]
    print(rL)
    print([n.S for n in nodeL])
    quit()

def choseAnc(nodeA, nodeB):
    """return ancestor node"""
    a = isDescendant(nodeA, nodeB)
    b = isDescendant(nodeB, nodeA)
    if a == 1 and b == 0:
        return nodeA
    elif a == 0 and b == 1:
        return nodeB
    else:
        quit(str(a) + '.' + str(b))

def isSpeciationNode(node):
    bool(node.Ev == 'S')

class Vividict(dict):
    # by http://stackoverflow.com/users/541136/aaron-hall
    # from http://stackoverflow.com/questions/635483/what-is-the-best-way-to-implement-nested-dictionaries-in-python/19829714#19829714
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

def placeDupl(t, spt):
    """Find per species tree branch number of duplications"""
    dup = Counter()
    for node in t.traverse():
        if node.Ev == 'D':
            #get upstream speciation
            nodeUp = up2sp(node)
            if nodeUp.Ev == 'D':
                nodeUp = -1 #the upmost node is a duplication node
            else:
                nodeUp = nodeUp.S
            nodesDown = []
            #get downstream speciation
            for spnode in node.traverse(is_leaf_fn = isSpeciationNode):
                if spnode.Ev == 'S':
                    nodesDown.append(spnode)
            # check if all downstream speciation nodes are the same
            allEqual = True
            for i in nodesDown:
                if nodesDown[0].S != i.S:
                    allEqual = False
                    break
            if allEqual:             # if all equal
                nodeDown = nodesDown[0].S
            else: # modify to deal with more than 2
                nodesL = []
                for i in nodesDown:
                   nodesL.append(spt.iter_search_nodes(S = i.S).next())
                nodesL = list(set(nodesL)) #remove duplicates, if any
                nodeDown = choseAncList(nodesL)
                nodeDown = nodeDown.S
            dup[(nodeUp,  nodeDown)] += 1
    return dup

def readPalette():
    f = open(colorFile, 'r')
    f.readline()
    cp = [i.rstrip() for i in f]
    return cp

def readPalette9():
    f = open(colorFile9, 'r')
    f.readline()
    cp = [i.rstrip() for i in f]
    return cp

def colorNodes( t, numbered ):
    """ Add styling to nodes of a PHYLDOG tree

    t: tree
    cp: list of colors (hex) corresponding to speciation events
    numbered: if true print also speciation numbers on nodes
    """
    cp = readPalette()
    for node in t.traverse():
        nstyle = NodeStyle()
        nstyle["size"] = 10
        if node.Ev == 'S':
            nstyle["fgcolor"] = cp[int(node.S) - 1]
            nstyle["shape"] = "circle"
            if numbered:
                node.add_face(TextFace(node.S), column = 0, position = "branch-right")
        elif node.Ev == 'D':
            nstyle["fgcolor"] = "black" #cblue4
            nstyle["shape"] = "square"
        else:
            quit()
        node.set_style(nstyle)
    return t

def addattr( t, bratr ):
    """Add branch attributes to tree faces from bratr.

    t: tree
    bratr: dictionary with attributes assigned to nodes
    """

    for node in t.traverse():
        if not node.is_root():
            branch = (node.up.S, node.S)
            if branch in bratr:
                if 'over' in bratr[branch]:
                    node.add_face(TextFace(bratr[branch]['over']), column = 0, position = "branch-top")
                if 'under' in bratr[branch]:
                    node.add_face(TextFace(bratr[branch]['under']), column = 0, position = "branch-bottom")

    return t
