from __future__ import print_function
from ete3 import Tree, NodeStyle, TextFace
from collections import Counter
import numpy as np
import re

def parseClientOut(files):
    """Parse PHYLDOG output files and extract gene tree likelihoods
    """
    loglk = {}
    for f in files:
        for line in f:
            line = line.rstrip()
            m = re.search("^Gene Family:.*?(\d+)\.opt total logLk:\s+(.*)\s+;.*$", line)
            if m:
                loglk[m.group(1)] = m.group(2)
    return(loglk)
    
def up2sp(node):
    """
    Go upwards in a tree until a speciation node is found or return root
    """
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
    f = open('/home/dsellis/projects/IES/src/colors.hex', 'r')
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
