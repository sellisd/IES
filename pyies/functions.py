from __future__ import print_function
from ete3 import Tree, NodeStyle, SeqMotifFace, TextFace
from collections import defaultdict

def isDescendant(nodeA, nodeB):
    """find if nodeB is descendant of nodeA"""
    if nodeB in nodeA.get_descendants():
        return 1
    else:
        return 0

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

class Vividict(dict):
    # by http://stackoverflow.com/users/541136/aaron-hall
    # from http://stackoverflow.com/questions/635483/what-is-the-best-way-to-implement-nested-dictionaries-in-python/19829714#19829714
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

def placeDupl(t, spt):
    """Find per species tree branch number of duplications"""
    dup = defaultdict(int)
    for node in t.traverse():
        if node.Ev == 'D':
            if len(node.children) == 2:
                if node.children[0].Ev == 'S' and node.children[1].Ev == 'S' and node.children[0].S == node.children[1].S:
                    # both speciation nodes of the same speciation event
                    if node.is_root():
                        print("root")
                    else:
                        if node.up.Ev == 'S':
                            dup[(node.up.S, node.children[0].S)] += 1
                        else:
                            print(node.S + 'parent node not speciation')
                            return
                else:
                    # which is the most ancestral speciation node in the species tree?
                    anc = choseAnc(spt.iter_search_nodes(S = node.children[0].S).next(), spt.iter_search_nodes(S = node.children[1].S).next())
                    dup[(node.up.S, anc.S)] +=1
            else:
                print(node.S + 'more than two children')
                return
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
    """ Add branch attributes to tree faces from bratr
    
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
