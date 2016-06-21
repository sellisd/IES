from __future__ import print_function
from ete3 import Tree, NodeStyle, SeqMotifFace, TextFace

class Vividict(dict):
    # by http://stackoverflow.com/users/541136/aaron-hall
    # from http://stackoverflow.com/questions/635483/what-is-the-best-way-to-implement-nested-dictionaries-in-python/19829714#19829714
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value


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
