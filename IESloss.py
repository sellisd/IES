#!/usr/bin/python
#figure IES loss
from __future__ import print_function
from ete3 import Tree, NodeStyle, TreeStyle
from pyies.userOptions import basePath
import os.path
import string

geneFamily = '4137'
analysis = '2'
phyldogResultsPath = "analysis/phyldogT" + analysis + "/results"
treeF = os.path.join(basePath, 'analysis', 'phyldogT' + analysis, 'results', geneFamily + ".ReconciledTree"
)
alnF = os.path.join(basePath, 'analysis', 'msas', 'filtered', 'cluster.' + geneFamily + '.nucl.fa')
t = Tree(treeF)
ts = TreeStyle()
ts.scale = 500
nstyle = NodeStyle()
nstyle["size"] = 0
for n in t.traverse():
   n.set_style(nstyle)

#seqs = SeqGroup(sequences = alnF, format = "fasta")

# for leaf in t:
#     geneId = leaf.name
#     seq = seqs.get_seq(geneId)
#     seq = seq.translate(None, string.ascii_lowercase) # keep only CDS
#     #iesmotif = [[1, len(seq), "line", 2, 5, None, None, None]]
#     #seqFace = SeqMotifFace(seq = seq, gap_format = "blank", seq_format = "line")
#     seqFace = SequenceFace(seq, seqtype = 'dna')
#     leaf.add_face(seqFace, 0, "aligned")
#
#t.show(tree_style = ts)
t.render("/Users/dsellis/projects/IES/src/figureData/tree." + geneFamily + ".svg", tree_style = ts)
