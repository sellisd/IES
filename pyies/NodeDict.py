import pandas as pd
class NodeDict:
    """Holds data for node naming schemes.
    a2b conversion functions return pd.Series
    """

    def __init__(self, file):
        self._df = pd.read_csv(file, sep = "\t")

    def rb2phyldog(self, geneFamily, nodeRb):
        return(self._df.rb[(self._df.phyldog.isin([2])) & (self._df.geneFamily == geneFamily)])

    def phyldog2pb(self, geneFamily, nodeP):
        return(self._df.phyldog[(self._df.rb.isin([2])) & (self._df.geneFamily == geneFamily)])
