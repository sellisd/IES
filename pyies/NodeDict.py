import pandas as pd
class NodeDict:
    """Holds data for node naming schemes.
    a2b conversion functions return pd.Series
    """

    def __init__(self, nodeDictionaryFile):
        """Constructor
        Args:
        nodeDictionaryFile string  Path and filename for node dictionary file
        returns Nothing
        """
        self._df = pd.read_csv(nodeDictionaryFile, sep = "\t", dtype = 'str')

    def rb2phyldog(self, geneFamily, nodeRb):
        """translate node notation
        Args:
        geneFamily string Gene family unique Id
        nodeRb     string Node Id in revBayes numbering
        returns    string revBayes NodeId or None if node not in database
        """
        nodeId = self._df.phyldog[(self._df.cluster == geneFamily) & (self._df.rb == nodeRb)]
        if nodeId.empty:
            return None
        else:
            return(nodeId.item())

    def phyldog2rb(self, geneFamily, nodeP):
        """translate node notation
        Args:
        geneFamily string Gene family unique Id
        nodeP      string Node Id in PHYLDOG numbering
        returns    string phyldog NodeId or None if node not in database
        """
        nodeId = self._df.rb[(self._df.cluster == geneFamily) & (self._df.phyldog == nodeP)]
        if nodeId.empty:
            return None
        else:
            return(nodeId.item())
