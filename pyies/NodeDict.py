import pandas as pd
class NodeDict:
    """Holds data for node naming schemes."""

    def __init__(file):
        self._df = pd.readsvn(file, sep = "\t")
