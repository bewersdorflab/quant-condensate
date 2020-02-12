
from PYME.IO.DataSources.BaseDataSource import BaseDataSource


class DataSource(BaseDataSource):
    """
    """
    moduleName = 'MeanNormalizedDataSource'

    def __init__(self, parentSource, mdh):
        self.source = parentSource
        self.mdh = mdh

        self._first_frame_mean = self.source.getSlice(0).mean()
        self.mdh['FirstFrameMean'] = self._first_frame_mean

    def getSlice(self, ind):
        frame = self.source.getSlice(ind)
        return frame * (self._first_frame_mean / frame.mean())

    def getSliceShape(self):
        return self.source.getSliceShape()

    def getNumSlices(self):
        return self.source.getNumSlices()

    def getEvents(self):
        return self.source.getEvents()

    def release(self):
        self.source.release()