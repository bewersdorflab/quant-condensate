
from PYME.IO.DataSources.BaseDataSource import BaseDataSource
from PYME.IO.MetaDataHandler import get_camera_roi_origin


class DataSource(BaseDataSource):
    """
    PYME's FlatfieldDataSource requires a flatfield, this one makes it optional
    """
    moduleName = 'FlatfieldDarkCorrectedDataSource'

    def __init__(self, parentSource, mdh, flatfield=None, dark=None):
        self.source = parentSource
        self.mdh = mdh
        self.mdh['IntensityUnits'] = 'ADU'

        x0, y0 = get_camera_roi_origin(mdh)
        x1 = x0 + mdh.getOrDefault('Camera.ROIWidth', self.source.shape[0]) + 1
        y1 = y0 + mdh.getOrDefault('Camera.ROIHeight', self.source.shape[1]) + 1

        self.flat = flatfield[x0:x1, y0:y1] if flatfield is not None else 1.
        if dark is None:
            self.dark = float(self.mdh.getEntry('Camera.ADOffset'))
        else:
            self.dark = dark[x0:x1, y0:y1].astype(float)

    def getSlice(self, ind):
        return (self.source.getSlice(ind).astype(float) - self.dark) * self.flat

    def getSliceShape(self):
        return self.source.getSliceShape()

    def getNumSlices(self):
        return self.source.getNumSlices()

    def getEvents(self):
        return self.source.getEvents()

    def release(self):
        self.source.release()