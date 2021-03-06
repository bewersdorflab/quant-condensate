import numpy as np
from PYME.IO import tabular, image, MetaDataHandler
from PYME.recipes.base import register_module, ModuleBase
from PYME.recipes.traits import Input, Output, DictStrAny, CStr, Int, Float
import logging

logger = logging.getLogger(__name__)

@register_module('MeanNormalizeToFirstFrame')
class MeanNormalizeToFirstFrame(ModuleBase):
    """
    Mean-normalize all frames to the first frame. Any (analog-digital) offset should therefore be subtracted first,
    meaning this should wrap a dark-corrected datasource.
    """
    input_image = Input('input')
    output_name = Output('mean_normalized')

    def execute(self, namespace):
        from PYME.IO.image import ImageStack
        from quant_condensate import MeanNormalizedDataSource
        image = namespace[self.input_image]

        mnd = MeanNormalizedDataSource.DataSource(image.data, image.mdh)

        im = ImageStack(mnd, titleStub=self.output_name)
        im.mdh.copyEntriesFrom(image.mdh)
        im.mdh['Parent'] = image.filename
        namespace[self.output_name] = im


@register_module('FilterSpikes')
class FilterSpikes(ModuleBase):
    """
    Using a rolling window along the time (/ z) dimension, identify spikes which are greatly above the median within
    that window and remove them by replacing the value with the median.

    Parameters
    ----------
    series: PYME.IO.image.ImageStack
    time_window_size: int
        Size of window to use in rolling-median and standard deviation calculations
    threshold_factor: float
        Multiplicative factor used to set the spike threshold, which is threshold * median-absolute deviation + median,
        all calculated within the window for an individual x, y, pixel.
    threshold_change: float
        Absolute change required in a single time-step for a spike candidate to be considered a spike

    Returns
    -------
    output: PYME.IO.image.ImageStack
        Spike-filtered copy of the input series


    Notes
    -----
    Currently only set up for single-color data
    """

    input = Input('input')

    time_window_size = Int(10)
    threshold_factor = Float(5)
    threshold_change = Float(370)

    process_frames_individually = False

    output = Output('filtered')

    def execute(self, namespace):
        from scipy.stats import median_absolute_deviation
        series = namespace[self.input]

        diff = np.diff(series.data[:,:,:,0]).squeeze()
        over_jump_threshold = np.zeros(series.data.shape[:-1], dtype=bool)
        over_jump_threshold[:, :, 1:] = diff > self.threshold_change

        output = np.copy(series.data[:,:,:,0].squeeze())  # only 1 color for now

        for ti in range(series.data.shape[2] - self.time_window_size):
            data = output[:, :, ti:ti+self.time_window_size]
            median = np.median(data, axis=2)
            spikes = np.logical_and(data > (self.threshold_factor * median_absolute_deviation(data, scale=1, axis=2) + median)[:, :, None],
                                    over_jump_threshold[:, :, ti:ti+self.time_window_size])
            spike_locs = np.nonzero(spikes)
            output[spike_locs[0], spike_locs[1], spike_locs[2] + ti] = median[spike_locs[0], spike_locs[1]]

        out = image.ImageStack(data=output)
        out.mdh = MetaDataHandler.NestedClassMDHandler()
        try:
            out.mdh.copyEntriesFrom(series.mdh)
        except AttributeError:
            pass
        out.mdh['Analysis.FilterSpikes.ThresholdFactor'] = self.threshold_factor
        out.mdh['Analysis.FilterSpikes.ThresholdChange'] = self.threshold_change
        out.mdh['Analysis.FilterSpikes.TimeWindowSize'] = self.time_window_size

        namespace[self.output] = out

@register_module('SlidingWindowMAD')
class SlidingWindowMAD(ModuleBase):
    """
    Using a rolling window along the time (/ z) dimension, calculate the median-absolute deviation (MAD)

    Parameters
    ----------
    series: PYME.IO.image.ImageStack
    time_window_size: int
        Size of window to use in rolling-median and standard deviation calculations

    Returns
    -------
    output: PYME.IO.image.ImageStack
        MAD calculated within the rolling window. Note that the window size is kept constant, so output will be a
        shorter series than the input.

    Notes
    -----
    Currently only set up for single-color data
    """

    input = Input('input')

    time_window_size = Int(10)

    process_frames_individually = False

    output = Output('MAD')

    def execute(self, namespace):
        from scipy.stats import median_absolute_deviation
        series = namespace[self.input]

        steps = range(series.data.shape[2] - self.time_window_size)
        output = np.empty((series.data.shape[0], series.data.shape[1], len(steps)),
                          dtype=series.data[:,:,0, 0].dtype)  # only 1 color for now

        for ti in steps:
            output[:,:,ti] = median_absolute_deviation(series.data[:, :, ti:ti+self.time_window_size], scale=1, axis=2)

        out = image.ImageStack(data=output)
        out.mdh = MetaDataHandler.NestedClassMDHandler()
        try:
            out.mdh.copyEntriesFrom(series.mdh)
        except AttributeError:
            pass
        out.mdh['Analysis.FilterSpikes.TimeWindowSize'] = self.time_window_size

        namespace[self.output] = out

@register_module('FlatAndDarkCorrect')
class FlatAndDarkCorrect(ModuleBase):
    input_image = Input('input')
    flatfield_filename = CStr('')
    darkmap_filename = CStr('')
    output_name = Output('corrected')

    def execute(self, namespace):
        from quant_condensate import FlatfieldDarkCorrectedDataSource
        from PYME.IO.image import ImageStack
        image = namespace[self.input_image]

        if self.flatfield_filename == '':
            flat = None
        else:
            flat = ImageStack(filename=self.flatfield_filename).data[:, :, 0].squeeze()

        if not self.darkmap_filename == '':
            dark = ImageStack(filename=self.darkmap_filename).data[:, :, 0].squeeze()
        else:
            dark = None

        ffd = FlatfieldDarkCorrectedDataSource.DataSource(image.data, image.mdh, flatfield=flat, dark=dark)

        im = ImageStack(ffd, titleStub=self.output_name)
        im.mdh.copyEntriesFrom(image.mdh)
        im.mdh['Parent'] = image.filename
        if self.darkmap_filename:
            im.mdh['FlatAndDarkCorrect.Darkmap'] = self.darkmap_filename
        if self.flatfield_filename:
            im.mdh['FlatAndDarkCorrect.Flatmap'] = self.flatfield_filename
        namespace[self.output_name] = im

@register_module('ClusteringByLabel')
class ClusteringByLabel(ModuleBase):
    """

    Parameters
    ----------
    input_name : Input
        PYME.IO.ImageStack
    mask : Input
        PYME.IO.ImageStack. Optional mask to only calculate metrics

    Returns
    -------
    output_name = Output


    Notes
    -----

    """

    input_name = Input('input')

    mask = Input('')
    excitation_start_frame = Int(10)

    output_vom = CStr('')
    output_mean_pre_excitation = CStr('')
    output_name = Output('cluster_metrics')

    def execute(self, namespace):

        series = namespace[self.input_name]

        # squeeze down from 4D
        data = series.data[:,:,:].squeeze()
        if self.mask == '':  # not the most memory efficient, but make a mask
            logger.debug('No mask provided to ClusteringByLabel, analyzing full image')
            mask = np.ones((data.shape[0], data.shape[1]), int)
        else:
            mask = namespace[self.mask].data[:,:,:].squeeze()
        # toss any negative labels, as well as the zero label (per PYME clustering schema).
        labels = sorted(list(set(np.clip(np.unique(mask), 0, None)) - {0}))
        print(labels)

        n_labels = len(labels)

        # calculate the Variance_t over Mean_t
        var = np.var(data[:,:,self.excitation_start_frame:], axis=2)
        mean = np.mean(data[:,:,self.excitation_start_frame:], axis=2)

        variance_over_mean = var / mean
        if np.isnan(variance_over_mean).any():
            logger.error('Variance over mean contains NaN, see %s' % series.filename)

        mean_pre_excitation = np.mean(data[:,:,:self.excitation_start_frame], axis=2)

        cluster_metric_mean = np.zeros(n_labels)
        mean_before_excitation = np.zeros(n_labels)

        for li in range(n_labels):
            # everything is 2D at this point
            label_mask = mask == labels[li]
            cluster_metric_mean[li] = np.mean(variance_over_mean[label_mask])
            mean_before_excitation[li] = np.mean(mean_pre_excitation[label_mask])


        res = tabular.DictSource({'variance_over_mean': cluster_metric_mean,
                                  'mean_intensity_over_first_10_frames': mean_before_excitation,
                                  'labels': np.array(labels)})
        try:
            res.mdh = series.mdh
        except AttributeError:
            res.mdh = None

        namespace[self.output_name] = res

        if self.output_vom != '':
            namespace[self.output_vom] = image.ImageStack(data=variance_over_mean, mdh=res.mdh)

        if self.output_mean_pre_excitation != '':
            namespace[self.output_mean_pre_excitation] = image.ImageStack(data=mean_pre_excitation, mdh=res.mdh)

@register_module('AddMetaData')
class AddMetaData(ModuleBase):
    """
    Hack to inject missing metadata into an image / point dataset
    Parameters
    ----------
    input_name : Input
        PYME.IO.ImageStack or PYME.IO.Tabular

    Returns
    -------
    output_name : Output


    Notes
    -----

    """

    input_name = Input('input')

    metadata_to_add = DictStrAny()

    output_name = Output('with_metadata')

    def execute(self, namespace):
        from PYME.IO.MetaDataHandler import CachingMDHandler

        inp = namespace[self.input_name]
        # md_dict = namespace[]
        mdh = CachingMDHandler(self.metadata_to_add)

        try:  # add to the existing handler if there is one
            inp.mdh.copyEntriesFrom(mdh)
        except:  # move on with our new handler
            inp.mdh = mdh
        namespace[self.output_name] = inp


@register_module('SwapColorAndSlice')    
class SwapColorAndSlice(ModuleBase):
    """swap slice (z/t) with color"""
    input_name = Input('input')
    output_name = Output('swapped')

    def execute(self, namespace):
        from quant_condensate.SwapColorAndSliceDataSource import DataSource
        from PYME.IO.MetaDataHandler import DictMDHandler
        from PYME.IO.image import ImageStack
        im = namespace[self.input_name]
        mdh = DictMDHandler()
        mdh.copyEntriesFrom(im.mdh)
        mdh['SwapColorAndSlice'] = True
        namespace[self.output_name] = ImageStack(DataSource(im.data), mdh=mdh)

@register_module('OffsetMeasure')
class OffsetMeasure(ModuleBase):
    """

    Notes
    -----
    - voxelsize must be the same for all inputs
    - all inputs must be single-channel
    - inputs can be labeled (multiple objects, e.g. label 0, label 1, etc.) but
    the masks will be thresholded as everything that's not zero.
    """
    input_chan0 = Input('chan0')
    input_mask0 = Input('mask0')
    input_chan1 = Input('chan1')
    input_mask1 = Input('mask1')
    output_name = Output('measures')

    _dtype = [
        ('offset', float), ('com0', '<3f8'), ('com1', '<3f8'),  # [nm]
        ('n_overlapping', int), ('n_0', int), ('n_1', int),
        ('n_total', int), ('fractional_volume_overlap', float),
        ('fractional_intensity_overlap', float), ('intensity_total', int),
        ('intensity0', int), ('intensity1', int)
    ]

    def execute(self, namespace):
        from scipy.ndimage import center_of_mass
        from PYME.IO.MetaDataHandler import DictMDHandler
        from PYME.IO import tabular
        
        chan0 = namespace[self.input_chan0]
        mdh = DictMDHandler()
        mdh.copyEntriesFrom(chan0.mdh)
        vx, vy, vz = chan0.voxelsize
        chan0 = np.stack([chan0.data[:,:,t,0].squeeze() for t in range(chan0.data.shape[2])], axis=2)
        mask0 = namespace[self.input_mask0]
        mask0 = np.stack([mask0.data[:,:,t,0].squeeze() for t in range(mask0.data.shape[2])], axis=2)
        mask0 = mask0 > 0
        chan1 = namespace[self.input_chan1]
        chan1 = np.stack([chan1.data[:,:,t,0].squeeze() for t in range(chan1.data.shape[2])], axis=2)
        mask1 = namespace[self.input_mask1]
        mask1 = np.stack([mask1.data[:,:,t,0].squeeze() for t in range(mask1.data.shape[2])], axis=2)
        mask1 = mask1 > 0

        com0 = center_of_mass(chan0, mask0)  # [px]
        com1 = center_of_mass(chan1, mask1)
        ox = vx * (com0[0] - com1[0])  # [nm]
        oy = vy * (com0[1] - com1[1])
        oz = vz * (com0[2] - com1[2])
        offset = np.sqrt((ox ** 2) + (oy ** 2) + (oz ** 2))

        n0 = mask0.sum()
        n1 = mask1.sum()
        n_total = n0 + n1
        mask_both = mask0 * mask1
        intensity_overlap = (mask_both * (chan0 + chan1)).sum()
        intensity0 = (chan0 * mask0).sum()
        intensity1 = (chan1 * mask1).sum()
        intensity_total = intensity0 + intensity1
        n_overlapping = np.sum(mask0 * mask1)

        
        out = np.empty((1,), dtype=self._dtype)
        out[0]['offset'] = offset
        out[0]['com0'] = com0
        out[0]['com1'] = com1
        out[0]['n_overlapping'] = n_overlapping
        out[0]['n_0'] = n0
        out[0]['n_1'] = n1
        out[0]['n_total'] = n_total
        out[0]['fractional_volume_overlap'] = n_overlapping / n_total
        out[0]['fractional_intensity_overlap'] = intensity_overlap / intensity_total
        out[0]['intensity_total'] = intensity_total
        out[0]['intensity0'] = intensity0
        out[0]['intensity1'] = intensity1

        out = tabular.RecArraySource(out)
        out.mdh = mdh
        namespace[self.output_name] = out


@register_module('StatisticsByLabel')
class StatisticsByLabel(ModuleBase):
    """
    
    NOTE: only operates on first colour channel of stack.

    Parameters
    ----------
    input_name : PYME.IO.ImageStack
    labels : PYME.IO.ImageStack
        labels had better be int, preferably positive. zero will be ignored 'unlabeled'
    Returns
    -------
    output_name = Output


    Notes
    -----

    """

    input_name = Input('input')
    input_labels = Input('labels')

    output_name = Output('label_metrics')

    def execute(self, namespace):
        from scipy import stats

        series = namespace[self.input_name]
        data = np.stack([series.data[:,:,t,0].squeeze() for t in range(series.data.shape[2])], axis=2)

        labels = namespace[self.input_labels].data
        labels = np.stack([labels[:,:,t,0].squeeze() for t in range(labels.shape[2])], axis=2)

        # drop zero label
        zero_counts = 0
        uni, n = np.unique(labels, return_counts=True)
        if np.any(uni < 0):
            raise ValueError('statistics by label does not support negative labels')
        if 0 in uni:
            zind = np.where(uni == 0)[0][0]
            zero_counts = n[zind]
            uni = np.delete(uni, zind)
            n = np.delete(n, zind)
        
        logger.debug('labels: %s' % (uni))
        
        n_labels = len(uni)
        var = np.empty(n_labels, dtype=float)
        mean = np.empty_like(var)
        median = np.empty_like(var)
        mode = np.empty_like(var)
        sum_ = np.empty_like(var)
        n_pixels = np.empty(n_labels, dtype=int)
        label = np.empty_like(n_pixels)
        
        I = np.argsort(labels.ravel())
        data = data.ravel()[I]

        start = zero_counts
        for li in range(n_labels):
            label_data = data[start:start + n[li]]

            var[li] = np.var(label_data)
            mean[li] = np.mean(label_data)
            median[li] = np.median(label_data)
            mode[li] = stats.mode(label_data, axis=None)[0][0]
            sum_[li] = label_data.sum()
            n_pixels[li] = len(label_data)
            label[li] = uni[li]

            start += n[li]

        # package up and ship-out results
        res = tabular.DictSource({
            'variance': var, 
            'mean': mean, 
            'median': median, 
            'mode': mode,
            'sum': sum_,
            'n_pixels': n,
            'label': label
        })
        try:
            res.mdh = series.mdh
        except:
            pass
            
        namespace[self.output_name] = res
