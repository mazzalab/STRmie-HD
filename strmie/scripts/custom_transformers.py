import numpy as np
from scipy import stats
import torch


class Log1p(object):

    def __init__(self, base=np.pi):
        self.base = base
    
    def __call__(self, sample):
        series, labels = sample
        series = np.log(1+series)/np.log(self.base)
        return series, labels



class RandomNoise(object):

    def __init__(self, noise=.1, vmin=None, vmax=None):
        self.noise = noise
        self.vmin = vmin
        self.vmax = vmax

    def __call__(self, sample):
        series, labels = sample
        scale = series.std(1)*self.noise
        series = series + np.random.normal([0,0], scale, size=(series.shape[-1], 2)).T
        if (not self.vmin is None) or (not self.vmax is None):
            series = series.clip(self.vmin, self.vmax)
        return series, labels



class RandomPower(object):

    def __init__(self, alpha=10):
        """
        Transform the series through a power function.
        The exponent is randomly extracted from a random gamma distribution.
        Only alpha parameter is required, while beta is setted in order to have: alpha*beta=1.
        This in order to force the distribution to have expected value 1.
        the smaller alpha the higher the variance of the distribution.
        The recommended value for alpha is 10
        """
        self.alpha = alpha
        self.beta = 1/self.alpha

    def __call__(self, sample):
        series, labels = sample
        exp = np.random.gamma(shape=self.alpha, scale=self.beta)
        series = series**exp
        return series, labels
    


class RandomExpansion(object):

    def __init__(self, scale=1):
        self.scale = scale

    def __call__(self, sample):
        series, labels = sample
        slope = np.random.normal(0, self.scale)
        if slope>0:
            series = series * np.linspace(0,1,series.shape[-1]) * slope
        if slope<0:
            series = series * (np.linspace(0,1,series.shape[-1]) * slope -slope)
        return series, labels



class RandomShift(object):

    def __init__(self, margin: tuple):
        self.margin = margin
    
    def __call__(self, sample):
        series, labels = sample
        shift = np.random.randint(self.margin[0], self.margin[1])
        if shift>0:
            start, end = series[:,:2]
            slope = (end - start)/shift
            interpolation = (start + torch.arange(shift+1)[:,None]*slope).T
            series = torch.roll(series, shift, dims=1)
            series[:,:shift+1] = interpolation
        if shift<0:
            start, end = series[:,-2:]
            slope = (end - start)/-shift
            interpolation = (start + torch.arange(-shift+1)[:,None]*slope).T
            series = torch.roll(series, shift, dims=1)
            series[:,shift-1:] = interpolation
        labels = labels + shift
        return series, labels



class ClipLabels(object):

    def __init__(self, lower=None, upper=None):
        self.lower = lower
        self.upper = upper
    
    def __call__(self, sample):
        series, labels = sample
        labels = torch.clip(labels, self.lower, self.upper)
        return series, labels


class MinMaxScaling(object):

    def __init__(self, lower=0, upper=1):
        self.lower = lower
        self.upper = upper
    
    def __call__(self, sample):
        series, labels = sample
        series = (series.T - series.min(axis=1)[0]) / (series.max(axis=1)[0] - series.min(axis=1)[0])
        series = series * (self.upper - self.lower) + self.lower
        return series.T, labels



class RandomStretch(object):

    def __init__(self, alpha=1):
        """
        Transform the series through stretching.
        The stretching coefficient is randomly extracted from a gamma distribution.
        Only alpha parameter is required, while beta is setted in order to have: alpha*beta=1.
        This in order to force the distribution to have expected value 1.
        the smaller alpha the higher the variance of the distribution.
        """
        self.alpha = alpha
        self.beta = 1/self.alpha
    
    def __call__(self, sample):
        series, labels = sample
        cast = series.dtype
        coef = np.random.gamma(shape=self.alpha, scale=self.beta)
        x0 = np.arange(series.shape[-1])
        xf = x0 * coef
        series = torch.from_numpy(
            np.array([np.interp(x0, xf, s) for s in series]),
        )
        series.type(cast)
        labels = labels * coef
        return series, labels


class PeaksPDF(object):

    def __init__(self, tolerance=1):
        self.tolerance = tolerance
    
    def __call__(self, sample):
        series, labels = sample
        dist = stats.norm(loc=labels, scale=self.tolerance)
        labels = torch.from_numpy(
            dist.pdf(np.arange(series.shape[1])[:,None]+1)
        ).T
        labels = labels/labels.sum(1)[:,None]
        return series, labels