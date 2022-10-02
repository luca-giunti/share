import numpy as np
import astropy.units as u

from gammapy.datasets import SpectrumDataset
from gammapy.modeling import Parameter, Parameters
from gammapy.modeling.models import DatasetModels, SpectralModel
from gammapy.maps import RegionGeom, RegionNDMap, MapAxis

from sherpa.astro.instrument import RSPModelPHA

class SherpaSpectrumDataset(SpectrumDataset):
    """
    Parameters
    ----------
    data: `sherpa.data.Data` or `sherpa.data.DataSimulFit`
        An instance of a sherpa.data data-derived class
    stat: str or sherpa.stats.Stat instance
        A fit statistic from sherpa.stats
    """
    
    tag = "SherpaSpectrumDataset"

    # TODO: implement a model wrapper to fit any kind of Gammapy model (for now only sherpa models work)
    # TODO: implement read methods
    # TODO: handle the REGION and GTI information somewhere
    # TODO: possibly add counts, bkg and IRFs as properties
    # TODO: handle mask safe and mask fit somehow
    # TODO: implement compute residuals
    
    def __init__(self, data, stat, **kwargs):
        self.data = data
        if isinstance(stat, str) and not stat in ui.list_stats():
            raise ValueError("Invalid stat")
        self._stat = stat
        super().__init__(**kwargs)
            
    @property
    def _ebins(self):
        lo, hi = self.data._get_ebins()
        mask = self.data.mask
        lo = lo[mask].squeeze()
        hi = hi[mask].squeeze()
        return lo, hi
            
    @property
    def _energy_axis_reco(self):
        """Grouped reconstructed energy axis."""
        lo, hi = self._ebins
        energy_edges = np.append(lo, hi[-1])
        energy_edges = energy_edges
        # TODO: is the unit always keV and the interp always "lin"?
        return MapAxis.from_energy_edges(energy_edges, unit=u.keV, name="energy", interp='lin')

    @property
    def _geom(self):
        # TODO: read the region or pass it in init
        return RegionGeom(region=None, axes=[self._energy_axis_reco])
        
    @property
    def counts(self):
        # TODO: read the units somehow. Hardcoded for now
        data = self.data.get_y(filter=True)
        return RegionNDMap.from_geom(
            self._geom, 
            data=data, 
            unit='keV-1 s-1', 
        )
    
    @counts.setter
    def counts(self, counts):
        pass
    
    @property
    def background(self):
        # TODO: read the units somehow. Hardcoded for now
        data = self.data.get_background().get_y(filter=True)
        return RegionNDMap.from_geom(
            self._geom, 
            data=data, 
            unit='keV-1 s-1', 
        )

    @background.setter
    def background(self, background):
        pass

    def npred_signal(self):   
        data = self.data.eval_model_to_fit(self._sherpa_model)
        data /= self.data.exposure
        xlo, xhi = self._ebins
        data /= (xhi - xlo)
        # TODO: read the units somehow. Hardcoded for now
        return RegionNDMap.from_geom(
            self._geom, 
            data=data, 
            unit='keV-1 s-1', 
        )

    @property
    def _sherpa_model(self):
        """Sherpa model convolved with the IRFs"""
        data = self.data
        arf = self.data.get_arf()
        rmf = self.data.get_rmf()
        model = self.models[0]
        return RSPModelPHA(arf, rmf, data, model)
        
    @property
    def models(self):
        return self._models
        
    @models.setter
    def models(self, models):
        """Models setter. For now only sherpa-derived models (`sherpa.models.model.Model` or 
        `sherpa.models.model.SimulFitModel`) are supported"""
        if models is not None:
            models = DatasetModels(models)
            models = models.select(datasets_names=self.name)
            for model in models:
                model._exposure = self.data.exposure
        self._models = models   
               
    def notice(self, lo, hi, ignore=False):
        # TODO: implement mask safe and mask fit accordingly
        self.data.notice(lo, hi, ignore)
        
    @property
    def stat(self):
        """Fit statistic"""
        return self._stat
    
    @stat.setter
    def stat(self, stat):
        """Fit statistic setter"""
        self._stat = stat

    @property
    def stat_type(self):
        """Fit statistic name"""
        return self._stat.name

    def _calc_stat(self):
        stat_sum, stat_array = self.stat.calc_stat(self.data, self._sherpa_model)
        return stat_sum, stat_array
    
    def stat_sum(self):
        """Total statistic given the current model parameters."""
        return float(self._calc_stat()[0])

    def stat_array(self):
        """Likelihood per bin given the current model parameters"""
        return self._calc_stat()[0].astype(float)
    

class SherpaSpectralModel(SpectralModel):
    """A wrapper for Sherpa spectral models.
    Parameters
    ----------
    sherpa_model :
        An instance of the models defined in `~sherpa.models` or `~sherpa.astro.xspec`.
    default_units : tuple
        Units of the input energy array and output model evaluation (find them in the sherpa/xspec docs!)
    """

    tag = ["SherpaSpectralModel", "sherpa", "xspec"]

    def __init__(
        self, sherpa_model, name=None
    ):
        self.datasets_names = None
        self._exposure = None
        self._sherpa_model = sherpa_model
        self._name = name
        self.default_parameters = self._wrap_parameters(self.pars)
        super().__init__()
    
    @property
    def ndim(self):
        return self.sherpa_model.ndim
    
    @property
    def pars(self):
        return self.sherpa_model.pars

    @property
    def is_discrete(self):
        return self.sherpa_model.is_discrete

    @property
    def name(self):
        if self._name:
            return self._name
        else:
            return self.sherpa_model.name

    @property
    def sherpa_model(self):
        if self._exposure:
            return self._exposure * self._sherpa_model
        else:
            return self._sherpa_model
        
    def _wrap_parameters(self, pars):
        "Wraps sherpa model parameters into a gammapy Parameters object"
        parameters = []
        for par in pars:
            is_norm = par.name in ["ampl", "norm", "K"]
            parameter = Parameter(
                name=par.name, 
                value=par.val, 
                frozen=par.frozen, 
                min=par.min,
                max=par.max,
                is_norm=is_norm,
            )
            parameters.append(parameter)
        return Parameters(parameters)
        
    def calc(self, pars, x, xhi=None, *args, **kwargs):
        pars = self.parameters.value # Horrible thing to do. It is due to the fact that
                                     # internally iminuit varies the the self.parameters.value
        return self.sherpa_model.calc(pars, x, xhi, *args, **kwargs)
