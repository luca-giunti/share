from regions import CircleSkyRegion
from scipy.stats import norm
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

from gammapy.analysis import Analysis
from gammapy.maps import Map
from gammapy.estimators import ExcessMapEstimator
from gammapy.datasets import MapDataset
from skimage.filters import apply_hysteresis_threshold
from gammapy.estimators.utils import find_peaks
from gammapy.makers import FoVBackgroundMaker
from gammapy.modeling.models import FoVBackgroundModel, Models

def _pretty_plot_1d_significance(result, correlation_radius, energy_edges, exclusion):
    # 1D cash significance distributions (sqrt TS)
    from scipy.stats import norm

    plt.figure(figsize=(8, 6))

    res = result["sqrt_ts"].copy().slice_by_idx({"energy" : 0})
       
    # Remove pixels with 0 counts AFTER CORRELATION, bcs the cash significance computation fails for them
    counts = result["counts"].copy().slice_by_idx({"energy" : 0})
    lowstat = counts.data < 1
    background = result["npred"].copy().slice_by_idx({"energy" : 0})
    lowbkg = background.data < 1
    
    condition = np.logical_and.reduce((~lowstat, ~lowbkg, np.isfinite(res.data)))
    significance_all = res.data[condition]
    condition = np.logical_and(condition, exclusion.data==1).reshape(condition.shape)
    significance_off = res.data[condition]

    bins = np.arange(-100, 100, 0.1)
    significance_all, edges = np.histogram(significance_all, bins=bins)
    significance_all = significance_all / np.sum(significance_all)
    significance_off, edges = np.histogram(significance_off, bins=bins)
    significance_off = significance_off / np.sum(significance_off)

    significances = [significance_all, significance_off]
    colors = ["r", "b"]
    labels = ["All bins", "OFF bins"]
    for significance, color, label in zip(significances, colors, labels): 
        plt.plot(
            edges[:-1], 
            significance, 
            drawstyle='steps-mid',
            color=color,
            label=label,
            alpha=0.5
        )

    # Now, fit the off distribution with a Gaussian
    def gaussian(x, a, mean, sigma):
        return a * np.exp(-((x - mean)**2 / (2 * sigma**2)))

    x = bins
    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(gaussian, edges[:-1], significance_off, [10000, 0, 1])
    plt.plot(
        edges[:-1], 
        gaussian(edges[:-1], *popt), 
        lw=2, 
        color="cyan",
        label=r"$\mu$ = {:.3f}, $\sigma$ = {:.3f}".format(popt[1], popt[2]) 
    ) 

    p = norm.pdf(x, 0, 1)
    p /= np.sum(p)
    plt.plot(x, p, lw=2, color="lawngreen", label=r"$\mu$ = 0, $\sigma$ = 1")

    plt.xlabel(r"Significance (Cash) [$\sigma$]")
    plt.yscale("log")
    plt.ylim(1e-5, 1e-1)
    plt.xticks(np.arange(-5, 9, 1))
    xmin, xmax = np.min(significance_all), np.max(significance_all)
    plt.xlim(-6, 10)

    plt.title(f"correlation: {correlation_radius}, energy range: {energy_edges}")
    plt.legend()
    plt.grid(which="both", ls=":")    
    plt.show()
    plt.close()
    
def _compare_exclusions(exclusions):
    plt.figure(figsize=(8, 6))
    empty = Map.from_geom(exclusions[0].geom)
    empty.plot()
    colors = ["k", "r", "g", "b", "orange", "purple"]
    for exclusion, color in zip(exclusions, colors):
        plt.contour(exclusion.data, levels=1, linewidths=2, colors=color)
    plt.show()
    plt.close()

def _plot_significance_map(result, exclusion=None):
    plt.figure(figsize=(8, 6))
    final = (result["sqrt_ts"].slice_by_idx({"energy" : 0}))
    if exclusion is not None:
        final *= exclusion
    final.plot(add_cbar=True, cmap="coolwarm", vmin=-5, vmax=5)
    cs = plt.contour(final.data, levels=(-15, -10, -7, -5, -3, 3, 5, 7, 10, 15), colors="k")
    plt.clabel(cs, colors="k")
    
def _compare_significance_maps(result_A, result_B):
    plt.figure(figsize=(8, 6))
    (result_A["sqrt_ts"] - result_B["sqrt_ts"]).plot(add_cbar=True, cmap="coolwarm", vmin=-5, vmax=5)
    plt.show()
    plt.close()

def run_data_reduction(analysis, exclusion, model_simu):
    if exclusion is not None:
        analysis.config.datasets.background = {"method": "fov_background", "exclusion": exclusion, "parameters": {"method": "scale"}}
    analysis.get_datasets()
    
    stacked = analysis.datasets["stacked"]
    bkg_model = FoVBackgroundModel(dataset_name=stacked.name)
    stacked.models = Models([model_simu, bkg_model])
    stacked.fake(random_state=0)
    stacked.models = Models([bkg_model])
    
    return stacked
    
############
def estimate_significance_minimal(
    dataset,
    correlation_radius,
    energy_edges,
    n_iter,
    sigma,
    negative=False
):    
    
    def run(dataset, correlation_radius, sigma, negative):
        estimator_excess = ExcessMapEstimator(
            correlation_radius=correlation_radius,
            n_sigma=1, 
            n_sigma_ul=3,
            selection_optional=None,
            energy_edges=energy_edges, 
        )
        result = estimator_excess.run(dataset)
        
        sources = find_peaks(
            result["sqrt_ts"],
            threshold=sigma,
            min_distance=correlation_radius,
        )
        if negative is True:
            result["sqrt_ts"].data = -result["sqrt_ts"].data
            sources.vstack(
                find_peaks(
                        result["sqrt_ts"],
                        threshold=sigma,
                        min_distance=correlation_radius,
                )
            )

        regions = []
        for source in sources:
            skydir = SkyCoord(source["ra"], source["dec"], unit="deg", frame="icrs")
            if dataset.counts.geom.to_image().contains(skydir):
                regions.append(CircleSkyRegion(skydir, 0.2 * u.deg))
        return dataset.counts.geom.to_image().region_mask(regions=regions, inside=False), result

    first_guess_exclusion, result = run(dataset, correlation_radius, sigma, negative)
    
    for nn in range(n_iter):
        print(f"###########")
        print(nn)
        print(f"###########")

        _plot_exclusion(first_guess_exclusion, "Before")
        dataset = _run_fovbkgmaker(dataset, first_guess_exclusion)
        exclusion_enlarged, result = run(dataset, correlation_radius, sigma, negative)
        _pretty_plot_1d_significance(result, correlation_radius, energy_edges, exclusion_enlarged)
            
        first_guess_exclusion = exclusion_enlarged
        
    print(f"###########")
    print("Final result")
    print(f"###########")

    plt.figure(figsize=(8, 6))
    final = (result["sqrt_ts"].slice_by_idx({"energy" : 0}) * exclusion_enlarged)
    final.plot(add_cbar=True, cmap="coolwarm", vmin=-5, vmax=5)
    cs = plt.contour(final.data, levels=(-15, -10, -7, -5, -3, 3, 5, 7, 10, 15), colors="k")
    plt.clabel(cs, colors="k")
    plt.title(f"{correlation_radius} {energy_edges}")
    plt.show()
    plt.close()
    
def estimate_significance_luca(
    dataset, 
    first_guess_exclusion,
    correlation_radii,
    energy_edgess,
    n_iter,
    low,
    high,
    negative=False
):    
    for nn in range(n_iter):
        print(f"###########")
        print(nn)
        print(f"###########")

        results = []
        for correlation_radius, energy_edges in zip(correlation_radii, energy_edgess):
            estimator_excess = ExcessMapEstimator(
                correlation_radius=correlation_radius,
                n_sigma=1, 
                n_sigma_ul=3,
                selection_optional=None,
                energy_edges=energy_edges, 
            )
            results.append(estimator_excess.run(dataset))

        def _reproject_exclusion_mask(geom, exclusion):
            """Reproject the exclusion on the dataset geometry"""
            mask_map = Map.from_geom(geom)
            coords = geom.get_coord()
            vals = exclusion.get_by_coord(coords)
            mask_map.data += vals
            return mask_map
        first_guess_exclusion = _reproject_exclusion_mask(dataset.counts.geom, first_guess_exclusion).slice_by_idx({"energy" : 0})

        plt.figure(figsize=(8, 6))
        for correlation_radius, energy_edges, result in zip(correlation_radii, energy_edgess, results):
            (result["sqrt_ts"].slice_by_idx({"energy" : 0})*first_guess_exclusion).plot(add_cbar=True, cmap="coolwarm", vmin=-5, vmax=5)
            plt.title(f"{correlation_radius} {energy_edges}")
            plt.show()
            plt.close()

        exclusion_enlarged = first_guess_exclusion.copy()
        # Enlarge exclusion mask
        for result, radius in zip(results, correlation_radii):
            mask_map_significance = result["sqrt_ts"].copy().slice_by_idx({"energy" : 0})
            mask_map_significance.data = np.nan_to_num(mask_map_significance.data)
            mask_map = mask_map_significance.copy()
            mask_map.data = ~apply_hysteresis_threshold(mask_map_significance.data, low=low, high=high)
            if negative == True:
                mask_map.data *= ~apply_hysteresis_threshold(-mask_map_significance.data, low=low, high=high)
            exclusion_enlarged *= mask_map

        plt.figure(figsize=(8, 6))
        first_guess_exclusion.plot()
        plt.contour(exclusion_enlarged.data, levels=1, colors="r", inewidths=2)
        plt.show()
        plt.close()

        for correlation_radius, energy_edges, result in zip(correlation_radii, energy_edgess, results):
            pretty_plot_1d_significance(result, correlation_radius, energy_edges, exclusion_enlarged)
            
        first_guess_exclusion = exclusion_enlarged
        
    print(f"###########")
    print("Final result")
    print(f"###########")

    for correlation_radius, energy_edges, result in zip(correlation_radii, energy_edgess, results):
        plt.figure(figsize=(8, 6))
        final = (result["sqrt_ts"].slice_by_idx({"energy" : 0}))
        final.plot(add_cbar=True, cmap="coolwarm", vmin=-5, vmax=5)
        cs = plt.contour(final.data, levels=(-15, -10, -7, -5, -3, 3, 5, 7, 10, 15), colors="k")
        plt.clabel(cs, colors="k")
        plt.title(f"{correlation_radius} {energy_edges}")
        plt.show()
        plt.close()