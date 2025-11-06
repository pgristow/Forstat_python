"""
STR (Short Tandem Repeat) specific analysis
"""
import numpy as np
from typing import Dict, List, Tuple
import pandas as pd
from collections import Counter

from forstat.utils.logger import get_logger

logger = get_logger(__name__)


def analyze_str_markers(genetic_data) -> Dict:
    """
    Comprehensive STR marker analysis

    Args:
        genetic_data: GeneticData object

    Returns:
        Dictionary with STR statistics
    """
    results = {}

    for locus in genetic_data.loci:
        logger.info(f"Analyzing STR marker {locus}")

        # Get all alleles for this locus
        alleles = genetic_data.get_alleles(locus)

        # Convert alleles to integers (repeat numbers)
        try:
            allele_repeats = [int(a) for a in alleles if a != '0']
        except ValueError:
            # Non-numeric alleles
            allele_repeats = []

        if not allele_repeats:
            results[locus] = {'error': 'No numeric allele data'}
            continue

        # Calculate statistics
        min_repeats = min(allele_repeats)
        max_repeats = max(allele_repeats)
        range_repeats = max_repeats - min_repeats
        mean_repeats = np.mean(allele_repeats)
        std_repeats = np.std(allele_repeats)

        # Allele size distribution
        allele_counts = Counter(allele_repeats)

        # Most common allele (mode)
        mode_allele = allele_counts.most_common(1)[0][0] if allele_counts else None

        results[locus] = {
            'n_alleles': len(alleles),
            'min_repeats': min_repeats,
            'max_repeats': max_repeats,
            'range': range_repeats,
            'mean_repeats': mean_repeats,
            'std_repeats': std_repeats,
            'mode_allele': mode_allele,
            'allele_distribution': dict(allele_counts)
        }

    return results


def calculate_stutter_ratios(
    peak_heights: Dict[str, float],
    true_allele: str
) -> Dict:
    """
    Calculate stutter ratios for STR analysis

    Stutter = PCR artifact that produces peaks at n-1 or n+1 repeats

    Args:
        peak_heights: Dictionary mapping allele -> peak height
        true_allele: The true allele (not stutter)

    Returns:
        Dictionary with stutter statistics
    """
    try:
        true_repeat = int(true_allele)
    except ValueError:
        return {'error': 'Non-numeric allele'}

    true_height = peak_heights.get(true_allele, 0)

    if true_height == 0:
        return {'error': 'True allele not found'}

    # n-1 stutter (most common)
    minus_1 = str(true_repeat - 1)
    minus_1_height = peak_heights.get(minus_1, 0)
    stutter_ratio_minus_1 = minus_1_height / true_height if true_height > 0 else 0

    # n+1 stutter (less common)
    plus_1 = str(true_repeat + 1)
    plus_1_height = peak_heights.get(plus_1, 0)
    stutter_ratio_plus_1 = plus_1_height / true_height if true_height > 0 else 0

    # n-2 stutter (rare)
    minus_2 = str(true_repeat - 2)
    minus_2_height = peak_heights.get(minus_2, 0)
    stutter_ratio_minus_2 = minus_2_height / true_height if true_height > 0 else 0

    return {
        'stutter_ratio_n_minus_1': stutter_ratio_minus_1,
        'stutter_ratio_n_plus_1': stutter_ratio_plus_1,
        'stutter_ratio_n_minus_2': stutter_ratio_minus_2,
        'true_allele': true_allele,
        'true_allele_height': true_height,
        'interpretation': interpret_stutter(stutter_ratio_minus_1)
    }


def interpret_stutter(stutter_ratio: float) -> str:
    """
    Interpret stutter ratio

    Args:
        stutter_ratio: Stutter/parent peak height ratio

    Returns:
        String interpretation
    """
    if stutter_ratio < 0.05:
        return "Minimal stutter (normal)"
    elif stutter_ratio < 0.15:
        return "Expected stutter range"
    elif stutter_ratio < 0.30:
        return "Elevated stutter (consider verification)"
    else:
        return "Very high stutter (potential artifact or true allele)"


def calculate_peak_height_ratios(
    heterozygote_heights: List[Tuple[float, float]]
) -> Dict:
    """
    Calculate peak height ratios for heterozygotes

    PHR = smaller peak / larger peak
    Used to detect allele dropout, degradation, or mixture

    Args:
        heterozygote_heights: List of (height1, height2) tuples for heterozygotes

    Returns:
        Dictionary with PHR statistics
    """
    if not heterozygote_heights:
        return {'error': 'No heterozygote data'}

    phrs = []
    for h1, h2 in heterozygote_heights:
        if h1 > 0 and h2 > 0:
            phr = min(h1, h2) / max(h1, h2)
            phrs.append(phr)

    if not phrs:
        return {'error': 'No valid peak height data'}

    mean_phr = np.mean(phrs)
    median_phr = np.median(phrs)
    min_phr = np.min(phrs)
    max_phr = np.max(phrs)

    # Quality assessment
    if mean_phr < 0.6:
        quality = "Poor - possible allele dropout or degradation"
    elif mean_phr < 0.7:
        quality = "Fair - acceptable but monitor"
    else:
        quality = "Good - balanced heterozygotes"

    return {
        'mean_phr': mean_phr,
        'median_phr': median_phr,
        'min_phr': min_phr,
        'max_phr': max_phr,
        'n_heterozygotes': len(phrs),
        'quality_assessment': quality
    }


def detect_null_alleles(genetic_data, locus: str, threshold: float = 0.05) -> Dict:
    """
    Detect potential null alleles (non-amplifying alleles)

    Null alleles indicated by:
    - Excess homozygosity
    - Deficiency of heterozygotes

    Args:
        genetic_data: GeneticData object
        locus: Locus name
        threshold: Threshold for null allele frequency

    Returns:
        Dictionary with null allele analysis
    """
    from forstat.analysis.population.allele_frequencies import calculate_heterozygosity

    het = calculate_heterozygosity(genetic_data.data, locus)

    Ho = het['Ho']
    He = het['He']
    Fis = het['Fis']

    # Estimate null allele frequency using Chakraborty method
    # r = (He - Ho) / (1 + He) where r is null allele frequency
    if He > 0:
        null_freq_estimate = (He - Ho) / (1 + He)
    else:
        null_freq_estimate = 0

    # Interpret
    if null_freq_estimate < threshold:
        interpretation = "No evidence of null alleles"
    elif null_freq_estimate < 0.10:
        interpretation = "Weak evidence of null alleles"
    elif null_freq_estimate < 0.20:
        interpretation = "Moderate evidence of null alleles"
    else:
        interpretation = "Strong evidence of null alleles"

    return {
        'null_allele_frequency': max(0, null_freq_estimate),  # Can't be negative
        'Ho': Ho,
        'He': He,
        'Fis': Fis,
        'interpretation': interpretation,
        'recommendation': 'Consider alternative markers' if null_freq_estimate > 0.10 else 'Marker appears reliable'
    }
