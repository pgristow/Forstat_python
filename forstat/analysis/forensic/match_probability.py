"""
Match Probability and Power of Discrimination calculations
"""
import numpy as np
from typing import Dict, Tuple
import pandas as pd
from collections import Counter

from forstat.utils.logger import get_logger

logger = get_logger(__name__)


def calculate_random_match_probability(genetic_data, locus: str) -> Dict:
    """
    Calculate random match probability for a locus

    Args:
        genetic_data: GeneticData object
        locus: Locus name

    Returns:
        Dictionary with match probability statistics
    """
    from forstat.analysis.population.allele_frequencies import calculate_allele_frequencies

    # Get allele frequencies
    allele_freqs = calculate_allele_frequencies(genetic_data.data, locus)

    if not allele_freqs:
        return {'error': 'No data available'}

    # Calculate genotype frequencies
    alleles = sorted(allele_freqs.keys())
    genotype_probs = {}

    for i, a1 in enumerate(alleles):
        for j in range(i, len(alleles)):
            a2 = alleles[j]
            p1 = allele_freqs[a1]
            p2 = allele_freqs[a2]

            if a1 == a2:
                # Homozygote: p^2
                prob = p1 ** 2
            else:
                # Heterozygote: 2pq
                prob = 2 * p1 * p2

            genotype_probs[f"{a1}/{a2}"] = prob

    # Random match probability (sum of squared genotype frequencies)
    PM = sum(p ** 2 for p in genotype_probs.values())

    # Power of Discrimination
    PD = 1 - PM

    # Power of Exclusion (PE) - probability of excluding a random individual
    # PE = H^2 where H is heterozygosity
    from forstat.analysis.population.allele_frequencies import calculate_heterozygosity
    het = calculate_heterozygosity(genetic_data.data, locus)
    PE = het['Ho'] ** 2

    # Typical Paternity Index (assumes mother known, tests alleged father)
    # Average PI = 1 / (2 * sum(pi^2))
    sum_p_squared = sum(p ** 2 for p in allele_freqs.values())
    typical_PI = 1 / (2 * sum_p_squared) if sum_p_squared > 0 else 0

    return {
        'PM': PM,  # Match probability
        'PD': PD,  # Power of discrimination
        'PE': PE,  # Power of exclusion
        'typical_PI': typical_PI,
        'n_alleles': len(allele_freqs),
        'genotype_probabilities': genotype_probs
    }


def calculate_combined_match_probability(genetic_data) -> Dict:
    """
    Calculate combined match probability across all loci

    Args:
        genetic_data: GeneticData object

    Returns:
        Dictionary with combined statistics
    """
    results = {}
    PM_values = []
    PD_values = []
    PE_values = []
    PI_values = []

    for locus in genetic_data.loci:
        logger.info(f"Calculating match probability for {locus}")
        locus_result = calculate_random_match_probability(genetic_data, locus)

        if 'error' not in locus_result:
            results[locus] = locus_result
            PM_values.append(locus_result['PM'])
            PD_values.append(locus_result['PD'])
            PE_values.append(locus_result['PE'])
            PI_values.append(locus_result['typical_PI'])

    # Combined match probability (product across loci)
    if PM_values:
        combined_PM = np.prod(PM_values)
        combined_PD = 1 - combined_PM

        # Combined PE (product rule)
        combined_PE = np.prod(PE_values)

        # Combined PI (product rule)
        combined_PI = np.prod(PI_values)

        # Calculate how rare (1 in X people)
        if combined_PM > 0:
            one_in_X = 1 / combined_PM
        else:
            one_in_X = float('inf')

        results['_combined'] = {
            'combined_PM': combined_PM,
            'combined_PD': combined_PD,
            'combined_PE': combined_PE,
            'combined_PI': combined_PI,
            'one_in_X': one_in_X,
            'n_loci': len(PM_values)
        }

    return results


def calculate_profile_probability(genetic_data, profile: Dict[str, Tuple[str, str]]) -> Dict:
    """
    Calculate the probability of a specific genotype profile

    Args:
        genetic_data: GeneticData object
        profile: Dictionary mapping locus -> (allele1, allele2)

    Returns:
        Dictionary with profile probability
    """
    from forstat.analysis.population.allele_frequencies import calculate_allele_frequencies

    profile_prob = 1.0
    locus_probs = {}

    for locus, (a1, a2) in profile.items():
        if locus not in genetic_data.loci:
            continue

        allele_freqs = calculate_allele_frequencies(genetic_data.data, locus)

        p1 = allele_freqs.get(a1, 0)
        p2 = allele_freqs.get(a2, 0)

        if a1 == a2:
            # Homozygote
            genotype_prob = p1 ** 2
        else:
            # Heterozygote
            genotype_prob = 2 * p1 * p2

        locus_probs[locus] = genotype_prob
        profile_prob *= genotype_prob

    if profile_prob > 0:
        one_in_X = 1 / profile_prob
    else:
        one_in_X = float('inf')

    return {
        'profile_probability': profile_prob,
        'one_in_X': one_in_X,
        'locus_probabilities': locus_probs
    }
