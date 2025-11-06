"""
Hardy-Weinberg Equilibrium testing
"""
import numpy as np
from scipy.stats import chi2
from typing import Dict, Tuple
import pandas as pd
from collections import Counter

from forstat.utils.logger import get_logger

logger = get_logger(__name__)


def hardy_weinberg_test(data: pd.DataFrame, locus: str) -> Dict:
    """
    Test Hardy-Weinberg Equilibrium for a locus using chi-square test

    Args:
        data: DataFrame with genetic data
        locus: Locus name

    Returns:
        Dictionary with test results
    """
    allele_col1 = f"{locus}_1"
    allele_col2 = f"{locus}_2"

    if allele_col1 not in data.columns or allele_col2 not in data.columns:
        logger.error(f"Locus {locus} not found in data")
        return {'error': 'Locus not found'}

    # Get genotypes
    genotypes = []
    for _, row in data.iterrows():
        a1 = str(row[allele_col1])
        a2 = str(row[allele_col2])

        # Skip missing data
        if a1 in ['0', '00', '000'] or a2 in ['0', '00', '000']:
            continue

        # Sort alleles to ensure consistency (e.g., 01-02 same as 02-01)
        genotype = tuple(sorted([a1, a2]))
        genotypes.append(genotype)

    if len(genotypes) < 10:
        return {
            'chi_square': None,
            'p_value': None,
            'df': None,
            'hwe_status': 'Insufficient data',
            'n_genotypes': len(genotypes)
        }

    # Count genotypes
    genotype_counts = Counter(genotypes)

    # Get allele frequencies
    alleles = []
    for g in genotypes:
        alleles.extend(g)
    allele_counts = Counter(alleles)
    total_alleles = len(alleles)
    allele_freqs = {allele: count / total_alleles for allele, count in allele_counts.items()}

    # Calculate expected genotype frequencies under HWE
    n_samples = len(genotypes)
    expected_counts = {}

    allele_list = sorted(allele_freqs.keys())
    for i, a1 in enumerate(allele_list):
        for a2 in allele_list[i:]:
            genotype = tuple(sorted([a1, a2]))
            p_a1 = allele_freqs[a1]
            p_a2 = allele_freqs[a2]

            if a1 == a2:
                # Homozygote: p^2
                expected = p_a1 ** 2 * n_samples
            else:
                # Heterozygote: 2pq
                expected = 2 * p_a1 * p_a2 * n_samples

            expected_counts[genotype] = expected

    # Chi-square test
    chi_square = 0
    for genotype, observed in genotype_counts.items():
        expected = expected_counts.get(genotype, 0)
        if expected > 0:
            chi_square += (observed - expected) ** 2 / expected

    # Degrees of freedom: n_genotypes - n_alleles
    df = len(expected_counts) - len(allele_freqs)
    if df < 1:
        df = 1

    # P-value
    p_value = 1 - chi2.cdf(chi_square, df)

    # Interpretation
    if p_value > 0.05:
        hwe_status = 'In equilibrium'
    else:
        hwe_status = 'Deviates from equilibrium'

    return {
        'chi_square': chi_square,
        'p_value': p_value,
        'df': df,
        'hwe_status': hwe_status,
        'n_genotypes': n_samples,
        'n_alleles': len(allele_freqs),
        'observed_counts': dict(genotype_counts),
        'expected_counts': expected_counts
    }


def test_all_loci_hwe(genetic_data) -> Dict:
    """
    Test Hardy-Weinberg Equilibrium for all loci

    Args:
        genetic_data: GeneticData object

    Returns:
        Dictionary with results per locus
    """
    results = {}

    for locus in genetic_data.loci:
        logger.info(f"Testing HWE for {locus}")
        results[locus] = hardy_weinberg_test(genetic_data.data, locus)

    return results
