"""
Genetic diversity metrics
"""
import numpy as np
from typing import Dict
import pandas as pd
from collections import Counter

from forstat.utils.logger import get_logger

logger = get_logger(__name__)


def calculate_diversity_indices(genetic_data, locus: str) -> Dict:
    """
    Calculate various diversity indices for a locus

    Args:
        genetic_data: GeneticData object
        locus: Locus name

    Returns:
        Dictionary with diversity metrics
    """
    from forstat.analysis.population.allele_frequencies import (
        calculate_allele_frequencies,
        calculate_heterozygosity
    )

    # Get allele frequencies
    allele_freqs = calculate_allele_frequencies(genetic_data.data, locus)

    if not allele_freqs:
        return {'error': 'No data available'}

    # Get heterozygosity
    het = calculate_heterozygosity(genetic_data.data, locus)

    # Number of alleles
    n_alleles = len(allele_freqs)

    # Shannon's diversity index
    shannon = -sum(p * np.log(p) for p in allele_freqs.values() if p > 0)

    # Effective number of alleles (inverse Simpson)
    simpson = sum(p ** 2 for p in allele_freqs.values())
    effective_n_alleles = 1 / simpson if simpson > 0 else 0

    # Allelic richness (simple count)
    allelic_richness = n_alleles

    # Polymorphic Information Content (PIC)
    pic = 0
    freq_list = list(allele_freqs.values())
    for i in range(len(freq_list)):
        for j in range(i + 1, len(freq_list)):
            pic += 2 * (freq_list[i] ** 2) * (freq_list[j] ** 2)
    pic = 1 - sum(p ** 2 for p in freq_list) - pic

    return {
        'n_alleles': n_alleles,
        'Ho': het['Ho'],
        'He': het['He'],
        'Fis': het['Fis'],
        'shannon_index': shannon,
        'effective_n_alleles': effective_n_alleles,
        'allelic_richness': allelic_richness,
        'PIC': pic,
        'allele_frequencies': allele_freqs
    }


def calculate_diversity_per_population(genetic_data) -> Dict:
    """
    Calculate diversity metrics per population

    Args:
        genetic_data: GeneticData object

    Returns:
        Dictionary with diversity per population
    """
    populations = genetic_data.data['population'].unique()
    results = {}

    for pop in populations:
        pop_data = genetic_data.data[genetic_data.data['population'] == pop]

        # Create temporary genetic data object for this population
        from forstat.data.models import GeneticData
        temp_data = GeneticData(
            title=f"{genetic_data.title} - Population {pop}",
            file_path=genetic_data.file_path,
            file_name=genetic_data.file_name,
            file_type=genetic_data.file_type,
            loci=genetic_data.loci,
            n_loci=genetic_data.n_loci,
            n_samples=len(pop_data),
            n_populations=1,
            data=pop_data,
            populations=[]
        )

        pop_results = {}
        for locus in genetic_data.loci:
            pop_results[locus] = calculate_diversity_indices(temp_data, locus)

        results[f"Population_{pop}"] = pop_results

    return results


def calculate_overall_diversity(genetic_data) -> Dict:
    """
    Calculate overall diversity statistics across all loci

    Args:
        genetic_data: GeneticData object

    Returns:
        Dictionary with overall diversity metrics
    """
    results = {}

    for locus in genetic_data.loci:
        logger.info(f"Calculating diversity for {locus}")
        results[locus] = calculate_diversity_indices(genetic_data, locus)

    # Calculate means across loci
    all_metrics = ['Ho', 'He', 'Fis', 'shannon_index', 'effective_n_alleles', 'PIC']
    averages = {}

    for metric in all_metrics:
        values = [r[metric] for r in results.values()
                 if metric in r and r[metric] is not None]
        if values:
            averages[f'mean_{metric}'] = np.mean(values)

    results['_summary'] = averages

    return results
