"""
Per-population statistics
"""
import numpy as np
from typing import Dict
from forstat.utils.logger import get_logger
from forstat.data.models import GeneticData

logger = get_logger(__name__)


def calculate_per_population_stats(genetic_data: GeneticData) -> Dict:
    """
    Calculate comprehensive statistics for each population separately

    Args:
        genetic_data: GeneticData object

    Returns:
        Dictionary with per-population statistics
    """
    from forstat.analysis.population.allele_frequencies import (
        calculate_allele_frequencies, calculate_heterozygosity
    )
    from forstat.analysis.population.hardy_weinberg import hardy_weinberg_test
    from forstat.analysis.forensic.match_probability import calculate_random_match_probability

    populations = genetic_data.data['population'].unique()
    results = {}

    for pop in populations:
        pop_data = genetic_data.data[genetic_data.data['population'] == pop]

        # Create temporary genetic data object for this population
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

        pop_results = {
            'n_samples': len(pop_data),
            'loci': {}
        }

        # Calculate per-locus statistics for this population
        for locus in genetic_data.loci:
            locus_stats = {}

            # Allele frequencies
            allele_freqs = calculate_allele_frequencies(pop_data, locus)
            locus_stats['allele_frequencies'] = allele_freqs
            locus_stats['n_alleles'] = len(allele_freqs)

            # Heterozygosity
            het = calculate_heterozygosity(pop_data, locus)
            locus_stats['Ho'] = het['Ho']
            locus_stats['He'] = het['He']
            locus_stats['Fis'] = het['Fis']

            # Hardy-Weinberg
            hwe = hardy_weinberg_test(pop_data, locus)
            locus_stats['hwe_p_value'] = hwe.get('p_value')
            locus_stats['hwe_status'] = hwe.get('hwe_status')

            # Match probability for this locus in this population
            mp = calculate_random_match_probability(temp_data, locus)
            locus_stats['PM'] = mp.get('PM', 0)
            locus_stats['PD'] = mp.get('PD', 0)

            pop_results['loci'][locus] = locus_stats

        results[f'Population_{pop}'] = pop_results

    return results


def calculate_overall_match_probability_per_population(genetic_data: GeneticData) -> Dict:
    """
    Calculate combined match probability for each population

    Args:
        genetic_data: GeneticData object

    Returns:
        Dictionary with combined MP/PD per population
    """
    populations = genetic_data.data['population'].unique()
    results = {}

    for pop in populations:
        pop_data = genetic_data.data[genetic_data.data['population'] == pop]

        # Create temporary genetic data object
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

        # Calculate combined match probability
        from forstat.analysis.forensic.match_probability import calculate_combined_match_probability
        mp_results = calculate_combined_match_probability(temp_data)

        results[f'Population_{pop}'] = {
            'n_samples': len(pop_data),
            'loci_results': {k: v for k, v in mp_results.items() if not k.startswith('_')},
            'combined': mp_results.get('_combined', {})
        }

    return results
