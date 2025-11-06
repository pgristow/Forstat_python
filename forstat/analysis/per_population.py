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
    Calculate combined match probability for each population with 1-in-X per locus

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

        # Add 1-in-X for each locus
        loci_with_one_in_x = {}
        for locus, locus_data in mp_results.items():
            if not locus.startswith('_') and isinstance(locus_data, dict):
                pm = locus_data.get('PM', 0)
                if pm > 0:
                    locus_data['one_in_X'] = 1 / pm
                else:
                    locus_data['one_in_X'] = float('inf')
                loci_with_one_in_x[locus] = locus_data

        results[f'Population_{pop}'] = {
            'n_samples': len(pop_data),
            'loci_results': loci_with_one_in_x,
            'combined': mp_results.get('_combined', {})
        }

    # Add overall (across all populations)
    results['_overall'] = calculate_overall_summary_mp(genetic_data)

    return results


def calculate_overall_summary_mp(genetic_data: GeneticData) -> Dict:
    """Calculate overall match probability across all populations"""
    from forstat.analysis.forensic.match_probability import calculate_combined_match_probability

    overall = calculate_combined_match_probability(genetic_data)

    # Add 1-in-X for each locus
    for locus, locus_data in overall.items():
        if not locus.startswith('_') and isinstance(locus_data, dict):
            pm = locus_data.get('PM', 0)
            if pm > 0:
                locus_data['one_in_X'] = 1 / pm
            else:
                locus_data['one_in_X'] = float('inf')

    return overall


def add_overall_summary(results: Dict, genetic_data: GeneticData) -> Dict:
    """
    Add overall summary statistics across all populations

    Args:
        results: Per-population results
        genetic_data: GeneticData object

    Returns:
        Results with _overall key added
    """
    populations = [k for k in results.keys() if k.startswith('Population_')]

    if not populations:
        return results

    # Get all loci from first population
    first_pop = populations[0]
    loci = list(results[first_pop].get('loci', {}).keys())

    overall_loci = {}

    # Calculate mean statistics across populations for each locus
    for locus in loci:
        ho_values = []
        he_values = []
        fis_values = []
        n_alleles_values = []
        pd_values = []

        for pop in populations:
            locus_data = results[pop].get('loci', {}).get(locus, {})
            if locus_data:
                ho_values.append(locus_data.get('Ho', 0))
                he_values.append(locus_data.get('He', 0))
                fis_values.append(locus_data.get('Fis', 0))
                n_alleles_values.append(locus_data.get('n_alleles', 0))
                pd_values.append(locus_data.get('PD', 0))

        if ho_values:
            overall_loci[locus] = {
                'Ho_mean': np.mean(ho_values),
                'He_mean': np.mean(he_values),
                'Fis_mean': np.mean(fis_values),
                'n_alleles_mean': np.mean(n_alleles_values),
                'PD_mean': np.mean(pd_values),
                'n_populations': len(populations)
            }

    results['_overall'] = {
        'loci': overall_loci,
        'n_populations': len(populations),
        'description': 'Mean statistics across all populations'
    }

    return results
