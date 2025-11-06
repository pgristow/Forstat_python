"""
Fixation Index (Fst) calculations for population differentiation
"""
import numpy as np
from typing import Dict, List
import pandas as pd
from collections import Counter

from forstat.utils.logger import get_logger

logger = get_logger(__name__)


def calculate_fst_weir_cockerham(genetic_data, locus: str) -> Dict:
    """
    Calculate Fst using Weir & Cockerham (1984) method

    Args:
        genetic_data: GeneticData object
        locus: Locus name

    Returns:
        Dictionary with Fst statistics
    """
    allele_col1 = f"{locus}_1"
    allele_col2 = f"{locus}_2"

    data = genetic_data.data

    # Get populations
    populations = data['population'].unique()
    n_pops = len(populations)

    if n_pops < 2:
        return {
            'fst': None,
            'error': 'Need at least 2 populations'
        }

    # Calculate allele frequencies per population
    pop_allele_freqs = {}
    pop_sample_sizes = {}

    for pop in populations:
        pop_data = data[data['population'] == pop]

        # Count alleles
        alleles = []
        for _, row in pop_data.iterrows():
            a1 = str(row[allele_col1])
            a2 = str(row[allele_col2])

            if a1 not in ['0', '00', '000']:
                alleles.append(a1)
            if a2 not in ['0', '00', '000']:
                alleles.append(a2)

        if not alleles:
            continue

        allele_counts = Counter(alleles)
        total = len(alleles)

        pop_allele_freqs[pop] = {allele: count / total for allele, count in allele_counts.items()}
        pop_sample_sizes[pop] = len(pop_data)

    if len(pop_allele_freqs) < 2:
        return {
            'fst': None,
            'error': 'Insufficient data'
        }

    # Get all unique alleles
    all_alleles = set()
    for freqs in pop_allele_freqs.values():
        all_alleles.update(freqs.keys())

    # Calculate overall allele frequencies
    total_samples = sum(pop_sample_sizes.values())
    overall_freqs = {}

    for allele in all_alleles:
        weighted_sum = 0
        for pop, freqs in pop_allele_freqs.items():
            freq = freqs.get(allele, 0)
            weight = pop_sample_sizes[pop]
            weighted_sum += freq * weight
        overall_freqs[allele] = weighted_sum / total_samples

    # Calculate Hs (within-population heterozygosity)
    Hs = 0
    for pop, freqs in pop_allele_freqs.items():
        pop_het = 1 - sum(freqs.get(a, 0) ** 2 for a in all_alleles)
        weight = pop_sample_sizes[pop] / total_samples
        Hs += pop_het * weight

    # Calculate Ht (total heterozygosity)
    Ht = 1 - sum(overall_freqs[a] ** 2 for a in all_alleles)

    # Calculate Fst
    if Ht > 0:
        fst = (Ht - Hs) / Ht
    else:
        fst = 0

    return {
        'fst': fst,
        'Ht': Ht,
        'Hs': Hs,
        'n_populations': n_pops,
        'n_alleles': len(all_alleles)
    }


def calculate_pairwise_fst(genetic_data, locus: str) -> Dict:
    """
    Calculate pairwise Fst between all population pairs

    Args:
        genetic_data: GeneticData object
        locus: Locus name

    Returns:
        Dictionary with pairwise Fst values
    """
    populations = genetic_data.data['population'].unique()
    n_pops = len(populations)

    if n_pops < 2:
        return {'error': 'Need at least 2 populations'}

    pairwise_fst = {}

    for i in range(n_pops):
        for j in range(i + 1, n_pops):
            pop1 = populations[i]
            pop2 = populations[j]

            # Create subset with only these two populations
            subset_data = genetic_data.data[genetic_data.data['population'].isin([pop1, pop2])]

            # Create temporary genetic data object
            from forstat.data.models import GeneticData
            temp_data = GeneticData(
                title=genetic_data.title,
                file_path=genetic_data.file_path,
                file_name=genetic_data.file_name,
                file_type=genetic_data.file_type,
                loci=genetic_data.loci,
                n_loci=genetic_data.n_loci,
                n_samples=len(subset_data),
                n_populations=2,
                data=subset_data,
                populations=[]
            )

            result = calculate_fst_weir_cockerham(temp_data, locus)
            pairwise_fst[f"Pop{pop1}_vs_Pop{pop2}"] = result.get('fst', None)

    return pairwise_fst


def calculate_fst_all_loci(genetic_data) -> Dict:
    """
    Calculate Fst for all loci

    Args:
        genetic_data: GeneticData object

    Returns:
        Dictionary with Fst per locus
    """
    results = {}

    for locus in genetic_data.loci:
        logger.info(f"Calculating Fst for {locus}")
        fst_result = calculate_fst_weir_cockerham(genetic_data, locus)

        if 'error' not in fst_result:
            pairwise = calculate_pairwise_fst(genetic_data, locus)
            fst_result['pairwise_fst'] = pairwise

        results[locus] = fst_result

    # Calculate average Fst across loci
    fst_values = [r['fst'] for r in results.values() if 'fst' in r and r['fst'] is not None]
    if fst_values:
        results['_average_fst'] = np.mean(fst_values)

    return results
