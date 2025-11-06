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


def calculate_per_population_fst(genetic_data: GeneticData) -> Dict:
    """
    Calculate Fst-related statistics for each population
    Shows Hs (within-population heterozygosity) and pairwise Fst values

    Args:
        genetic_data: GeneticData object

    Returns:
        Dictionary with per-population Fst statistics
    """
    from forstat.analysis.population.fst import calculate_fst_weir_cockerham
    from forstat.analysis.population.allele_frequencies import calculate_allele_frequencies
    from collections import Counter

    populations = genetic_data.data['population'].unique()
    results = {}

    for pop in populations:
        pop_data = genetic_data.data[genetic_data.data['population'] == pop]
        pop_results = {
            'n_samples': len(pop_data),
            'loci': {}
        }

        # Calculate per-locus statistics for this population
        for locus in genetic_data.loci:
            allele_col1 = f"{locus}_1"
            allele_col2 = f"{locus}_2"

            # Calculate Hs (within-population heterozygosity) for this population
            alleles = []
            for _, row in pop_data.iterrows():
                a1 = str(row[allele_col1])
                a2 = str(row[allele_col2])
                if a1 not in ['0', '00', '000']:
                    alleles.append(a1)
                if a2 not in ['0', '00', '000']:
                    alleles.append(a2)

            if alleles:
                allele_counts = Counter(alleles)
                total = len(alleles)
                allele_freqs = {allele: count / total for allele, count in allele_counts.items()}

                # Hs = 1 - sum(p^2) for this population
                Hs = 1 - sum(freq ** 2 for freq in allele_freqs.values())

                # Fis for this population (inbreeding coefficient)
                # Count heterozygotes
                n_het = 0
                n_total = 0
                for _, row in pop_data.iterrows():
                    a1 = str(row[allele_col1])
                    a2 = str(row[allele_col2])
                    if a1 not in ['0', '00', '000'] and a2 not in ['0', '00', '000']:
                        n_total += 1
                        if a1 != a2:
                            n_het += 1

                Ho = n_het / n_total if n_total > 0 else 0
                He = Hs
                Fis = (He - Ho) / He if He > 0 else 0

                pop_results['loci'][locus] = {
                    'Hs': Hs,
                    'Ho': Ho,
                    'He': He,
                    'Fis': Fis,
                    'n_alleles': len(allele_freqs)
                }
            else:
                pop_results['loci'][locus] = {
                    'Hs': 0,
                    'Ho': 0,
                    'He': 0,
                    'Fis': 0,
                    'n_alleles': 0
                }

        results[f'Population_{pop}'] = pop_results

    # Add overall Fst calculations
    from forstat.analysis.population.fst import calculate_fst_all_loci
    overall_fst = calculate_fst_all_loci(genetic_data)
    results['_overall'] = overall_fst

    return results


def add_overall_summary(results: Dict, genetic_data: GeneticData, analysis_type: str = 'general') -> Dict:
    """
    Add overall summary statistics across all populations

    Args:
        results: Per-population results
        genetic_data: GeneticData object
        analysis_type: Type of analysis ('general', 'hwe', 'heterozygosity', 'allele_frequencies')

    Returns:
        Results with _overall key added
    """
    from forstat.analysis.population.hardy_weinberg import hardy_weinberg_test
    from forstat.analysis.population.allele_frequencies import (
        calculate_allele_frequencies, calculate_heterozygosity
    )
    from forstat.analysis.forensic.match_probability import calculate_random_match_probability

    populations = [k for k in results.keys() if k.startswith('Population_')]

    if not populations:
        return results

    # Get all loci from first population
    first_pop = populations[0]
    loci = list(results[first_pop].get('loci', {}).keys())

    overall_loci = {}

    # For HWE, calculate actual HWE test on all samples (not averaged)
    if analysis_type == 'hwe':
        for locus in loci:
            # Calculate HWE on ALL samples combined
            hwe = hardy_weinberg_test(genetic_data.data, locus)

            # Also calculate mean Ho/He across populations for comparison
            ho_values = []
            he_values = []
            for pop in populations:
                locus_data = results[pop].get('loci', {}).get(locus, {})
                if locus_data:
                    ho_values.append(locus_data.get('Ho', 0))
                    he_values.append(locus_data.get('He', 0))

            overall_loci[locus] = {
                'Ho': hwe.get('Ho', 0) if 'Ho' in hwe else (np.mean(ho_values) if ho_values else 0),
                'He': hwe.get('He', 0) if 'He' in hwe else (np.mean(he_values) if he_values else 0),
                'Ho_mean': np.mean(ho_values) if ho_values else 0,
                'He_mean': np.mean(he_values) if he_values else 0,
                'hwe_p_value': hwe.get('p_value'),
                'hwe_status': hwe.get('hwe_status', 'N/A'),
                'chi_square': hwe.get('chi_square'),
                'df': hwe.get('df'),
                'n_genotypes': hwe.get('n_genotypes', 0),
                'n_populations': len(populations),
                'description': 'HWE test on all samples combined'
            }

    # For Heterozygosity, calculate actual values on all samples
    elif analysis_type == 'heterozygosity':
        for locus in loci:
            # Calculate heterozygosity on ALL samples
            het = calculate_heterozygosity(genetic_data.data, locus)
            allele_freqs = calculate_allele_frequencies(genetic_data.data, locus)

            # Also get mean PD across populations
            pd_values = []
            for pop in populations:
                locus_data = results[pop].get('loci', {}).get(locus, {})
                if locus_data:
                    pd_values.append(locus_data.get('PD', 0))

            overall_loci[locus] = {
                'Ho': het.get('Ho', 0),
                'He': het.get('He', 0),
                'Fis': het.get('Fis', 0),
                'n_alleles': len(allele_freqs),
                'PD': 1 - sum(p**2 for p in allele_freqs.values()),
                'PD_mean': np.mean(pd_values) if pd_values else 0,
                'n_populations': len(populations),
                'description': 'Calculated on all samples combined'
            }

    # For Allele Frequencies, calculate overall frequencies
    elif analysis_type == 'allele_frequencies':
        for locus in loci:
            # Calculate allele frequencies on ALL samples
            allele_freqs = calculate_allele_frequencies(genetic_data.data, locus)
            het = calculate_heterozygosity(genetic_data.data, locus)

            overall_loci[locus] = {
                'allele_frequencies': allele_freqs,
                'n_alleles': len(allele_freqs),
                'Ho': het.get('Ho', 0),
                'He': het.get('He', 0),
                'Fis': het.get('Fis', 0),
                'n_populations': len(populations),
                'description': 'Calculated on all samples combined'
            }

    # For general case, calculate means across populations
    else:
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
        'n_samples': len(genetic_data.data),
        'description': 'Statistics calculated on all samples combined' if analysis_type in ['hwe', 'heterozygosity', 'allele_frequencies'] else 'Mean statistics across all populations'
    }

    return results
