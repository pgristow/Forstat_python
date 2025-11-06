"""
Likelihood Ratio calculations for forensic genetics
"""
import numpy as np
from typing import Dict, Tuple, Optional
import pandas as pd

from forstat.utils.logger import get_logger

logger = get_logger(__name__)


def calculate_likelihood_ratio(
    genetic_data,
    evidence_profile: Dict[str, Tuple[str, str]],
    suspect_profile: Dict[str, Tuple[str, str]],
    theta: float = 0.01
) -> Dict:
    """
    Calculate Likelihood Ratio for DNA evidence

    LR = P(Evidence | Suspect is source) / P(Evidence | Random person is source)

    Args:
        genetic_data: GeneticData object
        evidence_profile: Evidence genotype profile {locus: (allele1, allele2)}
        suspect_profile: Suspect genotype profile {locus: (allele1, allele2)}
        theta: Population substructure coefficient (default 0.01)

    Returns:
        Dictionary with LR statistics
    """
    from forstat.analysis.population.allele_frequencies import calculate_allele_frequencies

    locus_LRs = {}
    combined_LR = 1.0

    for locus in evidence_profile.keys():
        if locus not in suspect_profile or locus not in genetic_data.loci:
            continue

        evidence_alleles = evidence_profile[locus]
        suspect_alleles = suspect_profile[locus]

        # Check if profiles match
        if sorted(evidence_alleles) != sorted(suspect_alleles):
            # No match - LR = 0
            locus_LRs[locus] = {
                'LR': 0,
                'match': False
            }
            combined_LR = 0
            continue

        # Get allele frequencies
        allele_freqs = calculate_allele_frequencies(genetic_data.data, locus)

        a1, a2 = evidence_alleles

        # P(Evidence | Suspect is source) = 1 if they match
        p_E_given_suspect = 1.0

        # P(Evidence | Random match) with theta correction
        p1 = allele_freqs.get(a1, 0.0001)  # Small value if allele not in database
        p2 = allele_freqs.get(a2, 0.0001)

        if a1 == a2:
            # Homozygote with theta correction
            p_E_given_random = (2 * theta + (1 - theta) * p1) * (3 * theta + (1 - 2 * theta) * p1) / ((1 + theta) * (1 + 2 * theta))
        else:
            # Heterozygote with theta correction
            p_E_given_random = 2 * (theta + (1 - theta) * p1) * (theta + (1 - theta) * p2) / ((1 + theta) * (1 + 2 * theta))

        if p_E_given_random > 0:
            LR = p_E_given_suspect / p_E_given_random
        else:
            LR = float('inf')

        locus_LRs[locus] = {
            'LR': LR,
            'match': True,
            'p_E_given_suspect': p_E_given_suspect,
            'p_E_given_random': p_E_given_random
        }

        combined_LR *= LR

    return {
        'combined_LR': combined_LR,
        'log10_LR': np.log10(combined_LR) if combined_LR > 0 else -np.inf,
        'locus_LRs': locus_LRs,
        'interpretation': interpret_LR(combined_LR)
    }


def interpret_LR(LR: float) -> str:
    """
    Interpret Likelihood Ratio according to ENFSI guidelines

    Args:
        LR: Likelihood ratio value

    Returns:
        String interpretation
    """
    log10_LR = np.log10(LR) if LR > 0 else -np.inf

    if log10_LR >= 6:
        return "Extremely strong support for prosecution hypothesis"
    elif log10_LR >= 4:
        return "Very strong support for prosecution hypothesis"
    elif log10_LR >= 2:
        return "Strong support for prosecution hypothesis"
    elif log10_LR >= 1:
        return "Moderate support for prosecution hypothesis"
    elif log10_LR >= 0:
        return "Limited support for prosecution hypothesis"
    elif log10_LR >= -1:
        return "Limited support for defense hypothesis"
    elif log10_LR >= -2:
        return "Moderate support for defense hypothesis"
    elif log10_LR >= -4:
        return "Strong support for defense hypothesis"
    else:
        return "Very strong support for defense hypothesis"


def calculate_mixture_LR(
    genetic_data,
    mixture_profile: Dict[str, list],
    poi_profile: Dict[str, Tuple[str, str]],
    n_contributors: int = 2
) -> Dict:
    """
    Calculate Likelihood Ratio for DNA mixture evidence

    Simplified calculation for mixed profiles

    Args:
        genetic_data: GeneticData object
        mixture_profile: Mixture alleles {locus: [allele1, allele2, ...]}
        poi_profile: Person of interest profile {locus: (allele1, allele2)}
        n_contributors: Number of contributors to mixture

    Returns:
        Dictionary with mixture LR statistics
    """
    from forstat.analysis.population.allele_frequencies import calculate_allele_frequencies

    locus_LRs = {}
    combined_LR = 1.0

    for locus in mixture_profile.keys():
        if locus not in poi_profile or locus not in genetic_data.loci:
            continue

        mixture_alleles = set(mixture_profile[locus])
        poi_alleles = poi_profile[locus]

        # Check if POI alleles are present in mixture
        poi_in_mixture = all(a in mixture_alleles for a in poi_alleles)

        if not poi_in_mixture:
            locus_LRs[locus] = {
                'LR': 0,
                'POI_in_mixture': False
            }
            combined_LR = 0
            continue

        # Simplified LR calculation
        # This is a basic implementation - real mixture analysis is much more complex
        allele_freqs = calculate_allele_frequencies(genetic_data.data, locus)

        # P(Mixture | POI is contributor)
        p_M_given_POI = 1.0  # Simplified

        # P(Mixture | Random person is contributor)
        # Probability that random person has alleles matching mixture
        p_random = 1.0
        for allele in poi_alleles:
            p_random *= allele_freqs.get(allele, 0.0001)

        if p_random > 0:
            LR = p_M_given_POI / p_random
        else:
            LR = float('inf')

        locus_LRs[locus] = {
            'LR': LR,
            'POI_in_mixture': True,
            'mixture_alleles': list(mixture_alleles),
            'POI_alleles': poi_alleles
        }

        combined_LR *= LR

    return {
        'combined_LR': combined_LR,
        'log10_LR': np.log10(combined_LR) if combined_LR > 0 else -np.inf,
        'locus_LRs': locus_LRs,
        'n_contributors': n_contributors,
        'note': 'Simplified mixture calculation - consult forensic expert for complex cases'
    }
