"""
Kinship and relatedness analysis
"""
import numpy as np
from typing import Dict, Tuple, Optional
import pandas as pd

from forstat.utils.logger import get_logger

logger = get_logger(__name__)


def calculate_paternity_index(
    genetic_data,
    child_profile: Dict[str, Tuple[str, str]],
    alleged_father_profile: Dict[str, Tuple[str, str]],
    mother_profile: Optional[Dict[str, Tuple[str, str]]] = None
) -> Dict:
    """
    Calculate Paternity Index (PI)

    PI = P(child's genotype | alleged father is true father) /
         P(child's genotype | random man is father)

    Args:
        genetic_data: GeneticData object
        child_profile: Child's genotype {locus: (allele1, allele2)}
        alleged_father_profile: Alleged father's genotype {locus: (allele1, allele2)}
        mother_profile: Mother's genotype (optional, improves accuracy)

    Returns:
        Dictionary with PI statistics
    """
    from forstat.analysis.population.allele_frequencies import calculate_allele_frequencies

    locus_PIs = {}
    combined_PI = 1.0
    exclusions = []

    for locus in child_profile.keys():
        if locus not in alleged_father_profile or locus not in genetic_data.loci:
            continue

        child_alleles = set(child_profile[locus])
        af_alleles = set(alleged_father_profile[locus])

        # Determine paternal allele from child
        if mother_profile and locus in mother_profile:
            mother_alleles = set(mother_profile[locus])
            # Paternal allele = child allele not from mother
            paternal_candidates = child_alleles - mother_alleles
            if not paternal_candidates:
                # Child could have gotten both alleles from mother
                paternal_candidates = child_alleles
        else:
            # Without mother, both child alleles are candidates
            paternal_candidates = child_alleles

        # Check if alleged father has any paternal allele
        shared_alleles = paternal_candidates & af_alleles

        if not shared_alleles:
            # Exclusion
            locus_PIs[locus] = {
                'PI': 0,
                'excluded': True,
                'child_alleles': list(child_alleles),
                'af_alleles': list(af_alleles)
            }
            exclusions.append(locus)
            combined_PI = 0
            continue

        # Calculate PI for this locus
        allele_freqs = calculate_allele_frequencies(genetic_data.data, locus)

        # Use most probable paternal allele
        paternal_allele = list(shared_alleles)[0]

        # P(transmission | alleged father is true father)
        # If father is heterozygous for this allele: 0.5, if homozygous: 1.0
        if alleged_father_profile[locus][0] == alleged_father_profile[locus][1]:
            p_transmission = 1.0
        else:
            p_transmission = 0.5

        # P(child has this allele | random man is father)
        p_random = allele_freqs.get(paternal_allele, 0.0001)

        if p_random > 0:
            PI = p_transmission / p_random
        else:
            PI = float('inf')

        locus_PIs[locus] = {
            'PI': PI,
            'excluded': False,
            'paternal_allele': paternal_allele,
            'p_transmission': p_transmission,
            'p_random': p_random
        }

        combined_PI *= PI

    # Calculate probability of paternity
    if combined_PI > 0:
        # Assuming prior probability of 0.5
        prob_paternity = combined_PI / (combined_PI + 1)
    else:
        prob_paternity = 0.0

    return {
        'combined_PI': combined_PI,
        'probability_of_paternity': prob_paternity,
        'locus_PIs': locus_PIs,
        'exclusions': exclusions,
        'is_excluded': len(exclusions) > 0,
        'interpretation': interpret_PI(combined_PI)
    }


def interpret_PI(PI: float) -> str:
    """
    Interpret Paternity Index

    Args:
        PI: Combined Paternity Index

    Returns:
        String interpretation
    """
    if PI == 0:
        return "Excluded as biological father"
    elif PI < 1:
        return "Not likely to be biological father"
    elif PI < 10:
        return "Weak support for paternity"
    elif PI < 100:
        return "Moderate support for paternity"
    elif PI < 1000:
        return "Strong support for paternity"
    elif PI < 10000:
        return "Very strong support for paternity"
    else:
        return "Extremely strong support for paternity"


def calculate_ibd_allele_sharing(
    profile1: Dict[str, Tuple[str, str]],
    profile2: Dict[str, Tuple[str, str]]
) -> Dict:
    """
    Calculate Identity-By-Descent (IBD) allele sharing

    Args:
        profile1: First individual's profile
        profile2: Second individual's profile

    Returns:
        Dictionary with IBD statistics
    """
    ibd_scores = {}
    total_alleles_compared = 0
    total_shared = 0

    for locus in profile1.keys():
        if locus not in profile2:
            continue

        alleles1 = set(profile1[locus])
        alleles2 = set(profile2[locus])

        # Count shared alleles
        shared = len(alleles1 & alleles2)
        total = len(alleles1) + len(alleles2)

        ibd_scores[locus] = {
            'shared_alleles': shared,
            'total_possible': 2,  # Maximum 2 shared for STRs
            'proportion': shared / 2
        }

        total_alleles_compared += 2
        total_shared += shared

    # Overall proportion of shared alleles
    if total_alleles_compared > 0:
        overall_proportion = total_shared / total_alleles_compared
    else:
        overall_proportion = 0

    # Interpret relationship
    relationship = interpret_relatedness(overall_proportion)

    return {
        'overall_proportion_shared': overall_proportion,
        'locus_scores': ibd_scores,
        'n_loci_compared': len(ibd_scores),
        'predicted_relationship': relationship
    }


def interpret_relatedness(proportion_shared: float) -> str:
    """
    Interpret relationship based on proportion of shared alleles

    Args:
        proportion_shared: Proportion of alleles shared (0-1)

    Returns:
        String describing predicted relationship
    """
    if proportion_shared >= 0.95:
        return "Identical twins or same individual"
    elif proportion_shared >= 0.65:
        return "Parent-offspring or full siblings"
    elif proportion_shared >= 0.40:
        return "Half-siblings or grandparent-grandchild"
    elif proportion_shared >= 0.25:
        return "First cousins or uncle/aunt-nephew/niece"
    elif proportion_shared >= 0.15:
        return "Second cousins or distant relatives"
    else:
        return "Unrelated"
