"""
Mitochondrial DNA analysis
"""
import numpy as np
from typing import Dict, List, Set
from collections import Counter

from forstat.utils.logger import get_logger

logger = get_logger(__name__)


def calculate_haplotype_diversity(sequences: List[str]) -> Dict:
    """
    Calculate haplotype diversity for mtDNA sequences

    Args:
        sequences: List of DNA sequences

    Returns:
        Dictionary with diversity statistics
    """
    if not sequences:
        return {'error': 'No sequences provided'}

    n = len(sequences)

    # Count unique haplotypes
    haplotype_counts = Counter(sequences)
    n_haplotypes = len(haplotype_counts)

    # Calculate haplotype diversity (H)
    # H = (n / (n - 1)) * (1 - sum(pi^2))
    # where pi is frequency of haplotype i
    sum_p_squared = sum((count / n) ** 2 for count in haplotype_counts.values())

    if n > 1:
        H = (n / (n - 1)) * (1 - sum_p_squared)
    else:
        H = 0

    # Variance of H
    var_H = ((2 / n) * (n - 2) / (n - 1)) * (1 - sum_p_squared) ** 2

    # Standard error
    se_H = np.sqrt(var_H)

    return {
        'haplotype_diversity': H,
        'variance': var_H,
        'standard_error': se_H,
        'n_sequences': n,
        'n_unique_haplotypes': n_haplotypes,
        'haplotype_frequencies': {hap: count / n for hap, count in haplotype_counts.items()}
    }


def calculate_nucleotide_diversity(sequences: List[str]) -> Dict:
    """
    Calculate nucleotide diversity (π) for mtDNA sequences

    π = average number of nucleotide differences per site between sequences

    Args:
        sequences: List of DNA sequences (must be aligned)

    Returns:
        Dictionary with nucleotide diversity statistics
    """
    if not sequences or len(sequences) < 2:
        return {'error': 'Need at least 2 sequences'}

    # Check if sequences are same length (aligned)
    seq_lengths = [len(s) for s in sequences]
    if len(set(seq_lengths)) > 1:
        return {'error': 'Sequences must be aligned (same length)'}

    n = len(sequences)
    seq_length = len(sequences[0])

    # Calculate pairwise differences
    total_differences = 0
    n_comparisons = 0

    for i in range(n):
        for j in range(i + 1, n):
            differences = sum(s1 != s2 for s1, s2 in zip(sequences[i], sequences[j]))
            total_differences += differences
            n_comparisons += 1

    # Average pairwise differences
    if n_comparisons > 0:
        avg_pairwise_diff = total_differences / n_comparisons
    else:
        avg_pairwise_diff = 0

    # Nucleotide diversity (per site)
    pi = avg_pairwise_diff / seq_length if seq_length > 0 else 0

    return {
        'nucleotide_diversity': pi,
        'average_pairwise_differences': avg_pairwise_diff,
        'n_sequences': n,
        'sequence_length': seq_length,
        'n_comparisons': n_comparisons
    }


def identify_polymorphic_sites(sequences: List[str]) -> Dict:
    """
    Identify polymorphic (variable) sites in mtDNA sequences

    Args:
        sequences: List of aligned DNA sequences

    Returns:
        Dictionary with polymorphic site information
    """
    if not sequences:
        return {'error': 'No sequences provided'}

    # Check alignment
    seq_lengths = [len(s) for s in sequences]
    if len(set(seq_lengths)) > 1:
        return {'error': 'Sequences must be aligned'}

    seq_length = len(sequences[0])
    polymorphic_sites = []
    mutations = {}

    # Check each position
    for pos in range(seq_length):
        bases_at_pos = [seq[pos] for seq in sequences if pos < len(seq)]
        unique_bases = set(bases_at_pos)

        # Remove gaps
        unique_bases.discard('-')
        unique_bases.discard('N')

        if len(unique_bases) > 1:
            # Polymorphic site
            polymorphic_sites.append(pos)

            # Count each variant
            base_counts = Counter(bases_at_pos)
            mutations[pos] = {
                'position': pos + 1,  # 1-indexed
                'variants': dict(base_counts),
                'n_variants': len(unique_bases)
            }

    # Calculate segregating sites (S)
    S = len(polymorphic_sites)

    # Watterson's theta (θw) - estimate of mutation rate
    # θw = S / a_n where a_n = sum(1/i) for i=1 to n-1
    n = len(sequences)
    if n > 1:
        a_n = sum(1 / i for i in range(1, n))
        theta_w = S / a_n if a_n > 0 else 0
    else:
        theta_w = 0

    return {
        'n_polymorphic_sites': S,
        'polymorphic_positions': polymorphic_sites,
        'mutations': mutations,
        'theta_watterson': theta_w,
        'proportion_polymorphic': S / seq_length if seq_length > 0 else 0
    }


def calculate_tajimas_d(sequences: List[str]) -> Dict:
    """
    Calculate Tajima's D statistic

    Tests for neutrality - deviations indicate selection or population change

    Args:
        sequences: List of aligned DNA sequences

    Returns:
        Dictionary with Tajima's D statistics
    """
    if len(sequences) < 4:
        return {'error': 'Need at least 4 sequences for Tajima\'s D'}

    n = len(sequences)

    # Calculate nucleotide diversity (π)
    pi_result = calculate_nucleotide_diversity(sequences)
    if 'error' in pi_result:
        return pi_result

    pi = pi_result['average_pairwise_differences']

    # Calculate segregating sites (S)
    poly_result = identify_polymorphic_sites(sequences)
    S = poly_result['n_polymorphic_sites']

    # Watterson's theta
    a1 = sum(1 / i for i in range(1, n))
    theta_w = S / a1 if a1 > 0 else 0

    # Tajima's D = (π - θw) / sqrt(Var(π - θw))
    # Calculate variance
    a2 = sum(1 / (i ** 2) for i in range(1, n))

    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n ** 2 + n + 3) / (9 * n * (n - 1))

    c1 = b1 - (1 / a1)
    c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1 ** 2))

    e1 = c1 / a1
    e2 = c2 / (a1 ** 2 + a2)

    var_d = e1 * S + e2 * S * (S - 1)
    std_d = np.sqrt(var_d) if var_d > 0 else 0

    # Calculate D
    if std_d > 0:
        D = (pi - theta_w) / std_d
    else:
        D = 0

    # Interpret
    if D < -2:
        interpretation = "Negative D: Excess rare variants (purifying selection or population expansion)"
    elif D > 2:
        interpretation = "Positive D: Excess intermediate-frequency variants (balancing selection or population bottleneck)"
    else:
        interpretation = "D not significant: Consistent with neutral evolution"

    return {
        'tajimas_d': D,
        'pi': pi,
        'theta_w': theta_w,
        'segregating_sites': S,
        'interpretation': interpretation
    }


def match_mtdna_sequences(query: str, reference: str, allow_mismatches: int = 0) -> Dict:
    """
    Compare mtDNA sequences for forensic matching

    Args:
        query: Query sequence
        reference: Reference sequence
        allow_mismatches: Number of allowed mismatches

    Returns:
        Dictionary with match statistics
    """
    if len(query) != len(reference):
        return {'error': 'Sequences must be same length'}

    # Count differences
    differences = []
    for i, (q, r) in enumerate(zip(query, reference)):
        if q != r and q != '-' and r != '-' and q != 'N' and r != 'N':
            differences.append({
                'position': i + 1,
                'query': q,
                'reference': r
            })

    n_differences = len(differences)

    # Determine match
    if n_differences == 0:
        match_status = "Perfect match"
    elif n_differences <= allow_mismatches:
        match_status = f"Match (within {allow_mismatches} mismatch tolerance)"
    else:
        match_status = "No match"

    # Calculate similarity
    total_compared = sum(1 for q, r in zip(query, reference) if q != '-' and r != '-')
    similarity = (total_compared - n_differences) / total_compared if total_compared > 0 else 0

    return {
        'match_status': match_status,
        'n_differences': n_differences,
        'differences': differences,
        'similarity': similarity,
        'percent_identity': similarity * 100
    }
