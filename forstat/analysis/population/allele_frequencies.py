"""
Allele frequency calculations
"""
from collections import Counter
from typing import Dict, List
import pandas as pd

from forstat.utils.logger import get_logger

logger = get_logger(__name__)


def calculate_allele_frequencies(data: pd.DataFrame, locus: str) -> Dict[str, float]:
    """
    Calculate allele frequencies for a locus

    Args:
        data: DataFrame with genetic data
        locus: Locus name

    Returns:
        Dictionary mapping allele -> frequency
    """
    allele_col1 = f"{locus}_1"
    allele_col2 = f"{locus}_2"

    if allele_col1 not in data.columns or allele_col2 not in data.columns:
        logger.error(f"Locus {locus} not found in data")
        return {}

    # Count all alleles
    alleles = []
    for _, row in data.iterrows():
        a1 = str(row[allele_col1])
        a2 = str(row[allele_col2])

        # Skip missing data
        if a1 != '0' and a1 != '00' and a1 != '000':
            alleles.append(a1)
        if a2 != '0' and a2 != '00' and a2 != '000':
            alleles.append(a2)

    if not alleles:
        return {}

    # Calculate frequencies
    total = len(alleles)
    counts = Counter(alleles)
    frequencies = {allele: count / total for allele, count in counts.items()}

    return frequencies


def calculate_heterozygosity(data: pd.DataFrame, locus: str) -> Dict[str, float]:
    """
    Calculate observed and expected heterozygosity

    Args:
        data: DataFrame with genetic data
        locus: Locus name

    Returns:
        Dictionary with Ho (observed) and He (expected)
    """
    allele_col1 = f"{locus}_1"
    allele_col2 = f"{locus}_2"

    # Observed heterozygosity
    heterozygotes = 0
    total_genotypes = 0

    for _, row in data.iterrows():
        a1 = str(row[allele_col1])
        a2 = str(row[allele_col2])

        # Skip missing data
        if a1 in ['0', '00', '000'] or a2 in ['0', '00', '000']:
            continue

        total_genotypes += 1
        if a1 != a2:
            heterozygotes += 1

    Ho = heterozygotes / total_genotypes if total_genotypes > 0 else 0

    # Expected heterozygosity (gene diversity)
    frequencies = calculate_allele_frequencies(data, locus)
    He = 1.0 - sum(f ** 2 for f in frequencies.values())

    return {
        'Ho': Ho,
        'He': He,
        'Fis': 1 - (Ho / He) if He > 0 else 0
    }


def calculate_summary_statistics(genetic_data) -> Dict:
    """
    Calculate summary statistics for all loci

    Args:
        genetic_data: GeneticData object

    Returns:
        Dictionary with statistics per locus
    """
    results = {}

    for locus in genetic_data.loci:
        logger.info(f"Calculating statistics for {locus}")

        frequencies = calculate_allele_frequencies(genetic_data.data, locus)
        heterozygosity = calculate_heterozygosity(genetic_data.data, locus)

        results[locus] = {
            'allele_frequencies': frequencies,
            'n_alleles': len(frequencies),
            'Ho': heterozygosity['Ho'],
            'He': heterozygosity['He'],
            'Fis': heterozygosity['Fis']
        }

    return results
