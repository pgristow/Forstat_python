"""
Data models for genetic data
"""
from dataclasses import dataclass, field
from typing import List, Dict, Optional
import pandas as pd


@dataclass
class GeneticData:
    """Container for genetic data"""

    title: str
    file_path: str
    file_name: str
    file_type: str

    # Loci information
    loci: List[str]
    n_loci: int

    # Sample information
    n_samples: int
    n_populations: int

    # Raw data
    data: pd.DataFrame
    populations: List[Dict] = field(default_factory=list)

    # Metadata
    metadata: Dict = field(default_factory=dict)

    def get_locus_data(self, locus: str) -> pd.DataFrame:
        """Get data for a specific locus"""
        cols = ['population', 'sample_name', f'{locus}_1', f'{locus}_2']
        return self.data[cols].copy()

    def get_population_data(self, pop_id: int) -> pd.DataFrame:
        """Get data for a specific population"""
        return self.data[self.data['population'] == pop_id].copy()

    def get_alleles(self, locus: str) -> List[str]:
        """Get all unique alleles for a locus"""
        alleles = set()
        alleles.update(self.data[f'{locus}_1'].unique())
        alleles.update(self.data[f'{locus}_2'].unique())
        alleles.discard('0')  # Remove missing data
        return sorted(list(alleles))

    def __str__(self):
        return (f"GeneticData: {self.file_name}\n"
                f"  Samples: {self.n_samples}\n"
                f"  Loci: {self.n_loci}\n"
                f"  Populations: {self.n_populations}")
