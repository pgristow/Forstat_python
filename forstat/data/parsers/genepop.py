"""
GenePop file format parser
"""
import re
from typing import Dict, List, Tuple, Optional
from pathlib import Path
import pandas as pd

from forstat.utils.logger import get_logger

logger = get_logger(__name__)


class GenepopParser:
    """Parser for GenePop format files"""

    def __init__(self, file_path: str):
        self.file_path = Path(file_path)
        self.title = ""
        self.loci = []
        self.populations = []
        self.data = None

    def parse(self) -> Dict:
        """
        Parse GenePop file and return structured data

        Returns:
            Dict containing:
                - title: str
                - loci: List[str]
                - populations: List[Dict]
                - n_populations: int
                - n_samples: int
                - n_loci: int
                - data: pandas.DataFrame
        """
        logger.info(f"Parsing GenePop file: {self.file_path}")

        try:
            with open(self.file_path, 'r', encoding='utf-8', errors='ignore') as f:
                lines = [line.strip() for line in f.readlines() if line.strip()]

            if not lines:
                raise ValueError("File is empty")

            # First line is the title/comment
            self.title = lines[0]
            logger.info(f"Title: {self.title}")

            # Parse loci names (until we hit "POP")
            idx = 1
            while idx < len(lines) and lines[idx].upper() != 'POP':
                locus = lines[idx].strip()
                if locus:
                    self.loci.append(locus)
                idx += 1

            logger.info(f"Found {len(self.loci)} loci: {self.loci}")

            # Parse populations and samples
            current_pop = []
            pop_number = 0

            for i in range(idx, len(lines)):
                line = lines[i].strip()

                if line.upper() == 'POP':
                    # Save previous population if it exists
                    if current_pop:
                        self.populations.append({
                            'pop_id': pop_number,
                            'samples': current_pop
                        })
                        pop_number += 1
                    current_pop = []
                else:
                    # Parse sample line
                    sample = self._parse_sample_line(line)
                    if sample:
                        current_pop.append(sample)

            # Add last population
            if current_pop:
                self.populations.append({
                    'pop_id': pop_number,
                    'samples': current_pop
                })

            logger.info(f"Found {len(self.populations)} populations")

            # Convert to DataFrame for easier analysis
            self._create_dataframe()

            # Create result dictionary
            result = {
                'title': self.title,
                'loci': self.loci,
                'populations': self.populations,
                'n_populations': len(self.populations),
                'n_samples': sum(len(pop['samples']) for pop in self.populations),
                'n_loci': len(self.loci),
                'data': self.data,
                'file_path': str(self.file_path),
                'file_name': self.file_path.name
            }

            logger.info(f"Successfully parsed: {result['n_samples']} samples, "
                       f"{result['n_loci']} loci, {result['n_populations']} populations")

            return result

        except Exception as e:
            logger.error(f"Error parsing GenePop file: {e}")
            raise

    def _parse_sample_line(self, line: str) -> Optional[Dict]:
        """Parse a single sample line"""
        if not line or line.upper() == 'POP':
            return None

        try:
            # Split on comma or tab
            parts = re.split(r'[,\t]', line, maxsplit=1)
            if len(parts) != 2:
                return None

            sample_name = parts[0].strip()
            genotypes_str = parts[1].strip()

            # Split genotypes by whitespace
            genotypes = genotypes_str.split()

            if len(genotypes) != len(self.loci):
                logger.warning(f"Sample {sample_name} has {len(genotypes)} genotypes, "
                             f"expected {len(self.loci)}")

            # Parse each genotype into two alleles
            parsed_genotypes = {}
            for i, (locus, genotype) in enumerate(zip(self.loci, genotypes)):
                alleles = self._parse_genotype(genotype)
                parsed_genotypes[locus] = alleles

            return {
                'sample_name': sample_name,
                'genotypes': parsed_genotypes
            }

        except Exception as e:
            logger.warning(f"Could not parse sample line: {line} - {e}")
            return None

    def _parse_genotype(self, genotype: str) -> Tuple[str, str]:
        """
        Parse genotype string into two alleles

        GenePop format: 0404 or 0000 (missing data)
        Returns: tuple of (allele1, allele2)
        """
        genotype = genotype.strip()

        # Handle missing data
        if genotype in ['0000', '000000', '00', '0', '-', 'NA', '']:
            return ('0', '0')

        # Standard format: 2-digit or 3-digit alleles
        if len(genotype) == 4:
            # 2-digit alleles: 0404 -> (04, 04)
            return (genotype[:2], genotype[2:])
        elif len(genotype) == 6:
            # 3-digit alleles: 010101 -> (010, 010)
            return (genotype[:3], genotype[3:])
        elif len(genotype) == 2:
            # Single digit each: 12 -> (1, 2)
            return (genotype[0], genotype[1])
        else:
            logger.warning(f"Unexpected genotype format: {genotype}")
            return ('0', '0')

    def _create_dataframe(self):
        """Create a pandas DataFrame from parsed data"""
        rows = []

        for pop in self.populations:
            pop_id = pop['pop_id']
            for sample in pop['samples']:
                row = {
                    'population': pop_id,
                    'sample_name': sample['sample_name']
                }

                # Add each locus as two columns (allele1, allele2)
                for locus in self.loci:
                    if locus in sample['genotypes']:
                        allele1, allele2 = sample['genotypes'][locus]
                        row[f"{locus}_1"] = allele1
                        row[f"{locus}_2"] = allele2
                    else:
                        row[f"{locus}_1"] = '0'
                        row[f"{locus}_2"] = '0'

                rows.append(row)

        self.data = pd.DataFrame(rows)
        logger.info(f"Created DataFrame with shape: {self.data.shape}")


def parse_genepop(file_path: str) -> Dict:
    """
    Convenience function to parse GenePop file

    Args:
        file_path: Path to GenePop file

    Returns:
        Dictionary with parsed data
    """
    parser = GenepopParser(file_path)
    return parser.parse()
