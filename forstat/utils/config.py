"""
Configuration management for Forstat application
"""
import os
import json
from pathlib import Path

class Config:
    """Application configuration manager"""

    # Application Info
    APP_NAME = "Forstat"
    APP_VERSION = "0.1.0"
    APP_TITLE = "Forstat - Forensic Statistics Analysis"

    # Directories
    BASE_DIR = Path(__file__).parent.parent.parent
    RESOURCES_DIR = BASE_DIR / "resources"
    ICONS_DIR = RESOURCES_DIR / "icons"
    THEMES_DIR = RESOURCES_DIR / "themes"
    TEMPLATES_DIR = RESOURCES_DIR / "templates"

    # User directories
    USER_DIR = Path.home() / ".forstat"
    OUTPUT_DIR = USER_DIR / "output"
    TEMP_DIR = USER_DIR / "temp"

    # File formats supported
    SUPPORTED_FORMATS = {
        'genepop': ['.gen', '.txt'],
        'excel': ['.xlsx', '.xls'],
        'csv': ['.csv'],
        'fasta': ['.fasta', '.fa', '.fna'],
        'vcf': ['.vcf']
    }

    # Analysis types
    ANALYSIS_TYPES = {
        'population': [
            'Hardy-Weinberg Equilibrium',
            'Fixation Index (Fst)',
            'Heterozygosity',
            'Allele Frequencies',
            'AMOVA'
        ],
        'forensic': [
            'Match Probability and Power of Discrimination',
            'Likelihood Ratio',
            'Paternity Index',
            'Kinship Analysis'
        ],
        'str': [
            'STR Marker Analysis'
        ],
        'mtdna': [
            'Haplotype Diversity',
            'Nucleotide Diversity',
            'Phylogenetic Analysis'
        ]
    }

    @classmethod
    def init_dirs(cls):
        """Initialize user directories"""
        cls.USER_DIR.mkdir(exist_ok=True)
        cls.OUTPUT_DIR.mkdir(exist_ok=True)
        cls.TEMP_DIR.mkdir(exist_ok=True)

    @classmethod
    def get_config_path(cls):
        """Get path to user configuration file"""
        return cls.USER_DIR / "config.json"

    @classmethod
    def load_user_config(cls):
        """Load user configuration"""
        config_path = cls.get_config_path()
        if config_path.exists():
            with open(config_path, 'r') as f:
                return json.load(f)
        return {}

    @classmethod
    def save_user_config(cls, config):
        """Save user configuration"""
        config_path = cls.get_config_path()
        with open(config_path, 'w') as f:
            json.dump(config, f, indent=4)
