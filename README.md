# Forstat - Forensic Statistics Application

A comprehensive Windows desktop application for forensic and population genetics analysis.

## Features

### Data Import
- **Multiple Format Support**: GenePop, Excel (.xlsx, .xls), CSV, FASTA
- **Drag & Drop Interface**: Easy file import
- **Data Validation**: Automatic format detection and validation
- **Data Preview**: View data before analysis

### Analysis Types

#### Population Genetics
- Hardy-Weinberg Equilibrium (HWE)
- Fixation Index (Fst)
- Heterozygosity (He, Ho)
- Allele Frequencies
- AMOVA (Analysis of Molecular Variance)
- Population Structure Analysis

#### Forensic Statistics
- Match Probability
- Likelihood Ratio (LR)
- Power of Discrimination (PD)
- Paternity Index (PI)
- Kinship Analysis

#### STR Analysis
- Allele Frequency Analysis
- Stutter Analysis
- Peak Height Ratio

#### Mitochondrial DNA
- Haplotype Diversity
- Nucleotide Diversity
- Phylogenetic Analysis

### Results & Export
- Interactive result tables
- Data visualization
- Export to Excel, PDF, and CSV
- Customizable reports

## Installation

### Prerequisites
- Python 3.9 or higher
- Windows 10 or higher

### Setup

1. Clone the repository:
```bash
git clone https://github.com/pgristow/Forstat_python.git
cd Forstat_python
```

2. Create a virtual environment:
```bash
python -m venv venv
venv\Scripts\activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### Running the Application

```bash
python main.py
```

### Basic Workflow

1. **Upload Data**
   - Click "Browse" or drag & drop your data file
   - Select the appropriate file format (or use auto-detect)
   - Preview your data
   - Click "Load Data"

2. **Run Analysis**
   - Select desired analyses from the available categories
   - Configure analysis parameters if needed
   - Click "Run Analysis"
   - Monitor progress

3. **View Results**
   - Review results in organized tabs
   - Add notes to specific analyses
   - Export results in your preferred format

## Project Structure

```
Forstat_python/
├── main.py                      # Application entry point
├── requirements.txt             # Python dependencies
├── setup.py                     # Installation script
├── forstat/                     # Main package
│   ├── gui/                     # GUI components
│   │   ├── main_window.py      # Main application window
│   │   ├── upload_page.py      # Data upload interface
│   │   ├── analysis_page.py    # Analysis selection
│   │   ├── output_page.py      # Results display
│   │   └── styles.py           # UI styling
│   ├── data/                    # Data handling
│   │   └── parsers/            # File format parsers
│   ├── analysis/               # Statistical analysis
│   │   ├── population/         # Population genetics
│   │   └── forensic/           # Forensic statistics
│   ├── reporting/              # Export functionality
│   └── utils/                  # Utilities
└── tests/                      # Tests and sample data
    └── sample_data/            # Sample files
```

## Development Status

**Version**: 0.1.0 (Alpha)

### Completed
- ✅ Project structure and configuration
- ✅ Modern GUI framework with PyQt6
- ✅ Upload page with drag & drop
- ✅ Analysis page with multiple categories
- ✅ Output page with export options
- ✅ Navigation and workflow management

### In Progress
- 🚧 Data parsers implementation
- 🚧 Statistical analysis modules
- 🚧 Result visualization
- 🚧 Export functionality

### Planned
- 📋 Additional analysis types
- 📋 Advanced visualization options
- 📋 Batch processing
- 📋 Database integration
- 📋 User preferences and themes

## Dependencies

Main dependencies:
- PyQt6 - GUI framework
- pandas - Data manipulation
- numpy - Numerical computing
- scipy - Scientific computing
- scikit-allel - Population genetics
- biopython - Sequence analysis
- matplotlib/plotly - Visualization
- reportlab - PDF generation

See `requirements.txt` for complete list.

## Building Standalone Executable

To create a standalone Windows executable:

```bash
pip install pyinstaller
pyinstaller --onefile --windowed --name Forstat main.py
```

The executable will be in the `dist/` directory.

## Contributing

This project is under active development. Contributions are welcome!

## License

MIT License - See LICENSE file for details

## Support

For issues, questions, or suggestions, please open an issue on GitHub.

## Authors

Forstat Development Team

## Acknowledgments

- Built for forensic and population genetics researchers
- Designed for Windows platform
- Modern UI inspired by Material Design principles
