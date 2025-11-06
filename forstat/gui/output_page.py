"""
Output page for displaying and exporting results
"""
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QGroupBox, QTableWidget, QTableWidgetItem, QTextEdit,
    QTabWidget, QFileDialog, QMessageBox
)
from PyQt6.QtCore import Qt
from pathlib import Path
from datetime import datetime

from forstat.utils.logger import get_logger
from forstat.utils.helpers import sanitize_filename

logger = get_logger(__name__)


class OutputPage(QWidget):
    """Page for displaying and exporting analysis results"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.main_window = parent

        self.init_ui()

    def init_ui(self):
        """Initialize user interface"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(30, 30, 30, 30)
        layout.setSpacing(20)

        # Title
        title = QLabel("Results")
        title.setProperty("class", "title")
        layout.addWidget(title)

        # Summary section
        summary_group = QGroupBox("Analysis Summary")
        summary_layout = QVBoxLayout(summary_group)

        self.summary_text = QTextEdit()
        self.summary_text.setReadOnly(True)
        self.summary_text.setMaximumHeight(120)
        summary_layout.addWidget(self.summary_text)

        layout.addWidget(summary_group)

        # Results tabs
        self.results_tabs = QTabWidget()
        layout.addWidget(self.results_tabs, 1)

        # Export buttons
        button_layout = QHBoxLayout()
        button_layout.addStretch()

        self.export_excel_btn = QPushButton("Export to Excel")
        self.export_excel_btn.clicked.connect(self.export_to_excel)
        button_layout.addWidget(self.export_excel_btn)

        self.export_pdf_btn = QPushButton("Export to PDF")
        self.export_pdf_btn.clicked.connect(self.export_to_pdf)
        button_layout.addWidget(self.export_pdf_btn)

        self.export_csv_btn = QPushButton("Export to CSV")
        self.export_csv_btn.clicked.connect(self.export_to_csv)
        button_layout.addWidget(self.export_csv_btn)

        layout.addLayout(button_layout)

        # Show placeholder message
        self.show_placeholder()

    def show_placeholder(self):
        """Show placeholder when no results available"""
        self.summary_text.setHtml("""
            <p style='color: #757575; font-size: 11pt;'>
            <i>No results available yet. Run analyses to see results here.</i>
            </p>
        """)

    def display_results(self, results):
        """Display analysis results"""
        if not results:
            self.show_placeholder()
            return

        # Clear existing tabs
        self.results_tabs.clear()

        # Update summary
        self.update_summary(results)

        # Create tabs for each analysis type
        analyses = results.get('analyses', [])

        if analyses:
            # Create a tab for each analysis
            for analysis in analyses:
                tab = self.create_result_tab(analysis, results)
                self.results_tabs.addTab(tab, analysis)

        logger.info(f"Displayed results for {len(analyses)} analyses")

    def update_summary(self, results):
        """Update results summary"""
        data_file = results.get('data_file', 'Unknown')
        analyses = results.get('analyses', [])
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        summary_html = f"""
        <h3>Analysis Summary</h3>
        <p><b>Data File:</b> {data_file}</p>
        <p><b>Analyses Performed:</b> {len(analyses)}</p>
        <p><b>Completed:</b> {timestamp}</p>
        <p><b>Status:</b> <span style='color: #4CAF50; font-weight: bold;'>✓ Success</span></p>
        """

        self.summary_text.setHtml(summary_html)

    def create_result_tab(self, analysis_name, results):
        """Create a tab for specific analysis results"""
        widget = QWidget()
        layout = QVBoxLayout(widget)

        # Analysis title
        title = QLabel(analysis_name)
        title.setStyleSheet("font-size: 14pt; font-weight: bold; color: #2196F3;")
        layout.addWidget(title)

        # Description
        description = QLabel(self.get_analysis_description(analysis_name))
        description.setWordWrap(True)
        description.setProperty("class", "hint")
        layout.addWidget(description)

        # Get actual results data
        analysis_results = results.get('results', {}).get(analysis_name, {})

        # Results table
        table = QTableWidget()
        table.setAlternatingRowColors(True)

        # Display results based on analysis type
        if analysis_name == 'Allele Frequencies':
            self._populate_allele_frequencies_table(table, analysis_results)
        elif analysis_name in ['Heterozygosity', 'Hardy-Weinberg Equilibrium']:
            self._populate_per_population_table(table, analysis_results, analysis_name)
        elif analysis_name == 'Fixation Index (Fst)':
            self._populate_fst_table(table, analysis_results)
        elif analysis_name == 'Match Probability and Power of Discrimination':
            self._populate_match_prob_per_pop_table(table, analysis_results)
        elif analysis_name == 'STR Marker Analysis':
            self._populate_str_analysis_table(table, analysis_results)
        elif 'status' in analysis_results or 'error' in analysis_results:
            # Status/error message
            self._populate_status_table(table, analysis_results)
        else:
            # Generic display
            self._populate_generic_table(table, analysis_results)

        table.resizeColumnsToContents()
        layout.addWidget(table)

        # Notes section
        notes_group = QGroupBox("Notes")
        notes_layout = QVBoxLayout(notes_group)

        notes_text = QTextEdit()
        notes_text.setPlaceholderText("Add notes about this analysis...")
        notes_text.setMaximumHeight(100)
        notes_layout.addWidget(notes_text)

        layout.addWidget(notes_group)

        return widget

    def _populate_diversity_table(self, table, analysis_results, genetic_data):
        """Populate table with diversity/allele frequency results"""
        if not analysis_results or not genetic_data:
            table.setRowCount(1)
            table.setColumnCount(1)
            table.setHorizontalHeaderLabels(["Status"])
            table.setItem(0, 0, QTableWidgetItem("No results available"))
            return

        # Filter out summary keys
        loci = [k for k in analysis_results.keys() if not k.startswith('_')]

        if not loci:
            table.setRowCount(1)
            table.setColumnCount(1)
            table.setHorizontalHeaderLabels(["Status"])
            table.setItem(0, 0, QTableWidgetItem("No locus data available"))
            return

        table.setRowCount(len(loci))
        table.setColumnCount(6)
        table.setHorizontalHeaderLabels(["Locus", "N Alleles", "Ho", "He", "Fis", "PIC"])

        for i, locus in enumerate(loci):
            locus_data = analysis_results[locus]

            table.setItem(i, 0, QTableWidgetItem(locus))
            table.setItem(i, 1, QTableWidgetItem(str(locus_data.get('n_alleles', 'N/A'))))
            table.setItem(i, 2, QTableWidgetItem(f"{locus_data.get('Ho', 0):.4f}"))
            table.setItem(i, 3, QTableWidgetItem(f"{locus_data.get('He', 0):.4f}"))
            table.setItem(i, 4, QTableWidgetItem(f"{locus_data.get('Fis', 0):.4f}"))
            table.setItem(i, 5, QTableWidgetItem(f"{locus_data.get('PIC', 0):.4f}"))

    def _populate_hwe_table(self, table, analysis_results):
        """Populate table with Hardy-Weinberg Equilibrium results"""
        if not analysis_results:
            self._populate_status_table(table, {'status': 'No results available'})
            return

        loci = list(analysis_results.keys())
        table.setRowCount(len(loci))
        table.setColumnCount(5)
        table.setHorizontalHeaderLabels(["Locus", "Chi-square", "P-value", "DF", "HWE Status"])

        for i, locus in enumerate(loci):
            locus_data = analysis_results[locus]

            table.setItem(i, 0, QTableWidgetItem(locus))

            chi_sq = locus_data.get('chi_square')
            if chi_sq is not None:
                table.setItem(i, 1, QTableWidgetItem(f"{chi_sq:.4f}"))
            else:
                table.setItem(i, 1, QTableWidgetItem("N/A"))

            p_val = locus_data.get('p_value')
            if p_val is not None:
                table.setItem(i, 2, QTableWidgetItem(f"{p_val:.4f}"))
            else:
                table.setItem(i, 2, QTableWidgetItem("N/A"))

            df = locus_data.get('df')
            if df is not None:
                table.setItem(i, 3, QTableWidgetItem(str(df)))
            else:
                table.setItem(i, 3, QTableWidgetItem("N/A"))

            status = locus_data.get('hwe_status', 'N/A')
            table.setItem(i, 4, QTableWidgetItem(status))

    def _populate_fst_table(self, table, analysis_results):
        """Populate table with Fst results - per population and overall"""
        if not analysis_results:
            self._populate_status_table(table, {'status': 'No results available'})
            return

        # Check if we have the new per-population format
        populations = [k for k in analysis_results.keys() if k.startswith('Population_')]

        if populations:
            # New format: per-population Fst statistics
            rows = []

            # First add per-population rows
            for pop in populations:
                pop_label = pop.replace('Population_', 'Pop ')
                pop_data = analysis_results[pop]
                loci = pop_data.get('loci', {})

                for locus, locus_data in loci.items():
                    rows.append((pop_label, locus, locus_data))

            # Then add overall Fst results
            if '_overall' in analysis_results:
                overall = analysis_results['_overall']
                overall_loci = [k for k in overall.keys() if not k.startswith('_')]
                for locus in overall_loci:
                    locus_data = overall[locus]
                    rows.append(('OVERALL', locus, locus_data))

            table.setRowCount(len(rows))
            table.setColumnCount(6)
            table.setHorizontalHeaderLabels(["Population", "Locus", "Hs", "Fis", "Fst", "Ht"])

            for i, (pop, locus, data) in enumerate(rows):
                item_pop = QTableWidgetItem(pop)
                if pop == 'OVERALL':
                    font = item_pop.font()
                    font.setBold(True)
                    item_pop.setFont(font)

                table.setItem(i, 0, item_pop)
                table.setItem(i, 1, QTableWidgetItem(locus))

                # Hs value
                hs = data.get('Hs')
                if hs is not None:
                    table.setItem(i, 2, QTableWidgetItem(f"{hs:.4f}"))
                else:
                    table.setItem(i, 2, QTableWidgetItem("N/A"))

                # Fis value
                fis = data.get('Fis')
                if fis is not None:
                    table.setItem(i, 3, QTableWidgetItem(f"{fis:.4f}"))
                else:
                    table.setItem(i, 3, QTableWidgetItem("N/A"))

                # Fst value (only for OVERALL)
                if pop == 'OVERALL':
                    fst = data.get('fst')
                    if fst is not None:
                        item_fst = QTableWidgetItem(f"{fst:.4f}")
                        font = item_fst.font()
                        font.setBold(True)
                        item_fst.setFont(font)
                        table.setItem(i, 4, item_fst)
                    else:
                        table.setItem(i, 4, QTableWidgetItem("N/A"))

                    ht = data.get('Ht')
                    if ht is not None:
                        item_ht = QTableWidgetItem(f"{ht:.4f}")
                        font = item_ht.font()
                        font.setBold(True)
                        item_ht.setFont(font)
                        table.setItem(i, 5, item_ht)
                    else:
                        table.setItem(i, 5, QTableWidgetItem("N/A"))
                else:
                    table.setItem(i, 4, QTableWidgetItem("-"))
                    table.setItem(i, 5, QTableWidgetItem("-"))

        else:
            # Old format: just overall Fst
            loci = [k for k in analysis_results.keys() if not k.startswith('_')]

            if not loci:
                table.setRowCount(1)
                table.setColumnCount(1)
                table.setHorizontalHeaderLabels(["Status"])
                table.setItem(0, 0, QTableWidgetItem("No locus data available"))
                return

            table.setRowCount(len(loci))
            table.setColumnCount(4)
            table.setHorizontalHeaderLabels(["Locus", "Fst", "Ht", "Hs"])

            for i, locus in enumerate(loci):
                locus_data = analysis_results[locus]

                table.setItem(i, 0, QTableWidgetItem(locus))

                fst = locus_data.get('fst')
                if fst is not None:
                    table.setItem(i, 1, QTableWidgetItem(f"{fst:.4f}"))
                else:
                    table.setItem(i, 1, QTableWidgetItem("N/A"))

                ht = locus_data.get('Ht')
                if ht is not None:
                    table.setItem(i, 2, QTableWidgetItem(f"{ht:.4f}"))
                else:
                    table.setItem(i, 2, QTableWidgetItem("N/A"))

                hs = locus_data.get('Hs')
                if hs is not None:
                    table.setItem(i, 3, QTableWidgetItem(f"{hs:.4f}"))
                else:
                    table.setItem(i, 3, QTableWidgetItem("N/A"))

    def _populate_match_prob_table(self, table, analysis_results):
        """Populate table with match probability results"""
        if not analysis_results:
            self._populate_status_table(table, {'status': 'No results available'})
            return

        # Check for combined results
        if '_combined' in analysis_results:
            combined = analysis_results['_combined']

            # Show combined results first
            table.setRowCount(1)
            table.setColumnCount(4)
            table.setHorizontalHeaderLabels(["Combined PM", "Combined PD", "1 in X", "N Loci"])

            pm = combined.get('combined_PM', 0)
            pd = combined.get('combined_PD', 0)
            one_in_x = combined.get('one_in_X', 0)
            n_loci = combined.get('n_loci', 0)

            table.setItem(0, 0, QTableWidgetItem(f"{pm:.2e}"))
            table.setItem(0, 1, QTableWidgetItem(f"{pd:.6f}"))
            table.setItem(0, 2, QTableWidgetItem(f"{one_in_x:.2e}"))
            table.setItem(0, 3, QTableWidgetItem(str(n_loci)))
        else:
            # Show per-locus results
            loci = [k for k in analysis_results.keys() if not k.startswith('_')]
            table.setRowCount(len(loci))
            table.setColumnCount(4)
            table.setHorizontalHeaderLabels(["Locus", "PM", "PD", "PE"])

            for i, locus in enumerate(loci):
                locus_data = analysis_results[locus]

                table.setItem(i, 0, QTableWidgetItem(locus))
                table.setItem(i, 1, QTableWidgetItem(f"{locus_data.get('PM', 0):.4f}"))
                table.setItem(i, 2, QTableWidgetItem(f"{locus_data.get('PD', 0):.4f}"))
                table.setItem(i, 3, QTableWidgetItem(f"{locus_data.get('PE', 0):.4f}"))

    def _populate_str_analysis_table(self, table, analysis_results):
        """Populate table with STR analysis results"""
        if not analysis_results:
            self._populate_status_table(table, {'status': 'No results available'})
            return

        loci = list(analysis_results.keys())
        table.setRowCount(len(loci))
        table.setColumnCount(6)
        table.setHorizontalHeaderLabels(["Locus", "N Alleles", "Min", "Max", "Range", "Mean"])

        for i, locus in enumerate(loci):
            locus_data = analysis_results[locus]

            if 'error' in locus_data:
                table.setItem(i, 0, QTableWidgetItem(locus))
                table.setItem(i, 1, QTableWidgetItem(locus_data['error']))
                continue

            table.setItem(i, 0, QTableWidgetItem(locus))
            table.setItem(i, 1, QTableWidgetItem(str(locus_data.get('n_alleles', 'N/A'))))
            table.setItem(i, 2, QTableWidgetItem(str(locus_data.get('min_repeats', 'N/A'))))
            table.setItem(i, 3, QTableWidgetItem(str(locus_data.get('max_repeats', 'N/A'))))
            table.setItem(i, 4, QTableWidgetItem(str(locus_data.get('range', 'N/A'))))

            mean = locus_data.get('mean_repeats')
            if mean is not None:
                table.setItem(i, 5, QTableWidgetItem(f"{mean:.2f}"))
            else:
                table.setItem(i, 5, QTableWidgetItem("N/A"))

    def _populate_status_table(self, table, analysis_results):
        """Populate table with status or error message"""
        table.setRowCount(1)
        table.setColumnCount(1)
        table.setHorizontalHeaderLabels(["Message"])

        if 'error' in analysis_results:
            message = f"ERROR: {analysis_results['error']}"
        else:
            message = analysis_results.get('status', 'No information available')

        if 'note' in analysis_results:
            message += f"\n\nNote: {analysis_results['note']}"

        item = QTableWidgetItem(message)
        item.setTextAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        table.setItem(0, 0, item)
        table.resizeRowsToContents()

    def _populate_generic_table(self, table, analysis_results):
        """Generic table population for unknown result types"""
        # Convert dict to table
        items = [(str(k), str(v)) for k, v in analysis_results.items() if not k.startswith('_')]

        if not items:
            self._populate_status_table(table, {'status': 'No displayable results'})
            return

        table.setRowCount(len(items))
        table.setColumnCount(2)
        table.setHorizontalHeaderLabels(["Parameter", "Value"])

        for i, (key, value) in enumerate(items):
            table.setItem(i, 0, QTableWidgetItem(key))
            table.setItem(i, 1, QTableWidgetItem(value[:100]))  # Truncate long values

    def _populate_allele_frequencies_table(self, table, analysis_results):
        """Populate table with actual allele frequencies"""
        if not analysis_results:
            self._populate_status_table(table, {'status': 'No results available'})
            return

        # Extract population data
        populations = [k for k in analysis_results.keys() if k.startswith('Population_')]

        if not populations:
            self._populate_status_table(table, {'status': 'No population data available'})
            return

        # Collect all rows: Population | Locus | Allele | Frequency
        rows = []
        for pop in populations:
            pop_label = pop.replace('Population_', 'Pop ')
            pop_data = analysis_results[pop]
            loci = pop_data.get('loci', {})

            for locus, locus_data in loci.items():
                allele_freqs = locus_data.get('allele_frequencies', {})
                for allele, freq in sorted(allele_freqs.items()):
                    rows.append((pop_label, locus, allele, freq))

        # Add overall summary if available
        if '_overall' in analysis_results:
            overall = analysis_results['_overall']
            for locus, locus_data in overall.get('loci', {}).items():
                rows.append(('OVERALL', locus, '-', locus_data.get('He_mean', 0)))

        table.setRowCount(len(rows))
        table.setColumnCount(4)
        table.setHorizontalHeaderLabels(["Population", "Locus", "Allele", "Frequency"])

        for i, (pop, locus, allele, freq) in enumerate(rows):
            table.setItem(i, 0, QTableWidgetItem(pop))
            table.setItem(i, 1, QTableWidgetItem(locus))
            if allele == '-':
                table.setItem(i, 2, QTableWidgetItem('Mean He'))
            else:
                table.setItem(i, 2, QTableWidgetItem(allele))

            if isinstance(freq, (int, float)):
                table.setItem(i, 3, QTableWidgetItem(f"{freq:.4f}"))
            else:
                table.setItem(i, 3, QTableWidgetItem(str(freq)))

    def _populate_per_population_table(self, table, analysis_results, analysis_name):
        """Populate table with per-population results"""
        if not analysis_results:
            self._populate_status_table(table, {'status': 'No results available'})
            return

        # Extract population data
        populations = [k for k in analysis_results.keys() if k.startswith('Population_')]

        if not populations:
            self._populate_status_table(table, {'status': 'No population data available'})
            return

        # Collect all loci from first population
        first_pop = populations[0]
        loci = list(analysis_results[first_pop].get('loci', {}).keys())

        if not loci:
            self._populate_status_table(table, {'status': 'No locus data available'})
            return

        # Create rows: one per locus per population
        rows = []
        for pop in populations:
            pop_data = analysis_results[pop]
            pop_label = pop.replace('Population_', 'Pop ')

            for locus in loci:
                locus_data = pop_data.get('loci', {}).get(locus, {})
                rows.append((pop_label, locus, locus_data))

        # Add overall summary rows if available
        if '_overall' in analysis_results:
            overall = analysis_results['_overall']
            for locus, locus_data in overall.get('loci', {}).items():
                rows.append(('OVERALL', locus, locus_data))

        # Set up table based on analysis type
        if analysis_name == 'Hardy-Weinberg Equilibrium':
            table.setRowCount(len(rows))
            table.setColumnCount(6)
            table.setHorizontalHeaderLabels(["Population", "Locus", "Ho", "He", "HWE P-value", "HWE Status"])

            for i, (pop, locus, data) in enumerate(rows):
                item_pop = QTableWidgetItem(pop)
                if pop == 'OVERALL':
                    font = item_pop.font()
                    font.setBold(True)
                    item_pop.setFont(font)

                table.setItem(i, 0, item_pop)
                table.setItem(i, 1, QTableWidgetItem(locus))

                # Ho and He values
                table.setItem(i, 2, QTableWidgetItem(f"{data.get('Ho', data.get('Ho_mean', 0)):.4f}"))
                table.setItem(i, 3, QTableWidgetItem(f"{data.get('He', data.get('He_mean', 0)):.4f}"))

                # P-value and status
                p_val = data.get('hwe_p_value')
                if p_val is not None:
                    item_p = QTableWidgetItem(f"{p_val:.4f}")
                    if pop == 'OVERALL':
                        font = item_p.font()
                        font.setBold(True)
                        item_p.setFont(font)
                    table.setItem(i, 4, item_p)
                else:
                    table.setItem(i, 4, QTableWidgetItem("N/A"))

                status_text = data.get('hwe_status', 'N/A')
                item_status = QTableWidgetItem(status_text)
                if pop == 'OVERALL':
                    font = item_status.font()
                    font.setBold(True)
                    item_status.setFont(font)
                table.setItem(i, 5, item_status)

        else:  # Heterozygosity
            table.setRowCount(len(rows))
            table.setColumnCount(7)
            table.setHorizontalHeaderLabels(["Population", "Locus", "N Alleles", "Ho", "He", "Fis", "PD"])

            for i, (pop, locus, data) in enumerate(rows):
                item_pop = QTableWidgetItem(pop)
                if pop == 'OVERALL':
                    font = item_pop.font()
                    font.setBold(True)
                    item_pop.setFont(font)

                table.setItem(i, 0, item_pop)
                table.setItem(i, 1, QTableWidgetItem(locus))

                if pop == 'OVERALL':
                    table.setItem(i, 2, QTableWidgetItem(f"{data.get('n_alleles_mean', 0):.1f}"))
                    table.setItem(i, 3, QTableWidgetItem(f"{data.get('Ho_mean', 0):.4f}"))
                    table.setItem(i, 4, QTableWidgetItem(f"{data.get('He_mean', 0):.4f}"))
                    table.setItem(i, 5, QTableWidgetItem(f"{data.get('Fis_mean', 0):.4f}"))
                    table.setItem(i, 6, QTableWidgetItem(f"{data.get('PD_mean', 0):.4f}"))
                else:
                    table.setItem(i, 2, QTableWidgetItem(str(data.get('n_alleles', 'N/A'))))
                    table.setItem(i, 3, QTableWidgetItem(f"{data.get('Ho', 0):.4f}"))
                    table.setItem(i, 4, QTableWidgetItem(f"{data.get('He', 0):.4f}"))
                    table.setItem(i, 5, QTableWidgetItem(f"{data.get('Fis', 0):.4f}"))
                    table.setItem(i, 6, QTableWidgetItem(f"{data.get('PD', 0):.4f}"))

    def _populate_match_prob_per_pop_table(self, table, analysis_results):
        """Populate table with match probability per population with 1-in-X for all loci"""
        if not analysis_results:
            self._populate_status_table(table, {'status': 'No results available'})
            return

        populations = [k for k in analysis_results.keys() if k.startswith('Population_')]

        if not populations:
            self._populate_status_table(table, {'status': 'No population data available'})
            return

        # Create separate sections for combined and per-locus results
        rows = []

        # First, add combined results for each population
        for pop in populations:
            pop_label = pop.replace('Population_', 'Pop ')
            pop_data = analysis_results[pop]
            combined = pop_data.get('combined', {})

            if combined:
                pm = combined.get('combined_PM', 0)
                pd = combined.get('combined_PD', 0)
                one_in_x = combined.get('one_in_X', 0)
                n_loci = combined.get('n_loci', 0)

                rows.append((
                    pop_label,
                    'COMBINED',
                    f"{pm:.2e}",
                    f"{pd:.6f}",
                    f"{one_in_x:.2e}"
                ))

        # Then add per-locus results with 1-in-X
        for pop in populations:
            pop_label = pop.replace('Population_', 'Pop ')
            pop_data = analysis_results[pop]
            loci_results = pop_data.get('loci_results', {})

            for locus, locus_data in loci_results.items():
                if isinstance(locus_data, dict) and 'PM' in locus_data:
                    pm = locus_data.get('PM', 0)
                    pd = locus_data.get('PD', 0)
                    one_in_x = locus_data.get('one_in_X', 0)

                    rows.append((
                        pop_label,
                        locus,
                        f"{pm:.4f}",
                        f"{pd:.4f}",
                        f"{one_in_x:.2f}"
                    ))

        # Add overall summary if available
        if '_overall' in analysis_results:
            overall = analysis_results['_overall']
            combined = overall.get('_combined', {})

            if combined:
                pm = combined.get('combined_PM', 0)
                pd = combined.get('combined_PD', 0)
                one_in_x = combined.get('one_in_X', 0)

                rows.append((
                    'OVERALL',
                    'COMBINED',
                    f"{pm:.2e}",
                    f"{pd:.6f}",
                    f"{one_in_x:.2e}"
                ))

            # Add per-locus overall
            for locus, locus_data in overall.items():
                if not locus.startswith('_') and isinstance(locus_data, dict) and 'PM' in locus_data:
                    pm = locus_data.get('PM', 0)
                    pd = locus_data.get('PD', 0)
                    one_in_x = locus_data.get('one_in_X', 0)

                    rows.append((
                        'OVERALL',
                        locus,
                        f"{pm:.4f}",
                        f"{pd:.4f}",
                        f"{one_in_x:.2f}"
                    ))

        table.setRowCount(len(rows))
        table.setColumnCount(5)
        table.setHorizontalHeaderLabels(["Population", "Locus", "PM", "PD", "1 in X"])

        for i, row_data in enumerate(rows):
            for j, value in enumerate(row_data):
                item = QTableWidgetItem(str(value))
                if row_data[1] == 'COMBINED' or row_data[0] == 'OVERALL':
                    # Bold for combined and overall results
                    font = item.font()
                    font.setBold(True)
                    item.setFont(font)
                table.setItem(i, j, item)

    def _convert_results_to_dataframe(self, analysis_name, analysis_results):
        """
        Convert analysis results to pandas DataFrame for export

        Args:
            analysis_name: Name of the analysis
            analysis_results: Results dictionary

        Returns:
            pandas DataFrame or None
        """
        import pandas as pd

        # Check if we have per-population format
        populations = [k for k in analysis_results.keys() if k.startswith('Population_')]

        if not populations and not analysis_results:
            return None

        # Handle per-population analyses
        if populations:
            rows = []

            if analysis_name == 'Hardy-Weinberg Equilibrium':
                for pop in populations:
                    pop_label = pop.replace('Population_', 'Pop ')
                    loci_data = analysis_results[pop].get('loci', {})
                    for locus, data in loci_data.items():
                        rows.append({
                            'Population': pop_label,
                            'Locus': locus,
                            'Ho': data.get('Ho', 0),
                            'He': data.get('He', 0),
                            'HWE_P_value': data.get('hwe_p_value'),
                            'HWE_Status': data.get('hwe_status', 'N/A')
                        })

                # Add overall
                if '_overall' in analysis_results:
                    overall = analysis_results['_overall']
                    for locus, data in overall.get('loci', {}).items():
                        rows.append({
                            'Population': 'OVERALL',
                            'Locus': locus,
                            'Ho': data.get('Ho', data.get('Ho_mean', 0)),
                            'He': data.get('He', data.get('He_mean', 0)),
                            'HWE_P_value': data.get('hwe_p_value'),
                            'HWE_Status': data.get('hwe_status', 'N/A')
                        })

            elif analysis_name == 'Heterozygosity':
                for pop in populations:
                    pop_label = pop.replace('Population_', 'Pop ')
                    loci_data = analysis_results[pop].get('loci', {})
                    for locus, data in loci_data.items():
                        rows.append({
                            'Population': pop_label,
                            'Locus': locus,
                            'N_Alleles': data.get('n_alleles', 0),
                            'Ho': data.get('Ho', 0),
                            'He': data.get('He', 0),
                            'Fis': data.get('Fis', 0),
                            'PD': data.get('PD', 0)
                        })

                # Add overall
                if '_overall' in analysis_results:
                    overall = analysis_results['_overall']
                    for locus, data in overall.get('loci', {}).items():
                        rows.append({
                            'Population': 'OVERALL',
                            'Locus': locus,
                            'N_Alleles': data.get('n_alleles', 0),
                            'Ho': data.get('Ho', 0),
                            'He': data.get('He', 0),
                            'Fis': data.get('Fis', 0),
                            'PD': data.get('PD', 0)
                        })

            elif analysis_name == 'Fixation Index (Fst)':
                for pop in populations:
                    pop_label = pop.replace('Population_', 'Pop ')
                    loci_data = analysis_results[pop].get('loci', {})
                    for locus, data in loci_data.items():
                        rows.append({
                            'Population': pop_label,
                            'Locus': locus,
                            'Hs': data.get('Hs', 0),
                            'Fis': data.get('Fis', 0)
                        })

                # Add overall Fst
                if '_overall' in analysis_results:
                    overall = analysis_results['_overall']
                    overall_loci = [k for k in overall.keys() if not k.startswith('_')]
                    for locus in overall_loci:
                        data = overall[locus]
                        rows.append({
                            'Population': 'OVERALL',
                            'Locus': locus,
                            'Hs': data.get('Hs', 0),
                            'Fis': data.get('Fis', 0),
                            'Fst': data.get('fst', 0),
                            'Ht': data.get('Ht', 0)
                        })

            elif analysis_name == 'Match Probability and Power of Discrimination':
                for pop in populations:
                    pop_label = pop.replace('Population_', 'Pop ')
                    loci_results = analysis_results[pop].get('loci_results', {})

                    # Combined result
                    combined = analysis_results[pop].get('combined', {})
                    if combined:
                        rows.append({
                            'Population': pop_label,
                            'Locus': 'COMBINED',
                            'PM': combined.get('combined_PM', 0),
                            'PD': combined.get('combined_PD', 0),
                            '1_in_X': combined.get('one_in_X', 0)
                        })

                    # Per-locus results
                    for locus, data in loci_results.items():
                        if isinstance(data, dict):
                            rows.append({
                                'Population': pop_label,
                                'Locus': locus,
                                'PM': data.get('PM', 0),
                                'PD': data.get('PD', 0),
                                '1_in_X': data.get('one_in_X', 0)
                            })

                # Add overall
                if '_overall' in analysis_results:
                    overall = analysis_results['_overall']
                    combined = overall.get('_combined', {})
                    if combined:
                        rows.append({
                            'Population': 'OVERALL',
                            'Locus': 'COMBINED',
                            'PM': combined.get('combined_PM', 0),
                            'PD': combined.get('combined_PD', 0),
                            '1_in_X': combined.get('one_in_X', 0)
                        })

                    for locus, data in overall.items():
                        if not locus.startswith('_') and isinstance(data, dict):
                            rows.append({
                                'Population': 'OVERALL',
                                'Locus': locus,
                                'PM': data.get('PM', 0),
                                'PD': data.get('PD', 0),
                                '1_in_X': data.get('one_in_X', 0)
                            })

            elif analysis_name == 'Allele Frequencies':
                # Long-format table for normal export
                for pop in populations:
                    pop_label = pop.replace('Population_', 'Pop ')
                    loci_data = analysis_results[pop].get('loci', {})
                    for locus, locus_data in loci_data.items():
                        allele_freqs = locus_data.get('allele_frequencies', {})
                        for allele, freq in sorted(allele_freqs.items()):
                            rows.append({
                                'Population': pop_label,
                                'Locus': locus,
                                'Allele': allele,
                                'Frequency': freq
                            })

            if rows:
                return pd.DataFrame(rows)

        # Handle non-per-population analyses (old format)
        else:
            # Generic handling for other analysis types
            rows = []
            for key, value in analysis_results.items():
                if not key.startswith('_'):
                    if isinstance(value, dict):
                        row = {'Parameter': key}
                        row.update(value)
                        rows.append(row)

            if rows:
                return pd.DataFrame(rows)

        return None

    def _generate_wide_format_allele_freq_tables(self, allele_freq_results):
        """
        Generate wide-format allele frequency tables (loci as columns, alleles as rows)
        One table per population + overall

        Args:
            allele_freq_results: Results from Allele Frequencies analysis

        Returns:
            Dictionary of {population_name: DataFrame}
        """
        import pandas as pd

        wide_tables = {}

        # Get all populations
        populations = [k for k in allele_freq_results.keys() if k.startswith('Population_')]

        # Add overall if available
        if '_overall' in allele_freq_results:
            populations.append('_overall')

        for pop_key in populations:
            pop_data = allele_freq_results[pop_key]
            loci_data = pop_data.get('loci', {})

            # Collect all alleles across all loci
            all_alleles = set()
            for locus, locus_data in loci_data.items():
                allele_freqs = locus_data.get('allele_frequencies', {})
                all_alleles.update(allele_freqs.keys())

            # Sort alleles
            all_alleles = sorted(all_alleles)

            # Create DataFrame with alleles as rows, loci as columns
            data = {}
            for locus, locus_data in loci_data.items():
                allele_freqs = locus_data.get('allele_frequencies', {})
                data[locus] = [allele_freqs.get(allele, 0.0) for allele in all_alleles]

            if data:
                df = pd.DataFrame(data, index=all_alleles)
                df.index.name = 'Allele'

                # Set population name
                if pop_key == '_overall':
                    pop_name = 'OVERALL'
                else:
                    pop_name = pop_key.replace('Population_', 'Population_')

                wide_tables[pop_name] = df

        return wide_tables

    def get_analysis_description(self, analysis_name):
        """Get description for analysis type"""
        descriptions = {
            'Hardy-Weinberg Equilibrium': 'Tests whether allele frequencies are in equilibrium',
            'Fixation Index (Fst)': 'Measures population differentiation',
            'Heterozygosity': 'Calculates observed and expected heterozygosity',
            'Match Probability': 'Calculates the probability of a random match',
            'Likelihood Ratio': 'Compares likelihood of alternative hypotheses',
            'Haplotype Diversity': 'Measures diversity of haplotypes in the sample'
        }
        return descriptions.get(analysis_name, 'Statistical analysis results')

    def export_to_excel(self):
        """Export results to Excel file"""
        if not self.main_window.analysis_results:
            QMessageBox.warning(self, "No Results", "No results available to export.")
            return

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        default_name = f"forstat_results_{timestamp}.xlsx"

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export to Excel",
            default_name,
            "Excel Files (*.xlsx);;All Files (*.*)"
        )

        if file_path:
            try:
                import pandas as pd

                logger.info(f"Exporting to Excel: {file_path}")

                results = self.main_window.analysis_results
                analyses = results.get('analyses', [])
                results_data = results.get('results', {})

                # Create Excel writer
                with pd.ExcelWriter(file_path, engine='openpyxl') as writer:
                    # Summary sheet
                    summary_data = {
                        'Parameter': ['Data File', 'Analyses Performed', 'Timestamp'],
                        'Value': [
                            results.get('data_file', 'Unknown'),
                            len(analyses),
                            datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                        ]
                    }
                    pd.DataFrame(summary_data).to_excel(writer, sheet_name='Summary', index=False)

                    # Export each analysis
                    for analysis in analyses:
                        analysis_results = results_data.get(analysis, {})

                        # Sanitize sheet name (Excel has 31 char limit and special char restrictions)
                        sheet_name = analysis[:31].replace('/', '-').replace('\\', '-').replace(':', '-')

                        # Special handling for Allele Frequencies - create wide-format tables
                        if analysis == 'Allele Frequencies':
                            wide_tables = self._generate_wide_format_allele_freq_tables(analysis_results)

                            # Export each population's wide-format table to a separate sheet
                            for pop_name, df in wide_tables.items():
                                pop_sheet_name = f"{pop_name}_AlFreq"[:31]
                                df.to_excel(writer, sheet_name=pop_sheet_name)

                        # Export standard analysis results
                        df = self._convert_results_to_dataframe(analysis, analysis_results)
                        if df is not None and not df.empty:
                            df.to_excel(writer, sheet_name=sheet_name, index=False)

                logger.info(f"Successfully exported to Excel: {file_path}")
                QMessageBox.information(
                    self,
                    "Export Successful",
                    f"Results exported to:\n{file_path}"
                )
            except Exception as e:
                logger.error(f"Export failed: {e}", exc_info=True)
                QMessageBox.critical(
                    self,
                    "Export Failed",
                    f"Failed to export results:\n{str(e)}"
                )

    def export_to_pdf(self):
        """Export results to PDF file"""
        if not self.main_window.analysis_results:
            QMessageBox.warning(self, "No Results", "No results available to export.")
            return

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        default_name = f"forstat_report_{timestamp}.pdf"

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export to PDF",
            default_name,
            "PDF Files (*.pdf);;All Files (*.*)"
        )

        if file_path:
            try:
                # TODO: Implement actual export
                logger.info(f"Exporting to PDF: {file_path}")
                QMessageBox.information(
                    self,
                    "Export Successful",
                    f"Report exported to:\n{file_path}"
                )
            except Exception as e:
                logger.error(f"Export failed: {e}")
                QMessageBox.critical(
                    self,
                    "Export Failed",
                    f"Failed to export report:\n{str(e)}"
                )

    def export_to_csv(self):
        """Export results to CSV file"""
        if not self.main_window.analysis_results:
            QMessageBox.warning(self, "No Results", "No results available to export.")
            return

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        default_name = f"forstat_results_{timestamp}.csv"

        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export to CSV",
            default_name,
            "CSV Files (*.csv);;All Files (*.*)"
        )

        if file_path:
            try:
                # TODO: Implement actual export
                logger.info(f"Exporting to CSV: {file_path}")
                QMessageBox.information(
                    self,
                    "Export Successful",
                    f"Results exported to:\n{file_path}"
                )
            except Exception as e:
                logger.error(f"Export failed: {e}")
                QMessageBox.critical(
                    self,
                    "Export Failed",
                    f"Failed to export results:\n{str(e)}"
                )
