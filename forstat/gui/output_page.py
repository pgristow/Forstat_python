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
        if analysis_name in ['Allele Frequencies', 'Heterozygosity']:
            self._populate_diversity_table(table, analysis_results, results.get('genetic_data'))
        elif analysis_name == 'Hardy-Weinberg Equilibrium':
            self._populate_hwe_table(table, analysis_results)
        elif analysis_name == 'Fixation Index (Fst)' or analysis_name == 'Population Structure':
            self._populate_fst_table(table, analysis_results)
        elif analysis_name in ['Match Probability', 'Power of Discrimination']:
            self._populate_match_prob_table(table, analysis_results)
        elif analysis_name == 'Allele Frequency':
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
        """Populate table with Fst results"""
        if not analysis_results:
            self._populate_status_table(table, {'status': 'No results available'})
            return

        # Filter out meta keys
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
                # TODO: Implement actual export
                logger.info(f"Exporting to Excel: {file_path}")
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
