"""
Analysis page for selecting and running analyses
"""
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QGroupBox, QCheckBox, QProgressBar, QTextEdit, QScrollArea,
    QMessageBox
)
from PyQt6.QtCore import Qt, pyqtSignal, QThread
from forstat.utils.config import Config
from forstat.utils.logger import get_logger

logger = get_logger(__name__)


class AnalysisPage(QWidget):
    """Page for selecting and running analyses"""

    analysis_complete = pyqtSignal(object)  # Signal when analysis is complete

    def __init__(self, parent=None):
        super().__init__(parent)
        self.main_window = parent
        self.selected_analyses = []

        self.init_ui()

    def init_ui(self):
        """Initialize user interface"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(30, 30, 30, 30)
        layout.setSpacing(20)

        # Title
        title = QLabel("Analysis")
        title.setProperty("class", "title")
        layout.addWidget(title)

        # Instructions
        instructions = QLabel(
            "Select the analyses you want to perform on your data"
        )
        instructions.setProperty("class", "hint")
        layout.addWidget(instructions)

        # Scrollable area for analysis options
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QScrollArea.Shape.NoFrame)

        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)
        scroll_layout.setSpacing(15)

        # Population genetics analyses
        pop_group = self.create_analysis_group(
            "Population Genetics",
            Config.ANALYSIS_TYPES['population']
        )
        scroll_layout.addWidget(pop_group)

        # Forensic analyses
        forensic_group = self.create_analysis_group(
            "Forensic Statistics",
            Config.ANALYSIS_TYPES['forensic']
        )
        scroll_layout.addWidget(forensic_group)

        # STR analyses
        str_group = self.create_analysis_group(
            "STR Analysis",
            Config.ANALYSIS_TYPES['str']
        )
        scroll_layout.addWidget(str_group)

        # mtDNA analyses
        mtdna_group = self.create_analysis_group(
            "Mitochondrial DNA",
            Config.ANALYSIS_TYPES['mtdna']
        )
        scroll_layout.addWidget(mtdna_group)

        scroll_layout.addStretch()

        scroll_area.setWidget(scroll_widget)
        layout.addWidget(scroll_area, 1)

        # Progress section
        self.progress_group = QGroupBox("Progress")
        progress_layout = QVBoxLayout(self.progress_group)

        self.progress_bar = QProgressBar()
        self.progress_bar.setValue(0)
        progress_layout.addWidget(self.progress_bar)

        self.progress_text = QTextEdit()
        self.progress_text.setReadOnly(True)
        self.progress_text.setMaximumHeight(100)
        progress_layout.addWidget(self.progress_text)

        self.progress_group.hide()  # Hidden until analysis starts
        layout.addWidget(self.progress_group)

        # Action buttons
        button_layout = QHBoxLayout()

        self.select_all_btn = QPushButton("Select All")
        self.select_all_btn.setProperty("class", "secondary")
        self.select_all_btn.clicked.connect(self.select_all)
        button_layout.addWidget(self.select_all_btn)

        self.clear_all_btn = QPushButton("Clear All")
        self.clear_all_btn.setProperty("class", "secondary")
        self.clear_all_btn.clicked.connect(self.clear_all)
        button_layout.addWidget(self.clear_all_btn)

        button_layout.addStretch()

        self.run_btn = QPushButton("Run Analysis")
        self.run_btn.setProperty("class", "success")
        self.run_btn.clicked.connect(self.run_analysis)
        button_layout.addWidget(self.run_btn)

        layout.addLayout(button_layout)

    def create_analysis_group(self, title, analyses):
        """Create a group of analysis checkboxes"""
        group = QGroupBox(title)
        layout = QVBoxLayout(group)

        self.checkboxes = getattr(self, 'checkboxes', [])

        for analysis in analyses:
            checkbox = QCheckBox(analysis)
            checkbox.stateChanged.connect(self.on_selection_changed)
            self.checkboxes.append(checkbox)
            layout.addWidget(checkbox)

        return group

    def on_selection_changed(self):
        """Handle analysis selection change"""
        self.selected_analyses = [
            cb.text() for cb in self.checkboxes if cb.isChecked()
        ]
        logger.info(f"Selected analyses: {len(self.selected_analyses)}")

    def select_all(self):
        """Select all analyses"""
        for cb in self.checkboxes:
            cb.setChecked(True)

    def clear_all(self):
        """Clear all selections"""
        for cb in self.checkboxes:
            cb.setChecked(False)

    def run_analysis(self):
        """Run selected analyses"""
        if not self.selected_analyses:
            QMessageBox.warning(
                self,
                "No Analyses Selected",
                "Please select at least one analysis to run."
            )
            return

        if not self.main_window.current_data:
            QMessageBox.warning(
                self,
                "No Data Loaded",
                "Please load data first from the Upload page."
            )
            return

        # Show progress group
        self.progress_group.show()
        self.progress_bar.setValue(0)
        self.progress_text.clear()

        # Disable run button during analysis
        self.run_btn.setEnabled(False)

        self.add_progress_message("Starting analysis...")
        logger.info(f"Running {len(self.selected_analyses)} analyses")

        # Simulate analysis (will be replaced with actual analysis)
        self.perform_mock_analysis()

    def perform_mock_analysis(self):
        """Perform actual analysis"""
        try:
            from PyQt6.QtWidgets import QApplication
            from forstat.analysis.population.allele_frequencies import calculate_summary_statistics

            total = len(self.selected_analyses)
            results_data = {}

            # Get genetic data
            genetic_data = self.main_window.current_data

            for i, analysis in enumerate(self.selected_analyses, 1):
                self.add_progress_message(f"Running: {analysis}")

                # Run actual analyses
                if analysis in ['Allele Frequencies', 'Heterozygosity']:
                    self.add_progress_message("  Calculating allele frequencies and heterozygosity...")
                    stats = calculate_summary_statistics(genetic_data)
                    results_data[analysis] = stats

                elif analysis == 'Hardy-Weinberg Equilibrium':
                    self.add_progress_message("  Testing Hardy-Weinberg equilibrium...")
                    # Placeholder for HWE test
                    results_data[analysis] = {'status': 'Not yet implemented'}

                else:
                    # Placeholder for other analyses
                    results_data[analysis] = {'status': 'Not yet implemented'}

                progress = int((i / total) * 100)
                self.progress_bar.setValue(progress)

                # Process events to update UI
                QApplication.processEvents()

            self.add_progress_message("Analysis complete!")
            self.progress_bar.setValue(100)

            # Re-enable run button
            self.run_btn.setEnabled(True)

            # Create results dictionary
            results = {
                'analyses': self.selected_analyses,
                'data_file': genetic_data.file_name,
                'genetic_data': genetic_data,
                'results': results_data
            }

            # Emit completion signal
            self.analysis_complete.emit(results)

            QMessageBox.information(
                self,
                "Analysis Complete",
                f"Successfully completed {len(self.selected_analyses)} analyses!\n\n"
                "View results in the Results page."
            )

        except Exception as e:
            logger.error(f"Analysis failed: {e}", exc_info=True)
            self.run_btn.setEnabled(True)
            self.add_progress_message(f"ERROR: {str(e)}")
            QMessageBox.critical(
                self,
                "Analysis Failed",
                f"An error occurred during analysis:\n{str(e)}"
            )

    def add_progress_message(self, message):
        """Add message to progress text"""
        self.progress_text.append(f"• {message}")
        # Scroll to bottom
        cursor = self.progress_text.textCursor()
        cursor.movePosition(cursor.MoveOperation.End)
        self.progress_text.setTextCursor(cursor)
