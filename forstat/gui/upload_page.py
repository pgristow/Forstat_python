"""
Upload page for data import
"""
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QTableWidget, QTableWidgetItem, QFileDialog, QGroupBox,
    QComboBox, QTextEdit, QMessageBox
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QDragEnterEvent, QDropEvent

from pathlib import Path
from forstat.utils.logger import get_logger
from forstat.utils.helpers import format_file_size, get_file_extension

logger = get_logger(__name__)


class UploadPage(QWidget):
    """Page for uploading and previewing data files"""

    data_loaded = pyqtSignal(object)  # Signal when data is successfully loaded

    def __init__(self, parent=None):
        super().__init__(parent)
        self.main_window = parent
        self.current_file = None
        self.current_data = None

        self.init_ui()

    def init_ui(self):
        """Initialize user interface"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(30, 30, 30, 30)
        layout.setSpacing(20)

        # Title
        title = QLabel("Upload Data")
        title.setProperty("class", "title")
        layout.addWidget(title)

        # Instructions
        instructions = QLabel(
            "Import your genetic data files in GenePop, Excel, CSV, or FASTA format"
        )
        instructions.setProperty("class", "hint")
        layout.addWidget(instructions)

        # File upload section
        upload_group = self.create_upload_section()
        layout.addWidget(upload_group)

        # File info section
        self.info_group = self.create_info_section()
        self.info_group.hide()  # Hidden until file is loaded
        layout.addWidget(self.info_group)

        # Data preview section
        self.preview_group = self.create_preview_section()
        self.preview_group.hide()  # Hidden until file is loaded
        layout.addWidget(self.preview_group, 1)

        # Action buttons
        button_layout = QHBoxLayout()
        button_layout.addStretch()

        self.clear_btn = QPushButton("Clear")
        self.clear_btn.setProperty("class", "secondary")
        self.clear_btn.clicked.connect(self.clear_data)
        self.clear_btn.setEnabled(False)
        button_layout.addWidget(self.clear_btn)

        self.load_btn = QPushButton("Load Data")
        self.load_btn.setProperty("class", "success")
        self.load_btn.clicked.connect(self.load_data)
        self.load_btn.setEnabled(False)
        button_layout.addWidget(self.load_btn)

        layout.addLayout(button_layout)

    def create_upload_section(self):
        """Create file upload section"""
        group = QGroupBox("Select File")
        layout = QVBoxLayout(group)

        # Drag and drop area
        self.drop_area = QLabel("Drag and drop file here\nor click 'Browse' to select")
        self.drop_area.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.drop_area.setMinimumHeight(150)
        self.drop_area.setStyleSheet("""
            QLabel {
                border: 3px dashed #2196F3;
                border-radius: 8px;
                background-color: #E3F2FD;
                color: #1976D2;
                font-size: 12pt;
            }
        """)
        self.drop_area.setAcceptDrops(True)
        self.drop_area.dragEnterEvent = self.drag_enter_event
        self.drop_area.dropEvent = self.drop_event
        layout.addWidget(self.drop_area)

        # Browse button and file type selector
        controls_layout = QHBoxLayout()

        self.file_type_combo = QComboBox()
        self.file_type_combo.addItems([
            "Auto-detect",
            "GenePop (.gen, .txt)",
            "Excel (.xlsx, .xls)",
            "CSV (.csv)",
            "FASTA (.fasta, .fa)"
        ])
        controls_layout.addWidget(QLabel("File Type:"))
        controls_layout.addWidget(self.file_type_combo)

        controls_layout.addStretch()

        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self.browse_file)
        controls_layout.addWidget(browse_btn)

        layout.addLayout(controls_layout)

        return group

    def create_info_section(self):
        """Create file information section"""
        group = QGroupBox("File Information")
        layout = QVBoxLayout(group)

        self.info_text = QTextEdit()
        self.info_text.setReadOnly(True)
        self.info_text.setMaximumHeight(100)
        layout.addWidget(self.info_text)

        return group

    def create_preview_section(self):
        """Create data preview section"""
        group = QGroupBox("Data Preview")
        layout = QVBoxLayout(group)

        self.preview_table = QTableWidget()
        self.preview_table.setAlternatingRowColors(True)
        layout.addWidget(self.preview_table)

        preview_hint = QLabel("Showing first 20 rows")
        preview_hint.setProperty("class", "hint")
        layout.addWidget(preview_hint)

        return group

    def drag_enter_event(self, event: QDragEnterEvent):
        """Handle drag enter event"""
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def drop_event(self, event: QDropEvent):
        """Handle drop event"""
        files = [u.toLocalFile() for u in event.mimeData().urls()]
        if files:
            self.load_file(files[0])

    def browse_file(self):
        """Open file browser dialog"""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Data File",
            "",
            "All Supported Files (*.gen *.txt *.xlsx *.xls *.csv *.fasta *.fa);;"
            "GenePop Files (*.gen *.txt);;"
            "Excel Files (*.xlsx *.xls);;"
            "CSV Files (*.csv);;"
            "FASTA Files (*.fasta *.fa);;"
            "All Files (*.*)"
        )

        if file_path:
            self.load_file(file_path)

    def load_file(self, file_path):
        """Load and preview file"""
        try:
            self.current_file = Path(file_path)
            logger.info(f"Loading file: {file_path}")

            # Update drop area
            self.drop_area.setText(f"Selected:\n{self.current_file.name}")

            # Show file information
            self.show_file_info()

            # Preview file content
            self.preview_file()

            # Enable buttons
            self.clear_btn.setEnabled(True)
            self.load_btn.setEnabled(True)

            self.info_group.show()
            self.preview_group.show()

        except Exception as e:
            logger.error(f"Error loading file: {e}")
            QMessageBox.critical(
                self,
                "Error",
                f"Failed to load file:\n{str(e)}"
            )

    def show_file_info(self):
        """Display file information"""
        if not self.current_file:
            return

        file_size = self.current_file.stat().st_size
        file_ext = get_file_extension(str(self.current_file))

        info_text = f"""
        <b>File Name:</b> {self.current_file.name}<br>
        <b>File Size:</b> {format_file_size(file_size)}<br>
        <b>File Type:</b> {file_ext}<br>
        <b>Location:</b> {self.current_file.parent}
        """

        self.info_text.setHtml(info_text)

    def preview_file(self):
        """Preview file content in table"""
        if not self.current_file:
            return

        try:
            # Read first few lines for preview
            with open(self.current_file, 'r', encoding='utf-8', errors='ignore') as f:
                lines = [f.readline() for _ in range(21)]  # Read 21 lines (header + 20 rows)

            if not lines:
                return

            # Simple preview - just show raw lines
            self.preview_table.clear()
            self.preview_table.setRowCount(min(20, len(lines)))
            self.preview_table.setColumnCount(1)
            self.preview_table.setHorizontalHeaderLabels(["Content"])

            for i, line in enumerate(lines[:20]):
                item = QTableWidgetItem(line.strip())
                self.preview_table.setItem(i, 0, item)

            self.preview_table.resizeColumnsToContents()

        except Exception as e:
            logger.error(f"Error previewing file: {e}")

    def load_data(self):
        """Load data and emit signal"""
        if not self.current_file:
            return

        try:
            # Import parsers
            from forstat.data.parsers.genepop import parse_genepop
            from forstat.data.models import GeneticData

            # Determine file type
            file_ext = get_file_extension(str(self.current_file))
            file_type = self.file_type_combo.currentText()

            # Parse based on file type
            if file_ext in ['.gen', '.txt'] or 'GenePop' in file_type:
                logger.info("Parsing as GenePop format")
                parsed_data = parse_genepop(str(self.current_file))

                # Create GeneticData object
                genetic_data = GeneticData(
                    title=parsed_data.get('title', ''),
                    file_path=parsed_data['file_path'],
                    file_name=parsed_data['file_name'],
                    file_type='GenePop',
                    loci=parsed_data['loci'],
                    n_loci=parsed_data['n_loci'],
                    n_samples=parsed_data['n_samples'],
                    n_populations=parsed_data['n_populations'],
                    data=parsed_data['data'],
                    populations=parsed_data['populations']
                )

                self.current_data = genetic_data
                self.data_loaded.emit(genetic_data)

                QMessageBox.information(
                    self,
                    "Success",
                    f"Data loaded successfully!\n\n"
                    f"Samples: {genetic_data.n_samples}\n"
                    f"Loci: {genetic_data.n_loci}\n"
                    f"Populations: {genetic_data.n_populations}\n\n"
                    f"You can now proceed to the Analysis page."
                )

                logger.info(f"Data loaded: {genetic_data}")

            else:
                QMessageBox.warning(
                    self,
                    "Unsupported Format",
                    f"File format {file_ext} is not yet supported.\n"
                    f"Currently supported: GenePop (.gen, .txt)"
                )
                return

        except Exception as e:
            logger.error(f"Error loading data: {e}", exc_info=True)
            QMessageBox.critical(
                self,
                "Error",
                f"Failed to load data:\n{str(e)}\n\n"
                f"Please check that the file is in correct GenePop format."
            )

    def clear_data(self):
        """Clear loaded data"""
        self.current_file = None
        self.current_data = None

        self.drop_area.setText("Drag and drop file here\nor click 'Browse' to select")
        self.info_group.hide()
        self.preview_group.hide()

        self.clear_btn.setEnabled(False)
        self.load_btn.setEnabled(False)

        logger.info("Data cleared")
