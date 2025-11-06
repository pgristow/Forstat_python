"""
Main window for Forstat application
"""
from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QStackedWidget, QLabel, QStatusBar, QMessageBox
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QIcon

from forstat.utils.config import Config
from forstat.utils.logger import get_logger

logger = get_logger(__name__)


class MainWindow(QMainWindow):
    """Main application window with page navigation"""

    def __init__(self):
        super().__init__()
        self.config = Config()
        self.config.init_dirs()

        self.current_data = None  # Store loaded data
        self.analysis_results = None  # Store analysis results

        self.init_ui()
        logger.info("Main window initialized")

    def init_ui(self):
        """Initialize the user interface"""
        self.setWindowTitle(Config.APP_TITLE)
        self.setMinimumSize(1200, 800)

        # Center the window on screen
        self.center_on_screen()

        # Create central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        # Main layout
        main_layout = QVBoxLayout(central_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # Header
        header = self.create_header()
        main_layout.addWidget(header)

        # Content area with navigation
        content_layout = QHBoxLayout()
        content_layout.setContentsMargins(0, 0, 0, 0)
        content_layout.setSpacing(0)

        # Navigation sidebar
        nav_widget = self.create_navigation()
        content_layout.addWidget(nav_widget)

        # Stacked widget for pages
        self.page_stack = QStackedWidget()
        content_layout.addWidget(self.page_stack, 1)

        main_layout.addLayout(content_layout)

        # Status bar
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready")

        # Load pages (will be imported)
        self.load_pages()

    def create_header(self):
        """Create application header"""
        header = QWidget()
        header.setObjectName("header")
        header.setStyleSheet(f"""
            #header {{
                background-color: #2196F3;
                padding: 15px;
            }}
        """)

        layout = QHBoxLayout(header)

        # App title
        title = QLabel(f"{Config.APP_NAME}")
        title.setStyleSheet("""
            color: white;
            font-size: 20pt;
            font-weight: bold;
        """)
        layout.addWidget(title)

        # Version
        version = QLabel(f"v{Config.APP_VERSION}")
        version.setStyleSheet("""
            color: rgba(255, 255, 255, 0.7);
            font-size: 10pt;
        """)
        layout.addWidget(version)

        layout.addStretch()

        # Help button
        help_btn = QPushButton("Help")
        help_btn.setStyleSheet("""
            QPushButton {
                background-color: rgba(255, 255, 255, 0.2);
                color: white;
                border: 2px solid white;
                border-radius: 4px;
                padding: 8px 16px;
            }
            QPushButton:hover {
                background-color: rgba(255, 255, 255, 0.3);
            }
        """)
        help_btn.clicked.connect(self.show_help)
        layout.addWidget(help_btn)

        return header

    def create_navigation(self):
        """Create navigation sidebar"""
        nav_widget = QWidget()
        nav_widget.setFixedWidth(200)
        nav_widget.setStyleSheet("""
            QWidget {
                background-color: #f5f5f5;
                border-right: 1px solid #e0e0e0;
            }
        """)

        layout = QVBoxLayout(nav_widget)
        layout.setContentsMargins(10, 20, 10, 10)
        layout.setSpacing(10)

        # Navigation buttons
        self.nav_buttons = []

        btn_upload = self.create_nav_button("1. Upload Data", 0)
        self.nav_buttons.append(btn_upload)
        layout.addWidget(btn_upload)

        btn_analysis = self.create_nav_button("2. Analysis", 1)
        btn_analysis.setEnabled(False)  # Disabled until data is loaded
        self.nav_buttons.append(btn_analysis)
        layout.addWidget(btn_analysis)

        btn_output = self.create_nav_button("3. Results", 2)
        btn_output.setEnabled(False)  # Disabled until analysis is run
        self.nav_buttons.append(btn_output)
        layout.addWidget(btn_output)

        layout.addStretch()

        return nav_widget

    def create_nav_button(self, text, page_index):
        """Create a navigation button"""
        btn = QPushButton(text)
        btn.setCheckable(True)
        btn.setStyleSheet("""
            QPushButton {
                text-align: left;
                padding: 15px;
                border: none;
                border-radius: 4px;
                font-size: 11pt;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #e0e0e0;
            }
            QPushButton:checked {
                background-color: #2196F3;
                color: white;
            }
            QPushButton:disabled {
                color: #bdbdbd;
                background-color: transparent;
            }
        """)
        btn.clicked.connect(lambda: self.switch_page(page_index))

        if page_index == 0:
            btn.setChecked(True)

        return btn

    def load_pages(self):
        """Load application pages"""
        # Import pages here to avoid circular imports
        try:
            from forstat.gui.upload_page import UploadPage
            from forstat.gui.analysis_page import AnalysisPage
            from forstat.gui.output_page import OutputPage

            # Create pages
            self.upload_page = UploadPage(self)
            self.analysis_page = AnalysisPage(self)
            self.output_page = OutputPage(self)

            # Add to stack
            self.page_stack.addWidget(self.upload_page)
            self.page_stack.addWidget(self.analysis_page)
            self.page_stack.addWidget(self.output_page)

            # Connect signals
            self.upload_page.data_loaded.connect(self.on_data_loaded)
            self.analysis_page.analysis_complete.connect(self.on_analysis_complete)

        except ImportError as e:
            logger.error(f"Failed to import pages: {e}")
            # Create placeholder pages for now
            for i in range(3):
                placeholder = QWidget()
                layout = QVBoxLayout(placeholder)
                label = QLabel(f"Page {i+1} - Under Construction")
                label.setAlignment(Qt.AlignmentFlag.AlignCenter)
                label.setStyleSheet("font-size: 16pt; color: #757575;")
                layout.addWidget(label)
                self.page_stack.addWidget(placeholder)

    def switch_page(self, index):
        """Switch to a different page"""
        # Update button states
        for i, btn in enumerate(self.nav_buttons):
            btn.setChecked(i == index)

        # Switch page
        self.page_stack.setCurrentIndex(index)
        logger.info(f"Switched to page {index}")

    def on_data_loaded(self, data):
        """Handle data loaded from upload page"""
        self.current_data = data
        self.nav_buttons[1].setEnabled(True)  # Enable analysis page
        self.status_bar.showMessage("Data loaded successfully")
        logger.info("Data loaded, analysis page enabled")

    def on_analysis_complete(self, results):
        """Handle analysis completion"""
        self.analysis_results = results
        self.nav_buttons[2].setEnabled(True)  # Enable output page

        # Display results in output page
        self.output_page.display_results(results)

        self.status_bar.showMessage("Analysis complete")
        self.switch_page(2)  # Switch to output page
        logger.info("Analysis complete, output page enabled")

    def center_on_screen(self):
        """Center the window on screen"""
        screen = self.screen().geometry()
        x = (screen.width() - self.width()) // 2
        y = (screen.height() - self.height()) // 2
        self.move(x, y)

    def show_help(self):
        """Show help dialog"""
        QMessageBox.information(
            self,
            "Help",
            f"""<h2>{Config.APP_NAME}</h2>
            <p><b>Version:</b> {Config.APP_VERSION}</p>
            <p>Forensic and Population Genetics Analysis Tool</p>
            <h3>How to use:</h3>
            <ol>
                <li><b>Upload Data:</b> Import your genetic data files</li>
                <li><b>Analysis:</b> Select and configure analyses</li>
                <li><b>Results:</b> View and export results</li>
            </ol>
            <p>Supported formats: GenePop, Excel, CSV, FASTA</p>
            """
        )

    def closeEvent(self, event):
        """Handle window close event"""
        reply = QMessageBox.question(
            self,
            'Exit',
            'Are you sure you want to exit?',
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No
        )

        if reply == QMessageBox.StandardButton.Yes:
            logger.info("Application closed by user")
            event.accept()
        else:
            event.ignore()
