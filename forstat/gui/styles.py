"""
Modern styling for Forstat application using PyQt6
"""

class ModernStyle:
    """Modern color scheme and styles"""

    # Color Palette
    PRIMARY = "#2196F3"  # Blue
    PRIMARY_DARK = "#1976D2"
    PRIMARY_LIGHT = "#BBDEFB"

    SECONDARY = "#FF5722"  # Deep Orange
    SECONDARY_DARK = "#E64A19"
    SECONDARY_LIGHT = "#FFCCBC"

    SUCCESS = "#4CAF50"  # Green
    WARNING = "#FF9800"  # Orange
    ERROR = "#F44336"  # Red
    INFO = "#00BCD4"  # Cyan

    # Neutral Colors
    BACKGROUND = "#FAFAFA"
    SURFACE = "#FFFFFF"
    TEXT_PRIMARY = "#212121"
    TEXT_SECONDARY = "#757575"
    BORDER = "#E0E0E0"
    HOVER = "#F5F5F5"

    # Dark Theme
    DARK_BACKGROUND = "#121212"
    DARK_SURFACE = "#1E1E1E"
    DARK_TEXT_PRIMARY = "#FFFFFF"
    DARK_TEXT_SECONDARY = "#B0B0B0"

    @staticmethod
    def get_stylesheet():
        """Get the main application stylesheet"""
        return f"""
            /* Main Window */
            QMainWindow {{
                background-color: {ModernStyle.BACKGROUND};
            }}

            /* Central Widget */
            QWidget {{
                background-color: {ModernStyle.BACKGROUND};
                color: {ModernStyle.TEXT_PRIMARY};
                font-family: 'Segoe UI', Arial, sans-serif;
                font-size: 10pt;
            }}

            /* Push Buttons */
            QPushButton {{
                background-color: {ModernStyle.PRIMARY};
                color: white;
                border: none;
                border-radius: 4px;
                padding: 10px 20px;
                font-weight: bold;
                min-width: 80px;
            }}

            QPushButton:hover {{
                background-color: {ModernStyle.PRIMARY_DARK};
            }}

            QPushButton:pressed {{
                background-color: {ModernStyle.PRIMARY_DARK};
                padding-top: 12px;
            }}

            QPushButton:disabled {{
                background-color: {ModernStyle.BORDER};
                color: {ModernStyle.TEXT_SECONDARY};
            }}

            /* Secondary Button */
            QPushButton[class="secondary"] {{
                background-color: {ModernStyle.SURFACE};
                color: {ModernStyle.PRIMARY};
                border: 2px solid {ModernStyle.PRIMARY};
            }}

            QPushButton[class="secondary"]:hover {{
                background-color: {ModernStyle.PRIMARY_LIGHT};
            }}

            /* Success Button */
            QPushButton[class="success"] {{
                background-color: {ModernStyle.SUCCESS};
            }}

            /* Warning Button */
            QPushButton[class="warning"] {{
                background-color: {ModernStyle.WARNING};
            }}

            /* Danger Button */
            QPushButton[class="danger"] {{
                background-color: {ModernStyle.ERROR};
            }}

            /* Line Edit */
            QLineEdit {{
                background-color: {ModernStyle.SURFACE};
                border: 2px solid {ModernStyle.BORDER};
                border-radius: 4px;
                padding: 8px;
                color: {ModernStyle.TEXT_PRIMARY};
            }}

            QLineEdit:focus {{
                border: 2px solid {ModernStyle.PRIMARY};
            }}

            /* Text Edit */
            QTextEdit, QPlainTextEdit {{
                background-color: {ModernStyle.SURFACE};
                border: 2px solid {ModernStyle.BORDER};
                border-radius: 4px;
                padding: 8px;
                color: {ModernStyle.TEXT_PRIMARY};
            }}

            /* Combo Box */
            QComboBox {{
                background-color: {ModernStyle.SURFACE};
                border: 2px solid {ModernStyle.BORDER};
                border-radius: 4px;
                padding: 8px;
                min-width: 150px;
            }}

            QComboBox:hover {{
                border: 2px solid {ModernStyle.PRIMARY};
            }}

            QComboBox::drop-down {{
                border: none;
                width: 30px;
            }}

            QComboBox::down-arrow {{
                image: none;
                border-left: 5px solid transparent;
                border-right: 5px solid transparent;
                border-top: 5px solid {ModernStyle.TEXT_SECONDARY};
                margin-right: 10px;
            }}

            /* Table Widget */
            QTableWidget {{
                background-color: {ModernStyle.SURFACE};
                border: 1px solid {ModernStyle.BORDER};
                border-radius: 4px;
                gridline-color: {ModernStyle.BORDER};
            }}

            QTableWidget::item {{
                padding: 5px;
            }}

            QTableWidget::item:selected {{
                background-color: {ModernStyle.PRIMARY_LIGHT};
                color: {ModernStyle.TEXT_PRIMARY};
            }}

            QHeaderView::section {{
                background-color: {ModernStyle.SURFACE};
                color: {ModernStyle.TEXT_PRIMARY};
                font-weight: bold;
                border: none;
                border-bottom: 2px solid {ModernStyle.PRIMARY};
                padding: 8px;
            }}

            /* Tab Widget */
            QTabWidget::pane {{
                border: 1px solid {ModernStyle.BORDER};
                border-radius: 4px;
                background-color: {ModernStyle.SURFACE};
            }}

            QTabBar::tab {{
                background-color: {ModernStyle.BACKGROUND};
                color: {ModernStyle.TEXT_SECONDARY};
                border: 1px solid {ModernStyle.BORDER};
                padding: 10px 20px;
                margin-right: 2px;
            }}

            QTabBar::tab:selected {{
                background-color: {ModernStyle.SURFACE};
                color: {ModernStyle.PRIMARY};
                border-bottom: 3px solid {ModernStyle.PRIMARY};
                font-weight: bold;
            }}

            QTabBar::tab:hover:!selected {{
                background-color: {ModernStyle.HOVER};
            }}

            /* Progress Bar */
            QProgressBar {{
                border: 2px solid {ModernStyle.BORDER};
                border-radius: 4px;
                background-color: {ModernStyle.SURFACE};
                text-align: center;
                height: 25px;
            }}

            QProgressBar::chunk {{
                background-color: {ModernStyle.PRIMARY};
                border-radius: 2px;
            }}

            /* Group Box */
            QGroupBox {{
                font-weight: bold;
                border: 2px solid {ModernStyle.BORDER};
                border-radius: 4px;
                margin-top: 10px;
                padding-top: 10px;
            }}

            QGroupBox::title {{
                subcontrol-origin: margin;
                subcontrol-position: top left;
                padding: 0 5px;
                color: {ModernStyle.PRIMARY};
            }}

            /* Scroll Bar */
            QScrollBar:vertical {{
                border: none;
                background: {ModernStyle.BACKGROUND};
                width: 12px;
                margin: 0px;
            }}

            QScrollBar::handle:vertical {{
                background: {ModernStyle.TEXT_SECONDARY};
                min-height: 20px;
                border-radius: 6px;
            }}

            QScrollBar::handle:vertical:hover {{
                background: {ModernStyle.PRIMARY};
            }}

            QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
                border: none;
                background: none;
            }}

            /* Labels */
            QLabel {{
                color: {ModernStyle.TEXT_PRIMARY};
            }}

            QLabel[class="title"] {{
                font-size: 24pt;
                font-weight: bold;
                color: {ModernStyle.PRIMARY};
            }}

            QLabel[class="subtitle"] {{
                font-size: 14pt;
                font-weight: bold;
                color: {ModernStyle.TEXT_PRIMARY};
            }}

            QLabel[class="hint"] {{
                color: {ModernStyle.TEXT_SECONDARY};
                font-size: 9pt;
            }}

            /* Menu Bar */
            QMenuBar {{
                background-color: {ModernStyle.SURFACE};
                border-bottom: 1px solid {ModernStyle.BORDER};
            }}

            QMenuBar::item {{
                padding: 8px 12px;
                background: transparent;
            }}

            QMenuBar::item:selected {{
                background: {ModernStyle.HOVER};
            }}

            QMenu {{
                background-color: {ModernStyle.SURFACE};
                border: 1px solid {ModernStyle.BORDER};
            }}

            QMenu::item {{
                padding: 8px 25px;
            }}

            QMenu::item:selected {{
                background-color: {ModernStyle.PRIMARY_LIGHT};
            }}

            /* Status Bar */
            QStatusBar {{
                background-color: {ModernStyle.SURFACE};
                border-top: 1px solid {ModernStyle.BORDER};
            }}

            /* Tooltips */
            QToolTip {{
                background-color: {ModernStyle.TEXT_PRIMARY};
                color: white;
                border: none;
                padding: 5px;
                border-radius: 3px;
            }}
        """
