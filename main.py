"""
Main entry point for Forstat application
"""
import sys
from PyQt6.QtWidgets import QApplication
from PyQt6.QtGui import QFont

from forstat.gui.main_window import MainWindow
from forstat.gui.styles import ModernStyle
from forstat.utils.config import Config
from forstat.utils.logger import setup_logger

def main():
    """Main application entry point"""
    # Initialize logger
    logger = setup_logger(
        name='forstat',
        log_dir=Config.USER_DIR / 'logs'
    )
    logger.info("="*60)
    logger.info(f"Starting {Config.APP_NAME} v{Config.APP_VERSION}")
    logger.info("="*60)

    # Create application
    app = QApplication(sys.argv)

    # Set application properties
    app.setApplicationName(Config.APP_NAME)
    app.setApplicationVersion(Config.APP_VERSION)
    app.setOrganizationName("Forstat Development Team")

    # Set application font
    font = QFont("Segoe UI", 10)
    app.setFont(font)

    # Apply stylesheet
    app.setStyleSheet(ModernStyle.get_stylesheet())

    # Create and show main window
    window = MainWindow()
    window.show()

    logger.info("Application window displayed")

    # Run application
    exit_code = app.exec()

    logger.info(f"Application exiting with code: {exit_code}")
    return exit_code

if __name__ == '__main__':
    sys.exit(main())
