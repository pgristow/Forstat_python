"""
Helper utility functions
"""
import os
from pathlib import Path
from typing import Optional

def format_file_size(size_bytes: int) -> str:
    """
    Format file size in bytes to human-readable format

    Args:
        size_bytes: Size in bytes

    Returns:
        str: Formatted size string
    """
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.2f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.2f} TB"

def get_file_extension(filepath: str) -> str:
    """
    Get file extension in lowercase

    Args:
        filepath: Path to file

    Returns:
        str: File extension (including dot)
    """
    return Path(filepath).suffix.lower()

def validate_file_exists(filepath: str) -> bool:
    """
    Check if file exists and is readable

    Args:
        filepath: Path to file

    Returns:
        bool: True if file exists and is readable
    """
    path = Path(filepath)
    return path.exists() and path.is_file() and os.access(path, os.R_OK)

def sanitize_filename(filename: str) -> str:
    """
    Sanitize filename by removing invalid characters

    Args:
        filename: Original filename

    Returns:
        str: Sanitized filename
    """
    invalid_chars = '<>:"/\\|?*'
    for char in invalid_chars:
        filename = filename.replace(char, '_')
    return filename

def truncate_string(text: str, max_length: int = 50, suffix: str = "...") -> str:
    """
    Truncate string to maximum length

    Args:
        text: String to truncate
        max_length: Maximum length
        suffix: Suffix to add when truncated

    Returns:
        str: Truncated string
    """
    if len(text) <= max_length:
        return text
    return text[:max_length - len(suffix)] + suffix
