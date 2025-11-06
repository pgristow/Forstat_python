@echo off
REM Build script for creating Windows executable

echo ========================================
echo Building Forstat Windows Executable
echo ========================================
echo.

REM Check if PyInstaller is installed
python -c "import PyInstaller" 2>nul
if errorlevel 1 (
    echo Installing PyInstaller...
    pip install pyinstaller
    echo.
)

REM Clean previous builds
echo Cleaning previous builds...
if exist build rmdir /s /q build
if exist dist rmdir /s /q dist
echo.

REM Build executable
echo Building executable...
pyinstaller build_exe.spec

if errorlevel 1 (
    echo.
    echo ========================================
    echo Build FAILED!
    echo ========================================
    pause
    exit /b 1
)

echo.
echo ========================================
echo Build SUCCESSFUL!
echo ========================================
echo.
echo Executable location: dist\Forstat.exe
echo.
pause
