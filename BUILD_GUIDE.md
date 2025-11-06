# Building Windows Executable

## Method 1: PyInstaller (Recommended)

### Quick Build

Single command to create executable:

```bash
# Install PyInstaller
pip install pyinstaller

# Create executable (one file, no console)
pyinstaller --onefile --windowed --name Forstat main.py
```

The executable will be in `dist/Forstat.exe`

### Using Build Script

For production builds with customization:

```bash
# Run the build script
build.bat
```

Or manually:

```bash
pyinstaller build_exe.spec
```

### Build Options Explained

- `--onefile`: Bundle everything into a single .exe file
- `--windowed`: No console window (for GUI apps)
- `--name Forstat`: Name of the executable
- `--icon=app.ico`: Add custom icon
- `--add-data`: Include additional files (data, resources)

### Advanced Configuration

Edit `build_exe.spec` for advanced options:
- Include/exclude specific modules
- Add hidden imports
- Customize file paths
- Set application metadata

## Method 2: cx_Freeze

Alternative to PyInstaller:

```bash
pip install cx_Freeze
python setup.py build_exe
```

## Method 3: py2exe

Windows-specific option:

```bash
pip install py2exe
python setup_py2exe.py
```

## Method 4: Nuitka

Compiles Python to C, then to executable (fastest performance):

```bash
pip install nuitka
python -m nuitka --standalone --windows-disable-console --enable-plugin=pyqt6 main.py
```

## Comparison

| Tool | Pros | Cons |
|------|------|------|
| **PyInstaller** | Easy, popular, good docs | Large file size (~100MB+) |
| **cx_Freeze** | Cross-platform, reliable | Complex configuration |
| **py2exe** | Windows-optimized | Windows only, outdated |
| **Nuitka** | Best performance, smaller | Longer compile time |

## File Size Optimization

Reduce executable size:

```bash
# Use UPX compression
pyinstaller --onefile --windowed --upx-dir=/path/to/upx main.py

# Exclude unnecessary modules
pyinstaller --exclude-module matplotlib main.py
```

## Installer Creation

After creating .exe, create an installer:

### Using Inno Setup (Free)

1. Download Inno Setup
2. Create installer script
3. Includes: start menu shortcuts, desktop icons, uninstaller

### Using NSIS (Free)

1. Download NSIS
2. Create .nsi script
3. Professional installer with wizard

### Using WiX Toolset (Free, Microsoft)

1. Download WiX
2. Create MSI installer
3. Enterprise-grade deployment

## Recommended Workflow

1. **Development**: Use Python directly
2. **Testing**: Create PyInstaller executable
3. **Distribution**:
   - Package with Inno Setup/NSIS
   - Or distribute standalone .exe
   - Optional: Code signing certificate

## Troubleshooting

### Missing DLL errors
- Use `--collect-all package-name`
- Or add to `hiddenimports` in spec file

### Import errors
- Add to `hiddenimports` list

### Large file size
- Use `--exclude-module` for unused packages
- Enable UPX compression

### Antivirus warnings
- Sign your executable with code signing certificate
- Or submit to antivirus companies for whitelisting
