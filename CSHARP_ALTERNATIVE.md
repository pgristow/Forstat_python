# C# Windows Application Alternative

If you want to create a native Windows application in C# instead of Python, here's what you need to know.

## Why Consider C#?

### Advantages
- ✅ Truly native Windows performance
- ✅ Smaller executable size
- ✅ Better Windows integration
- ✅ No antivirus false positives
- ✅ Professional look and feel
- ✅ Access to full .NET ecosystem
- ✅ Better memory management
- ✅ Easier deployment (no bundling issues)

### Disadvantages
- ❌ Complete rewrite required (significant time investment)
- ❌ Different ecosystem from Python scientific libraries
- ❌ Learning curve if not familiar with C#
- ❌ Less robust scientific computing libraries than Python

## Technology Stack for C#

### Option 1: WPF (Windows Presentation Foundation)
**Best for: Modern, rich UI applications**

```csharp
// Modern XAML-based UI
<Window x:Class="Forstat.MainWindow"
        xmlns="http://schemas.microsoft.com/winfx/2006/xaml/presentation"
        Title="Forstat" Height="800" Width="1200">
    <Grid>
        <!-- Your UI here -->
    </Grid>
</Window>
```

**Pros:**
- Beautiful, modern UI
- MVVM pattern support
- Data binding
- Vector graphics

### Option 2: WinForms
**Best for: Simpler, traditional applications**

```csharp
// Drag-and-drop designer
public partial class MainForm : Form
{
    public MainForm()
    {
        InitializeComponent();
    }
}
```

**Pros:**
- Easier to learn
- Rapid development
- Stable and mature

### Option 3: MAUI (Multi-platform App UI)
**Best for: Cross-platform (Windows, Mac, Linux)**

```csharp
// Modern, cross-platform
public class MainPage : ContentPage
{
    public MainPage()
    {
        Content = new StackLayout { /* ... */ };
    }
}
```

### Option 4: Avalonia UI
**Best for: Cross-platform with modern XAML**

Similar to WPF but runs on Windows, Mac, Linux.

## Scientific Computing in C#

### Math.NET Numerics
```csharp
using MathNet.Numerics;
using MathNet.Numerics.Statistics;

// Statistical calculations
var data = new[] { 1.0, 2.0, 3.0, 4.0, 5.0 };
var mean = Statistics.Mean(data);
var stdDev = Statistics.StandardDeviation(data);
```

### Accord.NET Framework
```csharp
using Accord.Statistics;
using Accord.Statistics.Analysis;

// Population genetics calculations
var hwe = new HardyWeinbergTest(observed, expected);
bool isEquilibrium = hwe.Significant;
```

### IronPython Integration
**Best of both worlds!**

```csharp
// Embed Python in C# app
using IronPython.Hosting;

var engine = Python.CreateEngine();
engine.Execute("import numpy as np");
dynamic np = engine.Runtime.ImportModule("numpy");
var array = np.array(new[] { 1, 2, 3 });
```

## Excel Export in C#

### EPPlus (Excel library)
```csharp
using OfficeOpenXml;

using (var package = new ExcelPackage())
{
    var worksheet = package.Workbook.Worksheets.Add("Results");
    worksheet.Cells["A1"].Value = "Sample ID";
    worksheet.Cells["B1"].Value = "Allele Frequency";

    package.SaveAs(new FileInfo("results.xlsx"));
}
```

## PDF Generation in C#

### iTextSharp or PdfSharp
```csharp
using PdfSharp.Pdf;
using PdfSharp.Drawing;

var document = new PdfDocument();
var page = document.AddPage();
var gfx = XGraphics.FromPdfPage(page);
gfx.DrawString("Forstat Report", font, XBrushes.Black,
    new XRect(0, 0, page.Width, page.Height),
    XStringFormats.TopCenter);
document.Save("report.pdf");
```

## Database Integration

### Entity Framework
```csharp
public class ForstatContext : DbContext
{
    public DbSet<Sample> Samples { get; set; }
    public DbSet<Analysis> Analyses { get; set; }
}

// Query data
var samples = context.Samples
    .Where(s => s.Population == "EUR")
    .ToList();
```

## Project Structure

```
ForstatCSharp/
├── Forstat.sln                    # Visual Studio solution
├── Forstat/                       # Main application
│   ├── App.xaml                   # Application config
│   ├── MainWindow.xaml            # Main window UI
│   ├── Views/                     # UI views
│   │   ├── UploadView.xaml
│   │   ├── AnalysisView.xaml
│   │   └── ResultsView.xaml
│   ├── ViewModels/                # MVVM view models
│   ├── Models/                    # Data models
│   ├── Services/                  # Business logic
│   │   ├── DataParserService.cs
│   │   ├── AnalysisService.cs
│   │   └── ExportService.cs
│   └── Utils/                     # Utilities
├── Forstat.Core/                  # Core library (portable)
│   ├── Genetics/
│   │   ├── PopulationGenetics.cs
│   │   ├── ForensicStats.cs
│   │   └── STRAnalysis.cs
│   └── Parsers/
│       ├── GenepopParser.cs
│       └── FastaParser.cs
└── Forstat.Tests/                 # Unit tests
```

## Building C# Executable

### Using Visual Studio
1. Right-click project → Publish
2. Select "Folder" target
3. Choose "Self-contained" deployment
4. Click "Publish"

Result: Single .exe file with all dependencies

### Using .NET CLI
```bash
dotnet publish -c Release -r win-x64 --self-contained true -p:PublishSingleFile=true
```

Output: Single .exe (~100MB), no dependencies needed

## Installation Distribution

### MSI Installer (Professional)
- Use WiX Toolset
- Windows Installer format
- Proper uninstall support

### ClickOnce (Easy updates)
- Built into Visual Studio
- Automatic updates
- Web deployment

### Microsoft Store
- Distribute via Windows Store
- Automatic updates
- Built-in payment

## Hybrid Approach: Python + C# Wrapper

**Best of both worlds!**

### Option A: C# GUI + Python Backend
```csharp
// Call Python from C#
using Python.Runtime;

class AnalysisEngine
{
    public void RunAnalysis()
    {
        using (Py.GIL())
        {
            dynamic scipy = Py.Import("scipy.stats");
            var result = scipy.ttest_ind(data1, data2);
        }
    }
}
```

### Option B: IronPython
```csharp
// Embed Python in C#
var engine = Python.CreateEngine();
var scope = engine.CreateScope();
engine.ExecuteFile("analysis.py", scope);
dynamic result = scope.GetVariable("result");
```

## Recommendation

### Choose Python + PyInstaller if:
- ✅ You want to keep existing code
- ✅ You need Python's scientific libraries
- ✅ Quick time to market
- ✅ You're comfortable with Python
- ✅ File size isn't critical

### Choose C# if:
- ✅ You want native Windows performance
- ✅ Professional enterprise deployment
- ✅ Better Windows integration
- ✅ Long-term maintenance
- ✅ Smaller file size
- ✅ You have time for rewrite

### Choose Hybrid (C# + Python) if:
- ✅ You want native UI performance
- ✅ You need Python libraries
- ✅ Best of both worlds
- ✅ Gradual migration path

## Time Estimates

| Approach | Development Time | Result |
|----------|-----------------|---------|
| **PyInstaller** | 1-2 hours | 100-200MB .exe |
| **C# Full Rewrite** | 3-6 months | 50MB .exe, native |
| **C# + Python Hybrid** | 1-2 months | 80-150MB .exe |

## My Recommendation for You

**Stick with Python + PyInstaller** because:
1. Your app is already built in Python
2. Scientific computing is better in Python
3. 2-3 hours to create distributable .exe
4. Users won't notice performance difference for this use case
5. Focus on features, not rewriting

Later, if needed:
- Add IronPython for better Windows integration
- Or gradually port performance-critical parts to C#

## Resources

### If you decide on C#:
- Visual Studio Community (Free)
- .NET 8.0 SDK
- Math.NET Numerics
- Accord.NET
- EPPlus for Excel
- WPF tutorial: https://docs.microsoft.com/wpf

### For hybrid approach:
- Python.NET: https://github.com/pythonnet/pythonnet
- IronPython: https://ironpython.net/
