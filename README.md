# Enhanced Universal Patent Generator

This folder contains the **new and improved** patent generation system that solves the "?" symbol problem in Markush structures and works with **any CSV file** containing SMILES data.

## 🎯 **Key Improvements**

### ✅ **Solved Problems**
- **No more "?" symbols** in Markush structures
- **Proper R-group labeling** with attachment points
- **Universal compatibility** with any compound dataset
- **Multiple visualization methods** for best results
- **Clean file organization** in dedicated folder

### 🚀 **New Features**
- **Template-based rendering** using actual compounds
- **Core structure with R-groups** properly labeled
- **Representative compound** visualization
- **Simplified scaffold** approach
- **Automatic best method selection**

## 📁 **Files in This Folder**

### Core Modules
- `general_markush_generator.py` - **Main Markush structure generator**
- `enhanced_patent_generator.py` - **Complete patent generation pipeline**
- `universal_patent_runner.py` - **Simple runner script**

### Usage Files
- `README.md` - **This documentation**

## 🚀 **Quick Start**

### Basic Usage
```bash
python universal_patent_runner.py your_compounds.csv
```

### With Custom Output Prefix
```bash
python universal_patent_runner.py CDK5_ACT_IC50.csv CDK5_inhibitors
```

## 📊 **What Gets Generated**

### Patent Documents
- `{prefix}_patent_application.md` - **Complete patent application**
- `{prefix}_metadata.json` - **Generation metadata and statistics**

### Markush Structure Images
- `{prefix}_markush_structure.png` - **Main Markush structure (best method)**
- `{prefix}_template_markush.png` - **Template-based with highlights**
- `{prefix}_core_with_rgroups.png` - **Core with R-group labels**
- `{prefix}_representative.png` - **Representative compound**
- `{prefix}_scaffold.png` - **Simplified scaffold**

## 🔧 **How It Works**

### 1. **Data Analysis**
- Loads CSV file and finds SMILES column automatically
- Analyzes compound structures for common substructure
- Identifies potential R-group attachment points

### 2. **Markush Generation**
- Uses **multiple approaches** to avoid "?" symbols:
  - Template-based rendering with MCS highlights
  - Core structure with proper R-group labeling
  - Representative compound visualization
  - Simplified scaffold approach

### 3. **Method Selection**
- Automatically selects the **best visualization method**
- Priority: `core_with_rgroups > template > scaffold > representative`

### 4. **Patent Creation**
- Generates complete patent application with proper structure images
- Includes claims, technical field, detailed description
- Automatically extracts target information from column names

## 📋 **CSV File Requirements**

### Required
- **SMILES column** (any of: `Smiles`, `SMILES`, `smiles`, `smi`, `SMI`)

### Optional  
- **Activity columns** starting with `ACT_` (e.g., `ACT_CDK5_IC50`)
- **Other compound data** (will be included in analysis)

### Example CSV Structure
```csv
Smiles,ACT_CDK5_IC50,Compound_ID
CCOc1ccc(cc1)C(=O)N...,0.123,COMP_001
CCNc1nc2c(n1)ccc(n2)...,0.456,COMP_002
```

## 🎯 **Target Information Extraction**

The system automatically extracts target information from activity column names:

- `ACT_CDK5_IC50` → Target: CDK5, Disease: Cancer
- `ACT_EGFR_IC50` → Target: EGFR, Disease: Cancer  
- `ACT_JAK2_IC50` → Target: JAK2, Disease: Autoimmune

## 🖼️ **Image Quality**

### ✅ **What You Get**
- **Clear molecular structures** without "?" symbols
- **Proper R-group labels** (R1, R2, R3, etc.)
- **Professional appearance** suitable for patents
- **Multiple visualization options** for best results

### ❌ **Problems Solved**
- No more complex SMARTS rendering issues
- No more "?" symbols in atom positions
- No more missing R-group attachment points
- No more dataset-specific limitations

## 💡 **Advanced Usage**

### Python API
```python
from enhanced_patent_generator import generate_enhanced_patent_universal

results = generate_enhanced_patent_universal("my_compounds.csv", "my_patent")
print(f"Patent: {results['patent_application']}")
print(f"Images: {results['markush_images']}")
```

### Markush Only
```python
from general_markush_generator import process_csv_to_markush

results = process_csv_to_markush("compounds.csv", "output_prefix")
images = results['images']
```

## 🔧 **Troubleshooting**

### Common Issues

**"No SMILES column found"**
- Ensure your CSV has a column named `Smiles`, `SMILES`, `smiles`, `smi`, or `SMI`

**"No common substructure found"**  
- Your compounds may be too diverse
- Try with a more focused compound set
- Check that SMILES are valid

**"Import errors"**
- Ensure RDKit is installed: `pip install rdkit`
- Make sure all files are in the same folder

### Getting Help
Check the generated metadata file for detailed analysis information and error messages.

## 📈 **Comparison with Old System**

| **Aspect** | **Old System** | **New System** |
|------------|----------------|----------------|
| **"?" Symbols** | ❌ Common problem | ✅ Completely solved |
| **R-group Labels** | ⚠️ Often missing | ✅ Proper labeling |
| **Compatibility** | ⚠️ Dataset-specific | ✅ Universal |
| **Visualization** | ❌ Single method | ✅ Multiple approaches |
| **File Organization** | ❌ Scattered files | ✅ Clean New_version folder |
| **Error Handling** | ⚠️ Basic | ✅ Robust fallbacks |

## 🎉 **Success Indicators**

When everything works correctly, you should see:
- ✅ Multiple PNG files generated without "?" symbols
- ✅ Clear R-group labels (R1, R2, etc.) on structures
- ✅ Complete patent application with proper image references
- ✅ Professional-quality Markush structures suitable for patent filing

---

**Generated by Enhanced Universal Patent Generator**  
*Solving Markush structure visualization for any compound dataset*