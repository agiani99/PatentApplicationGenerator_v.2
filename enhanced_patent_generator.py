"""
Enhanced Patent Generator with General Markush Structure Support

This module integrates the new general Markush structure generator
with the patent generation pipeline, providing a universal solution
for any compound dataset.
"""

import pandas as pd
import json
from datetime import datetime
from pathlib import Path
import logging
from typing import Dict, Optional, List, Any

# Import the general Markush generator
try:
    from .general_markush_generator import (
        GeneralMarkushGenerator, 
        process_csv_to_markush,
        get_enhanced_markush_images
    )
except ImportError:
    from general_markush_generator import (
        GeneralMarkushGenerator,
        process_csv_to_markush, 
        get_enhanced_markush_images
    )

# Import patent sections with prior art search
try:
    from patent_sections_integrated import generate_technical_field, generate_background_prior_art
    PATENT_SECTIONS_AVAILABLE = True
except ImportError:
    print("⚠️  patent_sections_integrated not found - prior art search will be skipped")
    PATENT_SECTIONS_AVAILABLE = False
    
logger = logging.getLogger(__name__)

class EnhancedPatentGenerator:
    """
    Enhanced patent generator using the new general Markush approach.
    
    This generator:
    1. Works with any CSV file containing SMILES
    2. Generates proper Markush structures with R-group labels
    3. Avoids the '?' symbol problem
    4. Provides multiple visualization methods
    """
    
    def __init__(self):
        self.compound_data = None
        self.markush_results = None
        self.target_info = None
        self.activity_analysis = None
        self.images = {}
        
    def generate_patent_from_csv(self, csv_path: str, output_prefix: Optional[str] = None) -> Dict[str, Any]:
        """
        Generate complete patent application from CSV file.
        
        Args:
            csv_path: Path to CSV file with SMILES and activity data
            output_prefix: Prefix for output files
            
        Returns:
            Dictionary with generated file paths and metadata
        """
        
        csv_file = Path(csv_path)
        if not csv_file.exists():
            raise FileNotFoundError(f"CSV file not found: {csv_path}")
        
        if output_prefix is None:
            output_prefix = csv_file.stem
        
        logger.info(f"🚀 Starting enhanced patent generation for {csv_path}")
        
        # Ensure output directory exists
        output_dir = Path("./New_version")
        output_dir.mkdir(exist_ok=True)
        
        # Step 1: Load and analyze compound data
        logger.info("📊 Step 1: Loading and analyzing compound data")
        self._load_compound_data(csv_path)
        
        # Step 2: Generate Markush structures
        logger.info("🧬 Step 2: Generating Markush structures")
        self._generate_markush_structures(csv_path, output_prefix)
        
        # Step 3: Analyze target and activity data
        logger.info("🎯 Step 3: Analyzing target and activity data")
        self._analyze_target_and_activity()
        
        # Step 4: Generate patent document
        logger.info("📄 Step 4: Generating patent document")
        patent_file = self._generate_patent_document(output_prefix)
        
        # Step 5: Generate metadata
        logger.info("📋 Step 5: Generating metadata")
        metadata_file = self._generate_metadata(output_prefix)
        
        # Compile results
        results = {
            'csv_file': str(csv_file),
            'output_prefix': output_prefix,
            'patent_application': patent_file,
            'metadata_file': metadata_file,
            'markush_images': self.images,
            'markush_analysis': self.markush_results,
            'target_info': self.target_info,
            'activity_analysis': self.activity_analysis,
            'generation_timestamp': datetime.now().isoformat(),
            'success': True
        }
        
        logger.info(f"🎉 Enhanced patent generation complete!")
        return results
    
    def _load_compound_data(self, csv_path: str):
        """Load and validate compound data from CSV."""
        
        try:
            self.compound_data = pd.read_csv(csv_path)
            logger.info(f"✅ Loaded {len(self.compound_data)} compounds from CSV")
            
            # Find SMILES column with flexible matching
            smiles_col = None
            
            # First, try exact matches
            exact_matches = ['Smiles', 'SMILES', 'smiles', 'smi', 'SMI']
            for col in exact_matches:
                if col in self.compound_data.columns:
                    smiles_col = col
                    break
            
            # If no exact match, try case-insensitive and stripped matches
            if not smiles_col:
                for col in self.compound_data.columns:
                    col_clean = col.strip().lower()
                    if col_clean in ['smiles', 'smi']:
                        smiles_col = col  # Use original column name
                        break
            
            # If still no match, look for columns containing 'smiles'
            if not smiles_col:
                for col in self.compound_data.columns:
                    col_clean = col.strip().lower()
                    if 'smiles' in col_clean or 'smi' in col_clean:
                        smiles_col = col  # Use original column name
                        break
            
            if not smiles_col:
                raise ValueError(f"No SMILES column found. Available columns: {list(self.compound_data.columns)}")
            
            valid_smiles = self.compound_data[smiles_col].dropna()
            logger.info(f"✅ Found {len(valid_smiles)} valid SMILES strings in column '{smiles_col}'")
            
        except Exception as e:
            logger.error(f"❌ Error loading compound data: {e}")
            raise
    
    def _generate_markush_structures(self, csv_path: str, output_prefix: str):
        """Generate Markush structures using the new general approach."""
        
        try:
            # Use the general Markush generator
            self.markush_results = process_csv_to_markush(csv_path, output_prefix=output_prefix)
            
            if "error" in self.markush_results:
                raise ValueError(f"Markush generation failed: {self.markush_results['error']}")
            
            # Get the image paths
            self.images = self.markush_results.get('images', {})
            
            logger.info(f"✅ Generated {len(self.images)} Markush structure images")
            
            # Log the generated images
            for method, path in self.images.items():
                if Path(path).exists():
                    logger.info(f"  ✅ {method}: {path}")
                else:
                    logger.warning(f"  ❌ {method}: {path} (not found)")
            
        except Exception as e:
            logger.error(f"❌ Error generating Markush structures: {e}")
            raise
    
    def _analyze_target_and_activity(self):
        """Analyze target information and activity data."""
        
        try:
            # Extract target information from column names
            self.target_info = self._extract_target_info()
            
            # Analyze activity data
            self.activity_analysis = self._analyze_activity_data()
            
            logger.info(f"✅ Target: {self.target_info.get('target_name', 'Unknown')}")
            logger.info(f"✅ Activity data: {self.activity_analysis.get('compounds_with_activity', 0)} compounds")
            
        except Exception as e:
            logger.error(f"❌ Error analyzing target/activity data: {e}")
            # Continue with defaults
            self.target_info = {'target_name': 'target protein', 'diseases': ['diseases']}
            self.activity_analysis = {'compounds_with_activity': len(self.compound_data) if self.compound_data is not None else 0}
    
    def _extract_target_info(self) -> Dict[str, Any]:
        """Extract target information from the data."""
        
        target_info = {
            'target_name': 'target protein',
            'target_type': 'Unknown',
            'organism': 'Human',
            'diseases': ['diseases'],
        }
        
        if self.compound_data is None:
            return target_info
        
        # Look for activity columns (ACT_*)
        activity_cols = [col for col in self.compound_data.columns if col.upper().startswith('ACT_')]
        
        if activity_cols:
            # Extract target name from first activity column
            target_name = activity_cols[0].split('_', 1)[1] if '_' in activity_cols[0] else activity_cols[0]
            target_info['target_name'] = target_name
            
            # Infer disease areas and target type based on name
            target_name_upper = target_name.upper()
            
            # Disease inference
            if any(term in target_name_upper for term in ['CDK', 'CYCLIN']):
                target_info['diseases'] = ['cancer', 'oncology']
                target_info['target_type'] = 'Kinase'
            elif any(term in target_name_upper for term in ['EGFR', 'HER', 'ERB']):
                target_info['diseases'] = ['cancer', 'solid tumors']
                target_info['target_type'] = 'Receptor Tyrosine Kinase'
            elif any(term in target_name_upper for term in ['VEGFR', 'KDR', 'FLT']):
                target_info['diseases'] = ['cancer', 'angiogenesis-related diseases']
                target_info['target_type'] = 'Receptor Tyrosine Kinase'
            elif any(term in target_name_upper for term in ['JAK', 'STAT']):
                target_info['diseases'] = ['autoimmune diseases', 'inflammatory diseases']
                target_info['target_type'] = 'Kinase'
            elif any(term in target_name_upper for term in ['KINASE', 'CDK']):
                target_info['target_type'] = 'Kinase'
        
        return target_info
    
    def _analyze_activity_data(self) -> Dict[str, Any]:
        """Analyze activity data in the dataset."""
        
        analysis = {
            'activity_column': None,
            'compounds_with_activity': 0,
            'activity_stats': {
                'min': 0.0, 'max': 0.0, 'median': 0.0,
                'mean': 0.0, 'std': 0.0, 'q1': 0.0, 'q3': 0.0
            },
            'top_compounds_count': 0
        }
        
        if self.compound_data is None:
            return analysis
        
        # Find activity columns
        activity_cols = [col for col in self.compound_data.columns if col.upper().startswith('ACT_')]
        
        if activity_cols:
            activity_col = activity_cols[0]
            analysis['activity_column'] = activity_col
            
            activity_data = self.compound_data[activity_col].dropna()
            analysis['compounds_with_activity'] = len(activity_data)
            
            if len(activity_data) > 0:
                analysis['activity_stats'] = {
                    'min': float(activity_data.min()),
                    'max': float(activity_data.max()),
                    'median': float(activity_data.median()),
                    'mean': float(activity_data.mean()),
                    'std': float(activity_data.std()),
                    'q1': float(activity_data.quantile(0.25)),
                    'q3': float(activity_data.quantile(0.75))
                }
                
                # Count highly active compounds (bottom 10% of values)
                threshold = activity_data.quantile(0.1)
                analysis['top_compounds_count'] = len(activity_data[activity_data <= threshold])
        
        return analysis
    
    def _generate_patent_document(self, output_prefix: str) -> str:
        """Generate the complete patent document."""
        
        output_dir = Path("./New_version")
        patent_file = output_dir / f"{output_prefix}_patent_application.md"
        
        try:
            # Copy main image to same directory as markdown for relative paths
            if self.images and 'template' in self.images:
                template_path = Path(self.images['template'])
                if template_path.exists():
                    target_image = output_dir / template_path.name
                    if not target_image.exists() or template_path != target_image:
                        import shutil
                        shutil.copy2(template_path, target_image)
                        logger.info(f"Copied image to markdown directory: {target_image.name}")
            
            # Generate all sections
            content = self._create_patent_content()
            
            # Write to file
            with open(patent_file, 'w', encoding='utf-8') as f:
                f.write(content)
            
            logger.info(f"✅ Patent document saved: {patent_file}")
            return str(patent_file)
            
        except Exception as e:
            logger.error(f"❌ Error generating patent document: {e}")
            raise
    
    def _create_patent_content(self) -> str:
        """Create the complete patent document content."""
        
        # Ensure target_info and activity_analysis have defaults
        if self.target_info is None:
            self.target_info = {'target_name': 'target protein', 'diseases': ['diseases'], 'target_type': 'protein'}
        
        if self.activity_analysis is None:
            self.activity_analysis = {'compounds_with_activity': 0, 'top_compounds_count': 0}
        
        target_name = self.target_info.get('target_name', 'target protein')
        diseases = self.target_info.get('diseases', ['diseases'])
        
        # Get the main Markush structure image - prioritize template method
        main_image = ''
        
        # First try template method (shows rings correctly)
        if 'template' in self.images and Path(self.images['template']).exists():
            # Use just the filename for relative path
            main_image = Path(self.images['template']).name
            logger.info("Using template method for main image (preserves ring structures)")
        
        # Fallback hierarchy if template not available
        elif 'main' in self.images and Path(self.images['main']).exists():
            main_image = Path(self.images['main']).name
        
        # If no good options, use any available but warn
        else:
            for key in ['core_with_rgroups', 'scaffold', 'representative']:
                if key in self.images and Path(self.images[key]).exists():
                    main_image = Path(self.images[key]).name
                    logger.warning(f"Using {key} method as fallback - may not show ring structures correctly")
                    break
        
        # Generate document sections with proper prior art search
        title = f"{target_name.upper()} INHIBITORS AND METHODS OF USE THEREOF"
        
        # Technical field section
        if PATENT_SECTIONS_AVAILABLE:
            technical_field = generate_technical_field(target_name)
            logger.info("✅ Generated technical field with patent_sections_integrated")
        else:
            technical_field = self._generate_basic_technical_field(target_name, diseases)
            logger.warning("⚠️  Using basic technical field - patent_sections_integrated not available")
        
        # Background and prior art section  
        if PATENT_SECTIONS_AVAILABLE:
            # Extract UniProt ID or use target name for prior art search
            uniprot_id = self._extract_uniprot_id()
            background = generate_background_prior_art(target_name, uniprot_id)
            logger.info(f"✅ Generated background with prior art search for {target_name} (UniProt: {uniprot_id})")
        else:
            background = self._generate_basic_background(target_name, diseases)
            logger.warning("⚠️  Using basic background - prior art search not available")
        
        summary = f"""## 3. Summary of the Invention

The present invention provides **novel {target_name} inhibitors** represented by the general **Formula I** below, along with pharmaceutical compositions and therapeutic methods.

### Formula I - General Markush Structure

![Formula I]({main_image})

**Figure 1: Formula I** - General structure showing the core scaffold with variable R-group positions for systematic optimization.

### 🚀 **Key Innovation Aspects**

| **Feature** | **Benefit** |
|-------------|-------------|
| 🎯 **Systematic Markush Approach** | Comprehensive R-group exploration |
| 🧪 **Core Scaffold Optimization** | Enhanced {target_name} binding |
| 💊 **Drug-like Properties** | Optimized ADMET characteristics |
| 🎭 **Improved Selectivity** | Reduced off-target effects |

### 📊 **Compound Library Overview**
- **Total Compounds**: {self.activity_analysis.get('compounds_with_activity', 0)}
- **Active Compounds**: {self.activity_analysis.get('top_compounds_count', 0)}
- **R-group Variations**: Multiple positions for systematic optimization
"""

        claims = f"""## 6. Claims

**1.** A compound of **Formula I**:

![Formula I]({main_image})

wherein R-groups represent variable substituents selected to optimize {target_name} activity and pharmaceutical properties.

**2.** The compound of claim 1, wherein the compound exhibits **{target_name} inhibitory activity** with an IC50 value of less than **10 μM**.

**3.** The compound of claim 1, wherein the compound exhibits **{target_name} inhibitory activity** with an IC50 value of less than **1 μM**.

**4.** A **pharmaceutical composition** comprising:
- A therapeutically effective amount of a compound of claim 1; and
- A pharmaceutically acceptable carrier.

**5.** A **method of treating** {', '.join(diseases)} comprising administering a therapeutically effective amount of a compound of claim 1.

**6.** A **method of inhibiting {target_name} activity** comprising contacting {target_name} with a compound of claim 1.
"""

        # Complete document
        full_document = f"""# {title}

{technical_field}

{background}

{summary}

## 4. Brief Description of Structures

The following figure illustrates the key structural elements of the invention:

### Markush Structure - Formula I

![Formula I - Markush Structure]({main_image})

**Figure 1:** General Markush structure showing the core scaffold with ring systems and variable R-group positions for systematic optimization.

> **Note:** The structure shows the common ring-containing core that is essential for {target_name} activity, with attachment points for R-group variations.

## 5. Detailed Description

The compounds of Formula I represent a novel class of **{target_name} inhibitors** with systematically optimized substituents providing enhanced potency, selectivity, and pharmaceutical properties.

### 📊 **Activity Summary**
- **Compounds Tested**: {self.activity_analysis.get('compounds_with_activity', 0)}
- **Target**: {target_name}
- **Therapeutic Areas**: {', '.join(diseases)}

{claims}

## 7. Abstract

Novel **{target_name} inhibitors** of Formula I for treating **{', '.join(diseases)}**. The compounds feature a systematic Markush approach with optimized R-group substituents providing superior activity, selectivity, and drug-like properties.

---

## Document Information

- **Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
- **Method**: Enhanced General Markush Generator with Prior Art Search
- **Images Generated**: {len(self.images)}
- **Target**: {target_name}
- **Prior Art Search**: {"✅ Included" if PATENT_SECTIONS_AVAILABLE else "❌ Skipped (patent_sections_integrated not available)"}

> **Note**: This patent application was generated using an enhanced approach that includes comprehensive prior art search and avoids visualization artifacts.
"""
        
        return full_document
    
    def _extract_uniprot_id(self) -> str:
        """Extract UniProt ID from activity columns or use default."""
        if self.compound_data is None:
            return "P24941"  # Default CDK5
        
        # Look for activity columns and extract target identifier
        activity_cols = [col for col in self.compound_data.columns if col.upper().strip().startswith('ACT_')]
        
        if activity_cols:
            # Extract the part after ACT_
            target_part = activity_cols[0].split('_', 1)[1] if '_' in activity_cols[0] else activity_cols[0]
            
            # Common target to UniProt ID mappings
            uniprot_map = {
                'CDK5': 'P24941',
                'KCNQ2': 'O43526', 
                'EGFR': 'P00533',
                'KDR': 'P35968',
                'VEGFR2': 'P35968',
                'JAK2': 'O60674',
                'CDK2': 'P24941'
            }
            
            # Try exact match first
            if target_part.upper() in uniprot_map:
                return uniprot_map[target_part.upper()]
            
            # Try partial match
            for target, uniprot in uniprot_map.items():
                if target in target_part.upper():
                    return uniprot
            
            # Return the target name itself (might work for some searches)
            return target_part
        
        return "P24941"  # Default fallback
    
    def _generate_basic_technical_field(self, target_name: str, diseases: List[str]) -> str:
        """Generate basic technical field when patent_sections_integrated is not available."""
        return f"""## 2. Technical Field

The present invention relates to **novel chemical compounds** that act as **{target_name} modulators**, pharmaceutical compositions containing such compounds, and **methods of using these compounds** for therapeutic purposes.

### 🎯 **Field of Application**
- **Target**: {target_name} ({self.target_info.get('target_type', 'protein') if self.target_info else 'protein'})
- **Therapeutic Areas**: {', '.join(diseases)}
- **Innovation**: Systematic Markush structure approach with optimized R-group substituents

### 📋 **Technical Domain**
The invention falls within the fields of:
- **Medicinal Chemistry**: Novel small molecule inhibitors
- **Pharmaceutical Sciences**: Drug discovery and development
- **Biochemistry**: Protein-drug interactions
- **Therapeutic Applications**: Treatment of {', '.join(diseases)}
"""
    
    def _generate_basic_background(self, target_name: str, diseases: List[str]) -> str:
        """Generate basic background when patent_sections_integrated is not available."""
        return f"""## 3. Background Art

### 🎯 **Target Background**

**{target_name}** is an important therapeutic target involved in {', '.join(diseases)}. Current therapeutic approaches have limitations including:

- **Limited efficacy** in certain patient populations
- **Side effects** due to off-target activity  
- **Resistance development** with existing therapies
- **Suboptimal pharmacokinetic properties**

### 🔬 **Prior Art Limitations**

Existing {target_name} modulators suffer from:

| **Limitation** | **Impact** |
|----------------|------------|
| **Poor selectivity** | Off-target effects |
| **Limited potency** | Higher doses required |
| **ADMET issues** | Poor drug-like properties |
| **Narrow SAR** | Limited optimization potential |

### 💡 **Need for Innovation**

There remains a significant need for **novel {target_name} inhibitors** with:
- ✅ **Enhanced potency** and selectivity
- ✅ **Improved pharmacokinetic properties**
- ✅ **Reduced side effects**
- ✅ **Systematic optimization potential**

> **Note**: Comprehensive prior art search requires patent_sections_integrated module. Please ensure this module is available for complete patent applications.
"""
    
    def _generate_metadata(self, output_prefix: str) -> str:
        """Generate metadata file."""
        
        output_dir = Path("./New_version")
        metadata_file = output_dir / f"{output_prefix}_metadata.json"
        
        # Ensure all attributes have defaults
        if self.target_info is None:
            self.target_info = {'target_name': 'target protein', 'diseases': ['diseases']}
        
        if self.activity_analysis is None:
            self.activity_analysis = {'compounds_with_activity': 0, 'activity_column': None}
        
        metadata = {
            'generation_info': {
                'timestamp': datetime.now().isoformat(),
                'method': 'Enhanced General Markush Generator',
                'output_prefix': output_prefix
            },
            'compound_data': {
                'total_compounds': len(self.compound_data) if self.compound_data is not None else 0,
                'compounds_with_activity': self.activity_analysis.get('compounds_with_activity', 0),
                'activity_column': self.activity_analysis.get('activity_column', None)
            },
            'markush_analysis': self.markush_results.get('analysis', {}) if self.markush_results else {},
            'target_info': self.target_info,
            'activity_analysis': self.activity_analysis,
            'generated_images': self.images,
            'success': True
        }
        
        try:
            with open(metadata_file, 'w', encoding='utf-8') as f:
                json.dump(metadata, f, indent=2, default=str)
            
            logger.info(f"✅ Metadata saved: {metadata_file}")
            return str(metadata_file)
            
        except Exception as e:
            logger.error(f"❌ Error saving metadata: {e}")
            raise

# Main function for easy usage

def generate_enhanced_patent_universal(csv_path: str, output_prefix: Optional[str] = None) -> Dict[str, Any]:
    """
    Universal function to generate enhanced patents from any CSV file.
    
    Args:
        csv_path: Path to CSV file with SMILES and optional activity data
        output_prefix: Output prefix (derived from filename if None)
        
    Returns:
        Dictionary with generation results
    """
    
    generator = EnhancedPatentGenerator()
    results = generator.generate_patent_from_csv(csv_path, output_prefix)
    
    print(f"\n🎉 ENHANCED PATENT GENERATION COMPLETE!")
    print(f"=" * 50)
    print(f"📁 Input: {results['csv_file']}")
    print(f"🎯 Target: {results['target_info'].get('target_name', 'Unknown')}")
    print(f"🧪 Compounds: {results['activity_analysis'].get('compounds_with_activity', 0)}")
    print(f"🖼️  Images: {len(results['markush_images'])}")
    print(f"📄 Patent: {results['patent_application']}")
    print(f"📋 Metadata: {results['metadata_file']}")
    
    return results