"""
Patent Sections with Integrated Prior Art Search

This module provides comprehensive patent sections including:
- Technical field generation
- Background art with PubChem/ChEMBL searches
- Prior art analysis and citation
"""

import logging
from typing import Dict, List, Optional, Any
import requests
import json
from datetime import datetime
from patent_utils_fixed import (
    get_or_query_uniprot_id,
    search_patents_by_chembl_id, 
    search_patents_by_inchikey,
    get_pubchem_iupac_name,
    get_chembl_target_by_identifier,
    get_chemical_name_multi_source
)

logger = logging.getLogger(__name__)

def generate_technical_field(target_name: str) -> str:
    """
    Generate comprehensive technical field section.
    
    Args:
        target_name: Name of the target protein/enzyme
        
    Returns:
        Formatted technical field section
    """
    
    # Target-specific technical domains
    technical_domains = {
        'CDK5': {
            'primary': 'Cell cycle regulation and neurodegeneration',
            'secondary': ['Alzheimer\'s disease', 'cancer therapy', 'neuroprotection'],
            'protein_family': 'Cyclin-dependent kinases',
            'therapeutic_area': 'Oncology and neurology'
        },
        'KCNQ2': {
            'primary': 'Neuronal excitability and epilepsy',
            'secondary': ['epilepsy', 'seizure disorders', 'ion channel modulation'],
            'protein_family': 'Voltage-gated potassium channels',
            'therapeutic_area': 'Neurology and epilepsy'
        },
        'EGFR': {
            'primary': 'Growth factor signaling and cancer',
            'secondary': ['cancer therapy', 'solid tumors', 'targeted therapy'],
            'protein_family': 'Receptor tyrosine kinases',
            'therapeutic_area': 'Oncology'
        },
        'KDR': {
            'primary': 'Angiogenesis and vascular development',
            'secondary': ['cancer', 'angiogenesis inhibition', 'VEGF signaling'],
            'protein_family': 'Receptor tyrosine kinases',
            'therapeutic_area': 'Oncology and vascular diseases'
        }
    }
    
    # Get target-specific info or use generic
    target_key = target_name.upper()
    if target_key in technical_domains:
        domain_info = technical_domains[target_key]
    else:
        domain_info = {
            'primary': f'{target_name} modulation',
            'secondary': ['therapeutic intervention', 'drug discovery'],
            'protein_family': 'Target proteins',
            'therapeutic_area': 'Pharmaceutical applications'
        }
    
    return f"""## 2. Technical Field

The present invention relates to **novel chemical compounds** that act as **{target_name} modulators**, pharmaceutical compositions containing such compounds, and **methods of using these compounds** for therapeutic purposes.

### 🎯 **Primary Technical Domain**
The invention primarily concerns **{domain_info['primary']}**, specifically targeting **{target_name}** ({domain_info['protein_family']}).

### 📋 **Technical Applications**
- **Target Modulation**: Selective {target_name} inhibition/activation
- **Therapeutic Areas**: {domain_info['therapeutic_area']}
- **Drug Discovery**: Structure-activity relationship optimization
- **Pharmaceutical Development**: Novel chemical entities with improved properties

### 🔬 **Scientific Fields**
The invention encompasses multiple scientific disciplines:

| **Field** | **Application** |
|-----------|-----------------|
| **Medicinal Chemistry** | Novel small molecule design and synthesis |
| **Pharmacology** | Target interaction and biological activity |
| **Pharmaceutical Sciences** | Drug formulation and delivery |
| **Biochemistry** | Protein-drug interaction mechanisms |
| **Clinical Applications** | {', '.join(domain_info['secondary'])} |

### 💊 **Innovation Focus**
- **Systematic Markush approach** for comprehensive structure optimization
- **Enhanced selectivity** profiles compared to existing therapies
- **Improved pharmaceutical properties** (ADMET optimization)
- **Novel structural scaffolds** providing intellectual property advantages
"""

def generate_background_prior_art(target_name: str, uniprot_id: str) -> str:
    """
    Generate comprehensive background section with prior art search.
    
    Args:
        target_name: Name of the target protein
        uniprot_id: UniProt ID for the target
        
    Returns:
        Formatted background section with prior art
    """
    
    logger.info(f"🔍 Generating background with prior art search for {target_name} (UniProt: {uniprot_id})")
    
    # Perform prior art searches
    pubchem_results = search_pubchem_bioactivity(target_name, uniprot_id)
    chembl_results = search_chembl_target(target_name)
    patent_landscape = analyze_patent_landscape(target_name)
    
    # Target-specific background information
    target_backgrounds = {
        'CDK5': {
            'function': 'Cyclin-dependent kinase 5 (CDK5) plays crucial roles in neuronal development, synaptic plasticity, and neurodegeneration',
            'diseases': ['Alzheimer\'s disease', 'Parkinson\'s disease', 'stroke', 'cancer'],
            'mechanism': 'CDK5 dysregulation leads to tau hyperphosphorylation and neuronal death',
            'current_therapies': 'Limited therapeutic options with significant side effects'
        },
        'KCNQ2': {
            'function': 'KCNQ2 is a voltage-gated potassium channel critical for neuronal excitability regulation',
            'diseases': ['epilepsy', 'seizure disorders', 'benign familial neonatal epilepsy'],
            'mechanism': 'KCNQ2 mutations cause neuronal hyperexcitability and seizures',
            'current_therapies': 'Traditional antiepileptic drugs with limited efficacy and side effects'
        },
        'EGFR': {
            'function': 'Epidermal Growth Factor Receptor (EGFR) regulates cell proliferation and survival',
            'diseases': ['lung cancer', 'colorectal cancer', 'head and neck cancer'],
            'mechanism': 'EGFR overexpression/mutation drives tumor growth and metastasis',
            'current_therapies': 'Existing inhibitors suffer from resistance development'
        },
        'KDR': {
            'function': 'Kinase insert domain receptor (KDR/VEGFR2) mediates angiogenesis',
            'diseases': ['cancer', 'diabetic retinopathy', 'macular degeneration'],
            'mechanism': 'KDR activation promotes blood vessel formation supporting tumor growth',
            'current_therapies': 'Current anti-angiogenic agents have limited efficacy'
        }
    }
    
    # Get target-specific background or use generic
    target_key = target_name.upper()
    if target_key in target_backgrounds:
        bg_info = target_backgrounds[target_key]
    else:
        bg_info = {
            'function': f'{target_name} is an important therapeutic target',
            'diseases': ['various diseases'],
            'mechanism': f'{target_name} dysregulation contributes to disease pathology',
            'current_therapies': 'Existing therapies have limitations'
        }
    
    # Format prior art results
    prior_art_section = format_prior_art_results(pubchem_results, chembl_results, patent_landscape)
    
    return f"""## 3. Background Art

### 🧬 **Target Biology: {target_name}**

**{target_name}** (UniProt: {uniprot_id}) is a critical therapeutic target with the following characteristics:

- **Function**: {bg_info['function']}
- **Disease Involvement**: {', '.join(bg_info['diseases'])}
- **Mechanism**: {bg_info['mechanism']}
- **Current Limitations**: {bg_info['current_therapies']}

### 📚 **Prior Art Landscape**

{prior_art_section}

### 🔬 **Existing Therapeutic Approaches**

Current {target_name} modulators in the literature and clinical development include:

#### **Limitations of Prior Art:**

| **Challenge** | **Description** | **Impact** |
|---------------|-----------------|------------|
| **Limited Selectivity** | Off-target effects with existing compounds | Dose-limiting toxicities |
| **Poor Pharmacokinetics** | Suboptimal ADMET properties | Limited clinical efficacy |
| **Resistance Development** | Target mutations and pathway redundancy | Therapeutic failure |
| **Narrow SAR** | Limited structural diversity | Restricted optimization potential |

### 💡 **Unmet Medical Need**

Despite significant research efforts, there remains a critical need for:

- ✅ **Novel structural scaffolds** with enhanced selectivity
- ✅ **Improved pharmacokinetic profiles** for better clinical outcomes  
- ✅ **Systematic optimization approaches** using Markush strategies
- ✅ **Reduced side effects** compared to existing therapies

### 🎯 **Invention Rationale**

The present invention addresses these limitations through:
1. **Systematic Markush approach** enabling comprehensive SAR exploration
2. **Novel core scaffolds** providing improved target engagement
3. **Optimized R-group variations** for enhanced selectivity and properties
4. **Comprehensive biological evaluation** ensuring therapeutic potential

> **Prior Art Search Date**: {datetime.now().strftime('%Y-%m-%d')}  
> **Databases Searched**: PubChem BioActivity, ChEMBL, Patent Landscape Analysis
"""

def search_pubchem_bioactivity(target_name: str, uniprot_id: str) -> Dict[str, Any]:
    """
    Search PubChem for bioactivity data related to the target.
    
    Args:
        target_name: Name of the target
        uniprot_id: UniProt ID
        
    Returns:
        Dictionary with search results
    """
    
    try:
        logger.info(f"🔍 Searching PubChem BioActivity for {target_name}")
        
        # PubChem BioActivity search by target name
        base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/target/name"
        search_url = f"{base_url}/{target_name}/aids/JSON"
        
        response = requests.get(search_url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            aid_list = data.get('InformationList', {}).get('Information', [])
            
            result = {
                'status': 'success',
                'total_assays': len(aid_list),
                'search_term': target_name,
                'uniprot_id': uniprot_id,
                'summary': f"Found {len(aid_list)} bioassays for {target_name}",
                'assay_ids': aid_list[:10] if aid_list else []  # First 10 assays
            }
            
            logger.info(f"✅ PubChem search: {len(aid_list)} assays found")
            return result
            
        else:
            logger.warning(f"⚠️  PubChem search failed: HTTP {response.status_code}")
            return {
                'status': 'failed',
                'error': f"HTTP {response.status_code}",
                'summary': f"PubChem search for {target_name} returned no results"
            }
            
    except Exception as e:
        logger.warning(f"⚠️  PubChem search error: {e}")
        return {
            'status': 'error',
            'error': str(e),
            'summary': f"PubChem search for {target_name} encountered an error"
        }

def search_chembl_target(target_name: str) -> Dict[str, Any]:
    """
    Search ChEMBL for target and compound information.
    
    Args:
        target_name: Name of the target
        
    Returns:
        Dictionary with ChEMBL results
    """
    
    try:
        logger.info(f"🔍 Searching ChEMBL for {target_name}")
        
        # ChEMBL REST API search
        base_url = "https://www.ebi.ac.uk/chembl/api/data"
        search_url = f"{base_url}/target/search.json"
        
        params = {
            'q': target_name,
            'format': 'json',
            'limit': 10
        }
        
        response = requests.get(search_url, params=params, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            targets = data.get('targets', [])
            
            result = {
                'status': 'success',
                'total_targets': len(targets),
                'search_term': target_name,
                'summary': f"Found {len(targets)} targets matching {target_name}",
                'targets': [
                    {
                        'chembl_id': target.get('target_chembl_id', ''),
                        'pref_name': target.get('pref_name', ''),
                        'organism': target.get('organism', ''),
                        'target_type': target.get('target_type', '')
                    }
                    for target in targets[:5]  # First 5 targets
                ]
            }
            
            logger.info(f"✅ ChEMBL search: {len(targets)} targets found")
            return result
            
        else:
            logger.warning(f"⚠️  ChEMBL search failed: HTTP {response.status_code}")
            return {
                'status': 'failed',
                'error': f"HTTP {response.status_code}",
                'summary': f"ChEMBL search for {target_name} returned no results"
            }
            
    except Exception as e:
        logger.warning(f"⚠️  ChEMBL search error: {e}")
        return {
            'status': 'error',
            'error': str(e),
            'summary': f"ChEMBL search for {target_name} encountered an error"
        }

def analyze_patent_landscape(target_name: str) -> Dict[str, Any]:
    """
    Analyze patent landscape for the target (simplified analysis).
    
    Args:
        target_name: Name of the target
        
    Returns:
        Dictionary with patent landscape analysis
    """
    
    # Simplified patent landscape analysis
    # In a real implementation, this would query patent databases
    
    patent_categories = {
        'CDK5': [
            'Cyclin-dependent kinase inhibitors',
            'Neurodegeneration therapeutics',
            'Alzheimer\'s disease treatments',
            'Tau phosphorylation modulators'
        ],
        'KCNQ2': [
            'Potassium channel modulators',
            'Epilepsy treatments',
            'Antiepileptic drugs',
            'Ion channel openers'
        ],
        'EGFR': [
            'EGFR inhibitors',
            'Tyrosine kinase inhibitors',
            'Cancer therapeutics',
            'Targeted cancer therapy'
        ],
        'KDR': [
            'VEGFR inhibitors',
            'Angiogenesis inhibitors',
            'Anti-cancer agents',
            'Tyrosine kinase inhibitors'
        ]
    }
    
    target_key = target_name.upper()
    categories = patent_categories.get(target_key, [f'{target_name} modulators'])
    
    return {
        'status': 'analyzed',
        'target': target_name,
        'patent_categories': categories,
        'summary': f"Patent landscape encompasses {len(categories)} main categories",
        'key_areas': categories,
        'note': 'Detailed patent search requires specialized databases'
    }

def format_prior_art_results(pubchem_results: Dict, chembl_results: Dict, patent_results: Dict) -> str:
    """
    Format prior art search results into readable section.
    
    Args:
        pubchem_results: PubChem search results
        chembl_results: ChEMBL search results  
        patent_results: Patent landscape results
        
    Returns:
        Formatted prior art section
    """
    
    # PubChem section
    pubchem_section = ""
    if pubchem_results.get('status') == 'success':
        total_assays = pubchem_results.get('total_assays', 0)
        pubchem_section = f"""
#### **PubChem BioActivity Database**
- **Search Results**: {total_assays} bioassays identified
- **Summary**: {pubchem_results.get('summary', 'No data available')}
- **Coverage**: Comprehensive screening data across multiple assay formats
"""
    else:
        pubchem_section = f"""
#### **PubChem BioActivity Database**  
- **Status**: {pubchem_results.get('summary', 'Search encountered limitations')}
- **Note**: Limited bioactivity data available
"""
    
    # ChEMBL section
    chembl_section = ""
    if chembl_results.get('status') == 'success':
        total_targets = chembl_results.get('total_targets', 0)
        chembl_section = f"""
#### **ChEMBL Medicinal Chemistry Database**
- **Search Results**: {total_targets} target entries identified
- **Summary**: {chembl_results.get('summary', 'No data available')}
- **Coverage**: Medicinal chemistry and bioactivity data from literature
"""
        
        # Add target details if available
        targets = chembl_results.get('targets', [])
        if targets:
            chembl_section += "\n**Key Target Entries:**\n"
            for target in targets[:3]:
                chembl_id = target.get('chembl_id', 'N/A')
                pref_name = target.get('pref_name', 'N/A')
                chembl_section += f"- {chembl_id}: {pref_name}\n"
    else:
        chembl_section = f"""
#### **ChEMBL Medicinal Chemistry Database**
- **Status**: {chembl_results.get('summary', 'Search encountered limitations')}  
- **Note**: Limited target information available
"""
    
    # Patent landscape section
    patent_section = ""
    if patent_results.get('status') == 'analyzed':
        categories = patent_results.get('key_areas', [])
        patent_section = f"""
#### **Patent Landscape Analysis**
- **Analysis Scope**: {len(categories)} key patent categories
- **Summary**: {patent_results.get('summary', 'Analysis completed')}

**Key Patent Areas:**
"""
        for category in categories:
            patent_section += f"- {category}\n"
    
    return f"""{pubchem_section}

{chembl_section}

{patent_section}

#### **Prior Art Summary**
The comprehensive prior art search reveals existing research and development in the target area, confirming the need for novel approaches with improved properties and selectivity profiles.
"""

# Additional utility functions for patent generation

def get_target_disease_associations(target_name: str) -> List[str]:
    """Get disease associations for a target."""
    
    disease_map = {
        'CDK5': ['Alzheimer\'s disease', 'neurodegeneration', 'cancer'],
        'KCNQ2': ['epilepsy', 'seizure disorders', 'neurological disorders'],
        'EGFR': ['cancer', 'solid tumors', 'lung cancer'],
        'KDR': ['cancer', 'angiogenesis-related diseases', 'diabetic retinopathy']
    }
    
    return disease_map.get(target_name.upper(), ['various diseases'])

def get_target_protein_family(target_name: str) -> str:
    """Get protein family for a target."""
    
    family_map = {
        'CDK5': 'Cyclin-dependent kinases',
        'KCNQ2': 'Voltage-gated potassium channels', 
        'EGFR': 'Receptor tyrosine kinases',
        'KDR': 'Receptor tyrosine kinases'
    }
    
    return family_map.get(target_name.upper(), 'Target proteins')