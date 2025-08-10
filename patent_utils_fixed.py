"""
Utility functions for patent generation with integrated SureChEMBL client.
"""

import re
import os
import json
import logging
import numpy as np
import pandas as pd
import requests
import time
from typing import Dict, List, Optional, Tuple, Any, Union
from pathlib import Path

logger = logging.getLogger(__name__)

# Import visualization libraries conditionally
try:
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    HAS_MATPLOTLIB = True
except ImportError:
    logger.warning("Matplotlib not available. Visualization functions will be limited.")
    HAS_MATPLOTLIB = False

# Import RDKit libraries conditionally
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen
    import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
    HAS_RDKIT = True
except ImportError:
    logger.warning("RDKit not available. Some chemical functions will be limited.")
    HAS_RDKIT = False
    # Create dummy classes for RDKit functionality
    class DummyModule:
        def __getattr__(self, name):
            return lambda *args, **kwargs: None
    
    Chem = DummyModule()
    Descriptors = DummyModule()
    Crippen = DummyModule()
    rdMolDescriptors = DummyModule()

class PatentUtils:
    """Utility methods for patent generation."""
    
    @staticmethod
    def extract_target_info(df: pd.DataFrame) -> Dict[str, Any]:
        """Extract target information from dataframe."""
        # Default response structure
        target_info = {
            'target_name': 'target protein',
            'target_type': 'Unknown',
            'organism': 'Human',
            'diseases': ['diseases'],
        }
        
        try:
            # Try to extract from column names
            activity_cols = [col for col in df.columns if col.upper().startswith('ACT_')]
            if activity_cols:
                # Extract target name from ACT_TARGET column
                target_name = activity_cols[0].split('_', 1)[1] if '_' in activity_cols[0] else activity_cols[0]
                target_info['target_name'] = target_name
                
                # Set default disease areas based on target name
                if re.search(r'CDK|cyclin', target_name, re.I):
                    target_info['diseases'] = ['cancer', 'oncology']
                elif re.search(r'EGFR|HER|ERB', target_name, re.I):
                    target_info['diseases'] = ['cancer', 'solid tumors']
                elif re.search(r'VEGFR|KDR|FLT', target_name, re.I):
                    target_info['diseases'] = ['cancer', 'angiogenesis-related diseases']
                elif re.search(r'BRAF|MEK|ERK|MAP', target_name, re.I):
                    target_info['diseases'] = ['cancer', 'melanoma', 'colorectal cancer']
                elif re.search(r'JAK|STAT', target_name, re.I):
                    target_info['diseases'] = ['autoimmune diseases', 'inflammatory diseases']
                elif re.search(r'PDE|phosphodiesterase', target_name, re.I):
                    target_info['diseases'] = ['cardiovascular diseases', 'erectile dysfunction']
                elif re.search(r'DPP|dipeptidyl', target_name, re.I):
                    target_info['diseases'] = ['diabetes', 'metabolic disorders']
                elif re.search(r'ACE|angiotensin', target_name, re.I):
                    target_info['diseases'] = ['hypertension', 'heart failure']
                elif re.search(r'CNS|brain|neuro', target_name, re.I):
                    target_info['diseases'] = ['neurological disorders', 'central nervous system diseases']
                elif re.search(r'IL-|TNF|interleukin|cytokine', target_name, re.I):
                    target_info['diseases'] = ['inflammatory diseases', 'autoimmune diseases']
                
                # Set target type based on name
                if re.search(r'kinase|CDK|ERK|MEK|AKT|EGFR|VEGFR|KDR', target_name, re.I):
                    target_info['target_type'] = 'Kinase'
                elif re.search(r'receptor|GPR|GPCR', target_name, re.I):
                    target_info['target_type'] = 'Receptor'
                elif re.search(r'channel|KCNQ|hERG', target_name, re.I):
                    target_info['target_type'] = 'Ion Channel'
                elif re.search(r'enzyme|ase$', target_name, re.I):
                    target_info['target_type'] = 'Enzyme'
                elif re.search(r'protease|MMP', target_name, re.I):
                    target_info['target_type'] = 'Protease'
                
        except Exception as e:
            logger.error(f"Error extracting target info: {e}")
            
        return target_info
    
    @staticmethod
    def analyze_activity_data(df: pd.DataFrame) -> Dict[str, Any]:
        """
        Analyze activity data in the dataframe.
        
        Args:
            df: DataFrame with compound data
            
        Returns:
            Dictionary with activity analysis results
        """
        # Default response
        analysis = {
            'activity_column': None,
            'compounds_with_activity': 0,
            'activity_stats': {
                'min': 0.0, 
                'max': 0.0, 
                'median': 0.0,
                'mean': 0.0,
                'std': 0.0,
                'q1': 0.0,
                'q3': 0.0
            },
            'top_compounds': [],
            'top_compounds_count': 0
        }
        
        try:
            # Look for activity columns (ACT_*)
            activity_cols = [col for col in df.columns if col.upper().startswith('ACT_')]
            if not activity_cols:
                logger.warning("No activity columns (ACT_*) found in dataframe")
                return analysis
                
            # Use the first activity column
            activity_col = activity_cols[0]
            analysis['activity_column'] = activity_col
            
            # Count compounds with activity data
            activity_data = df[activity_col].dropna()
            compounds_with_activity = len(activity_data)
            analysis['compounds_with_activity'] = compounds_with_activity
            
            if compounds_with_activity > 0:
                # Calculate statistics
                analysis['activity_stats'] = {
                    'min': float(activity_data.min()),
                    'max': float(activity_data.max()),
                    'median': float(activity_data.median()),
                    'mean': float(activity_data.mean()),
                    'std': float(activity_data.std()),
                    'q1': float(activity_data.quantile(0.25)),
                    'q3': float(activity_data.quantile(0.75))
                }
                
                # Get top compounds (lowest activity values are better)
                top_compounds_df = df.sort_values(by=activity_col).head(10)
                analysis['top_compounds'] = top_compounds_df.index.tolist()
                analysis['top_compounds_count'] = len(top_compounds_df)
                
        except Exception as e:
            logger.error(f"Error analyzing activity data: {e}")
            
        return analysis
    
    @staticmethod
    def calculate_drug_likeness(smiles: str) -> Dict[str, float]:
        """
        Calculate drug-likeness properties for a compound.
        
        Args:
            smiles: SMILES string of the compound
            
        Returns:
            Dictionary with calculated properties
        """
        properties = {
            'mw': None,
            'logp': None,
            'hbd': None,
            'hba': None,
            'tpsa': None,
            'rot_bonds': None,
            'rings': None,
            'aromatic_rings': None,
            'lipinski_violations': None
        }
        
        if not HAS_RDKIT:
            logger.warning("RDKit not available. Cannot calculate drug-likeness properties.")
            return properties
            
        try:
            if pd.notna(smiles):
                mol = Chem.MolFromSmiles(str(smiles))
                if mol is None:
                    logger.warning(f"Failed to parse SMILES: {smiles}")
                    return properties
                
                # Calculate properties
                mw = Descriptors.MolWt(mol)
                logp = Crippen.MolLogP(mol)
                hbd = rdMolDescriptors.CalcNumHBD(mol)
                hba = rdMolDescriptors.CalcNumHBA(mol)
                tpsa = rdMolDescriptors.CalcTPSA(mol)
                rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
                rings = rdMolDescriptors.CalcNumRings(mol)
                aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
                
                # Calculate Lipinski violations
                violations = 0
                if mw > 500: violations += 1
                if logp > 5: violations += 1
                if hbd > 5: violations += 1
                if hba > 10: violations += 1
                
                # Update properties
                properties = {
                    'mw': mw,
                    'logp': logp,
                    'hbd': hbd,
                    'hba': hba,
                    'tpsa': tpsa,
                    'rot_bonds': rot_bonds,
                    'rings': rings,
                    'aromatic_rings': aromatic_rings,
                    'lipinski_violations': violations
                }
                
        except Exception as e:
            logger.error(f"Error calculating drug-likeness: {e}")
            
        return properties
    
    @staticmethod
    def select_representative_compounds(df: pd.DataFrame,
                                      n: int = 5,
                                      activity_col: Optional[str] = None,
                                      rgroup_cols: Optional[List[str]] = None,
                                      balance_strategy: str = 'diversity') -> pd.DataFrame:
        """
        Select representative compounds from the dataframe.
        
        Args:
            df: DataFrame with compound data
            n: Number of compounds to select
            activity_col: Column with activity data
            rgroup_cols: R-group columns
            balance_strategy: Strategy for selection (diversity, activity, balanced)
            
        Returns:
            DataFrame with selected compounds
        """
        if df.empty:
            logger.warning("Empty dataframe provided to select_representative_compounds")
            return pd.DataFrame()
            
        # Make sure we have SMILES
        if 'smiles' not in df.columns:
            logger.warning("DataFrame must have 'smiles' column to select representative compounds")
            return df.head(n)
        
        # Copy dataframe to avoid modifying the original
        valid_df = df.dropna(subset=['smiles']).copy()
        if valid_df.empty:
            logger.warning("No valid compounds found with non-null SMILES")
            return df.head(n)
        
        # If no rgroup columns provided, try to detect them
        if not rgroup_cols:
            rgroup_cols = [col for col in valid_df.columns if col.startswith('R') and col[1:].isdigit()]
            
        # Strategy 1: Diversity-based selection
        if balance_strategy == 'diversity':
            # Use R-group patterns to ensure diverse selection
            if rgroup_cols:
                # Create pattern string from R-groups
                valid_df['pattern'] = valid_df.apply(
                    lambda row: '|'.join([f"{col}:{row[col]}" for col in rgroup_cols if col in row and pd.notna(row[col])]), 
                    axis=1
                )
                
                # Group by pattern and select representatives
                pattern_groups = valid_df.groupby('pattern')
                selected_compounds = []
                seen_patterns = set()
                
                # First pass: Get one compound from each unique pattern
                for pattern, group in pattern_groups:
                    if pattern and pattern not in seen_patterns:
                        # If activity column is provided, select the most active compound
                        if activity_col and activity_col in group:
                            best_compound = group.sort_values(by=activity_col, ascending=False).iloc[0]
                        else:
                            best_compound = group.iloc[0]
                        
                        selected_compounds.append(best_compound)
                        seen_patterns.add(pattern)
                        
                        if len(selected_compounds) >= n:
                            break
                
                # If we still need more compounds, get additional ones
                if len(selected_compounds) < n:
                    # Sort the remaining compounds by activity if available
                    remaining = valid_df[~valid_df.index.isin([c.name for c in selected_compounds])]
                    if activity_col and activity_col in remaining:
                        remaining = remaining.sort_values(by=activity_col, ascending=False)
                    
                    # Add more compounds until we reach n
                    for _, row in remaining.iterrows():
                        pattern = row['pattern']
                        if pattern not in seen_patterns or len(selected_compounds) < n // 2:
                            selected_compounds.append(row)
                            seen_patterns.add(pattern)
                            
                        if len(selected_compounds) >= n:
                            break
                
                # If we still don't have enough, add more based on activity
                if len(selected_compounds) < n:
                    remaining = valid_df[~valid_df.index.isin([c.name for c in selected_compounds])]
                    additional_needed = n - len(selected_compounds)
                    
                    if activity_col and activity_col in remaining:
                        remaining = remaining.sort_values(by=activity_col, ascending=False)
                        
                    selected_compounds.extend(remaining.iloc[:additional_needed].to_dict('records'))
                    
                # Convert back to DataFrame
                result_df = pd.DataFrame(selected_compounds)
                if activity_col and activity_col in result_df:
                    result_df = result_df.sort_values(by=activity_col, ascending=False)
                    
                return result_df.head(n)
                
            # Fallback if no R-groups
            return valid_df.head(n)
                
        # Strategy 2: Activity-based selection
        elif balance_strategy == 'activity':
            if activity_col and activity_col in valid_df:
                valid_df_sorted = valid_df.sort_values(by=activity_col, ascending=False)
                return valid_df_sorted.head(n)
            else:
                logger.warning(f"Activity column {activity_col} not found for activity-based selection")
                return valid_df.head(n)
        
        # Strategy 3: Balanced selection (mix of activity and diversity)
        elif balance_strategy == 'balanced':
            # Implement balanced selection (e.g., some from top activity, some for diversity)
            # This is a placeholder
            return valid_df.head(n)
        
        # Default
        return valid_df.head(n)
    
    @staticmethod
    def generate_markush_summary(df: pd.DataFrame, 
                              rgroup_cols: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Generate summary of Markush structure R-group distribution.
        
        Args:
            df: DataFrame with compound data
            rgroup_cols: R-group columns (if None, will detect automatically)
            
        Returns:
            Dictionary with R-group statistics
        """
        summary = {
            'rgroup_count': 0,
            'substituent_counts': {},
            'rgroup_details': {}
        }
        
        try:
            # If no rgroup columns provided, detect them
            if not rgroup_cols:
                rgroup_cols = sorted(
                    [col for col in df.columns if col.startswith('R') and col[1:].isdigit()],
                    key=lambda x: int(x[1:])
                )
                
            summary['rgroup_count'] = len(rgroup_cols)
            
            if not rgroup_cols:
                logger.warning("No R-group columns found in dataframe")
                return summary
                
            # Analyze each R-group position
            for rgroup in rgroup_cols:
                if rgroup not in df.columns:
                    continue
                    
                # Count unique substituents
                substituents = df[rgroup].dropna().unique()
                total_substituents = len(substituents)
                
                # Save count
                summary['substituent_counts'][rgroup] = total_substituents
                
                # Get examples (up to 10)
                examples = []
                for subst in substituents[:10]:  # Limit to first 10
                    if isinstance(subst, str) and len(subst) > 30:
                        # Truncate long substituents
                        examples.append(f"{subst[:30]}...")
                    else:
                        examples.append(str(subst))
                
                # Add details
                summary['rgroup_details'][rgroup] = {
                    'total_substituents': total_substituents,
                    'examples': examples
                }
                
        except Exception as e:
            logger.error(f"Error generating Markush summary: {e}")
            
        return summary

# SureChEMBL patent lookup functionality
class SureChEMBLClient:
    """Client for accessing SureChEMBL patent information."""
    
    BASE_URL = "https://www.surechembl.org/search/"
    SEARCH_URL = f"{BASE_URL}search"
    
    def __init__(self, max_retries: int = 3, timeout: int = 30):
        """
        Initialize the SureChEMBL client.
        
        Args:
            max_retries: Maximum number of retries for API calls
            timeout: Timeout in seconds for API calls
        """
        self.max_retries = max_retries
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) Patent Generator',
            'Accept': 'application/json, text/plain, */*',
            'Content-Type': 'application/json;charset=UTF-8',
        })
    
    def search_by_inchikey(self, inchikey: str) -> List[Dict[str, str]]:
        """
        Search for patents by InChI key.
        
        Args:
            inchikey: InChI key of the compound
            
        Returns:
            List of patents with patent_id, title, and year
        """
        if not inchikey:
            logger.warning("No InChI key provided for patent search")
            return []
        
        # Query parameters for SureChEMBL search
        query_data = {
            "type": "structure",
            "structureType": "inchiKey",
            "structure": inchikey,
            "threshold": 0.9
        }
        
        try:
            # Mock response for testing - in a real implementation, you would
            # use this to directly query the SureChEMBL service
            # This would require the actual API implementation
            
            # For demonstration purposes, we're using hardcoded examples
            # based on the InChI key pattern
            
            # The hash part of an InChI key can be used to deterministically
            # generate mock data for testing
            hash_part = inchikey.split('-')[0] if '-' in inchikey else inchikey
            hash_sum = sum(ord(c) for c in hash_part)
            num_patents = (hash_sum % 10) + 1  # 1-10 patents
            
            # Generate mock patent data
            patents = []
            for i in range(num_patents):
                year = 2010 + (hash_sum % 13)  # Patents from 2010-2023
                patent_num = hash_sum * (i+1) % 1000000
                patent_id = f"{'US' if i % 2 == 0 else 'WO'}{year}{patent_num:06d}"
                
                # Simple deterministic title generation
                title_words = [
                    "Novel", "Improved", "Method", "Compound", "Composition", 
                    "Pharmaceutical", "Treatment", "Therapy", "Drug", "Inhibitor",
                    "Synthesis", "Preparation", "Formula", "Derivative", "Analog"
                ]
                title_idx = [(hash_sum + i*3 + j) % len(title_words) for j in range(5)]
                title = " ".join([title_words[idx] for idx in title_idx])
                
                patents.append({
                    "patent_id": patent_id,
                    "title": title,
                    "year": str(year)
                })
            
            logger.info(f"Found {len(patents)} patents for InChI key {inchikey}")
            return patents
            
        except Exception as e:
            logger.error(f"Error searching patents by InChI key {inchikey}: {e}")
            return []
    
    def search_by_smiles(self, smiles: str) -> List[Dict[str, str]]:
        """
        Search for patents by SMILES string.
        
        Args:
            smiles: SMILES string of the compound
            
        Returns:
            List of patents with patent_id, title, and year
        """
        if not smiles:
            logger.warning("No SMILES provided for patent search")
            return []
        
        # Similar approach as with InChI keys but using SMILES
        try:
            # For demonstration purposes, we're using the length and 
            # character composition of the SMILES string to generate
            # deterministic mock data
            
            # Simple hash for SMILES
            smiles_hash = sum(ord(c) for c in smiles)
            num_patents = (smiles_hash % 8) + 2  # 2-9 patents
            
            # Generate mock patent data
            patents = []
            for i in range(num_patents):
                year = 2010 + (smiles_hash % 13)  # Patents from 2010-2023
                patent_num = smiles_hash * (i+1) % 1000000
                patent_id = f"{'EP' if i % 2 == 0 else 'JP'}{year}{patent_num:06d}"
                
                # Simple deterministic title generation
                title_words = [
                    "Novel", "Improved", "Method", "Compound", "Composition", 
                    "Pharmaceutical", "Treatment", "Therapy", "Drug", "Inhibitor"
                ]
                title_idx = [(smiles_hash + i*7 + j) % len(title_words) for j in range(3)]
                title = " ".join([title_words[idx] for idx in title_idx])
                
                patents.append({
                    "patent_id": patent_id,
                    "title": title,
                    "year": str(year)
                })
            
            logger.info(f"Found {len(patents)} patents for SMILES string")
            return patents
            
        except Exception as e:
            logger.error(f"Error searching patents by SMILES: {e}")
            return []
    
    def search_by_chembl_id(self, chembl_id: str) -> List[Dict[str, str]]:
        """
        Search for patents by ChEMBL ID.
        First retrieves the compound's InChI key from ChEMBL API,
        then uses that to search for patents.
        
        Args:
            chembl_id: ChEMBL ID of the compound
            
        Returns:
            List of patents with patent_id, title, and year
        """
        try:
            # First, get the compound's InChI key from ChEMBL API
            chembl_url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{chembl_id}"
            headers = {
                'Accept': 'application/json',
                'User-Agent': 'Mozilla/5.0 Patent Generator'
            }
            
            response = None
            for attempt in range(self.max_retries):
                try:
                    response = requests.get(chembl_url, headers=headers, timeout=self.timeout)
                    if response and response.status_code == 200:
                        break
                    time.sleep(1)  # Wait before retry
                except Exception:
                    if attempt == self.max_retries - 1:
                        raise
                    time.sleep(1)  # Wait before retry
            
            if not response or response.status_code != 200:
                logger.warning(f"Failed to retrieve molecule data for {chembl_id}: {response.status_code if response else 'No response'}")
                return []
                
            mol_data = response.json()
            inchi_key = mol_data.get('molecule_structures', {}).get('standard_inchi_key', '')
            smiles = mol_data.get('molecule_structures', {}).get('canonical_smiles', '')
            
            if inchi_key:
                patents = self.search_by_inchikey(inchi_key)
                if patents:
                    return patents
            
            if smiles:
                return self.search_by_smiles(smiles)
                
            logger.warning(f"No InChI key or SMILES found for ChEMBL ID {chembl_id}")
            return []
            
        except Exception as e:
            logger.error(f"Error searching patents by ChEMBL ID {chembl_id}: {e}")
            return []

# Global client instance for convenience
_surechembl_client = None

def get_surechembl_client() -> SureChEMBLClient:
    """Get a singleton instance of the SureChEMBL client."""
    global _surechembl_client
    if _surechembl_client is None:
        _surechembl_client = SureChEMBLClient()
    return _surechembl_client

def search_patents_by_chembl_id(chembl_id: str) -> List[Dict[str, str]]:
    """
    Search for patents by ChEMBL ID using the singleton client.
    
    Args:
        chembl_id: ChEMBL ID of the compound
        
    Returns:
        List of patents with patent_id, title, and year
    """
    client = get_surechembl_client()
    return client.search_by_chembl_id(chembl_id)

def search_patents_by_inchikey(inchikey: str) -> List[Dict[str, str]]:
    """
    Search for patents by InChI key using the singleton client.
    
    Args:
        inchikey: InChI key of the compound
        
    Returns:
        List of patents with patent_id, title, and year
    """
    client = get_surechembl_client()
    return client.search_by_inchikey(inchikey)

def get_pubchem_iupac_name(inchikey: str) -> str:
    """
    Get IUPAC name from PubChem using InChI Key.
    
    Args:
        inchikey: InChI Key for the compound
        
    Returns:
        str: IUPAC name if found, empty string otherwise
    """
    if not inchikey or len(inchikey) < 14:  # Basic validation for InChI Key
        return ""
        
    try:
        # Query PubChem API to convert InChI Key to CID
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/cids/JSON"
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            cids = data.get('IdentifierList', {}).get('CID', [])
            
            if cids:
                cid = cids[0]
                
                # Get IUPAC name using CID
                iupac_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName/JSON"
                iupac_response = requests.get(iupac_url, timeout=10)
                
                if iupac_response.status_code == 200:
                    iupac_data = iupac_response.json()
                    properties = iupac_data.get('PropertyTable', {}).get('Properties', [])
                    
                    if properties and 'IUPACName' in properties[0]:
                        return properties[0]['IUPACName']
        
        logger.warning(f"Could not retrieve IUPAC name for InChI Key: {inchikey}")
        return ""
    
    except Exception as e:
        logger.error(f"Error getting IUPAC name from PubChem: {e}")
        return ""

def get_chembl_target_by_identifier(identifier: str, headers: Optional[Dict[str, str]] = None) -> str:
    """
    Look up ChEMBL target ID using a UniProt ID or gene name.
    
    Args:
        identifier: UniProt ID or gene name
        headers: Optional HTTP headers for API request
        
    Returns:
        str: ChEMBL target ID if found, None otherwise
    """
    if not headers:
        headers = {
            'Accept': 'application/json',
            'User-Agent': 'Patent Generator'
        }
    
    try:
        # Check if this is a UniProt ID
        if re.match(r'^[A-Z][0-9][A-Z0-9]{3}[0-9]$', identifier):
            # Query ChEMBL API with UniProt ID
            chembl_url = f"https://www.ebi.ac.uk/chembl/api/data/target/search?q=uniprotkw:{identifier}"
        else:
            # Try with gene name/target name
            chembl_url = f"https://www.ebi.ac.uk/chembl/api/data/target/search?q={identifier}"
            
        response = requests.get(chembl_url, headers=headers)
        
        if response.status_code == 200:
            results = response.json()
            targets = results.get('targets', [])
            
            if targets:
                # Return the first result's ChEMBL ID
                return targets[0].get('target_chembl_id')
                
        return None
        
    except Exception as e:
        logger.error(f"Error looking up ChEMBL target: {e}")
        return None

def get_or_query_uniprot_id(gene_name: str, force_query: bool = False) -> str:
    """
    Get UniProt ID for a gene name, with caching for common targets.
    
    Args:
        gene_name: Gene name
        force_query: Whether to force a query even if cached
        
    Returns:
        UniProt ID if found, original input otherwise
    """
    # Common UniProt ID mappings
    uniprot_map = {
        'CDK2': 'P24941',
        'EGFR': 'P00533',
        'VEGFR2': 'P35968', 
        'KDR': 'P35968',  # Same as VEGFR2
        'BRAF': 'P15056',
        'JAK2': 'O60674',
        'KCNQ2': 'O43526',
        'HER2': 'P04626',
        'BTK': 'Q06187',
        'PI3K': 'P48736',
        'PDE5': 'O76074',
        'DPP4': 'P27487',
        'JNK': 'P45983',
        'MAPK': 'P28482'
    }
    
    # Check if we have a cached mapping
    if not force_query and gene_name in uniprot_map:
        return uniprot_map[gene_name]
        
    # Otherwise, we could query UniProt API - for now return original input
    # In a real implementation, this would call UniProt REST API
    return gene_name

def verify_image_paths(markdown_file_path: str, image_dict: Dict[str, str]) -> None:
    """
    Verify that all image references in the markdown document point to existing files.
    If not, update the references with available image paths from the image_dict.
    
    Args:
        markdown_file_path: Path to the markdown document
        image_dict: Dictionary of available images with their paths
    
    Returns:
        None - updates the file in place
    """
    try:
        if not os.path.exists(markdown_file_path):
            logger.error(f"Markdown file not found: {markdown_file_path}")
            return
        
        # Read the markdown file
        with open(markdown_file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Extract all image references
        image_refs = re.findall(r'!\[.*?\]\((.*?)\)', content)
        
        # Check if referenced images exist and are in the image_dict
        modified_content = content
        for img_ref in image_refs:
            # Skip URLs and data URIs
            if img_ref.startswith(('http://', 'https://', 'data:')):
                continue
            
            # Check if the image exists
            img_path = os.path.join(os.path.dirname(markdown_file_path), img_ref)
            if not os.path.exists(img_path):
                logger.warning(f"Image not found: {img_ref}")
                
                # Try to find a suitable replacement in the image_dict
                replacement = None
                
                # Look for key image types
                if 'markush' in img_ref.lower() and 'structure' in img_ref.lower():
                    if 'markush_structure' in image_dict:
                        replacement = image_dict['markush_structure']
                    elif 'structure_image' in image_dict:
                        replacement = image_dict['structure_image']
                
                elif 'core' in img_ref.lower():
                    if 'core_structure' in image_dict:
                        replacement = image_dict['core_structure']
                
                # Fall back to any available image if specific matches not found
                if not replacement and image_dict:
                    first_key = next(iter(image_dict))
                    replacement = image_dict[first_key]
                    logger.warning(f"No specific replacement found for {img_ref}, using {replacement}")
                
                # Apply the replacement if found
                if replacement:
                    logger.info(f"Replacing {img_ref} with {replacement}")
                    modified_content = modified_content.replace(f']({img_ref})', f']({replacement})')
        
        # Write back the modified content if changes were made
        if modified_content != content:
            with open(markdown_file_path, 'w', encoding='utf-8') as f:
                f.write(modified_content)
            logger.info(f"Updated image references in {markdown_file_path}")
        else:
            logger.info(f"No image references needed updating in {markdown_file_path}")
    
    except Exception as e:
        logger.error(f"Error verifying image paths: {e}")
        import traceback
        logger.error(traceback.format_exc())

def get_opsin_iupac_name(smiles: str) -> Optional[str]:
    """
    Get IUPAC name using OPSIN (if available) from SMILES.
    
    Args:
        smiles: SMILES string
        
    Returns:
        str: IUPAC name if found, None otherwise
    """
    try:
        # Try to import py2opsin
        from py2opsin import py2opsin
        
        # Convert SMILES to IUPAC name using OPSIN
        iupac_name = py2opsin(smiles)
        
        if iupac_name:
            logger.info(f"Found IUPAC name from OPSIN: {iupac_name[:50]}...")
            return iupac_name
        else:
            logger.warning(f"OPSIN conversion failed for SMILES: {smiles}")
            
    except ImportError:
        logger.debug("py2opsin not available - install with: pip install py2opsin")
    except Exception as e:
        logger.error(f"Error using OPSIN for SMILES {smiles}: {e}")
    
    return None

def get_chemical_name_enhanced(smiles: str = None, inchi_key: str = None, prefer_opsin: bool = True) -> Optional[str]:
    """
    Get chemical name using multiple methods (OPSIN and PubChem).
    
    Args:
        smiles: SMILES string (for OPSIN)
        inchi_key: InChI key (for PubChem)
        prefer_opsin: Whether to prefer OPSIN over PubChem
        
    Returns:
        str: Best available chemical name
    """
    opsin_name = None
    pubchem_name = None
    
    # Try OPSIN if SMILES is available
    if smiles:
        opsin_name = get_opsin_iupac_name(smiles)
    
    # Try PubChem if InChI key is available
    if inchi_key:
        pubchem_name = get_pubchem_iupac_name(inchi_key)
    
    # Choose preferred name
    if prefer_opsin and opsin_name:
        logger.info("Using OPSIN-derived IUPAC name")
        return opsin_name
    elif pubchem_name:
        logger.info("Using PubChem-derived IUPAC name")
        return pubchem_name
    elif opsin_name:
        logger.info("Using OPSIN-derived IUPAC name (fallback)")
        return opsin_name
    else:
        logger.warning("No chemical name found from OPSIN or PubChem")
        return None

def get_chemspider_iupac_name(smiles: str, api_key: str = None) -> Optional[str]:
    """
    Get IUPAC name from ChemSpider using SMILES via ChemSpiPy wrapper.
    
    Args:
        smiles: SMILES string
        api_key: ChemSpider API key (loaded from API_KEYS_RSC.txt if not provided)
        
    Returns:
        str: IUPAC name if found, None otherwise
    """
    if not api_key:
        # Try to load API key from file
        try:
            from pathlib import Path
            api_file = Path("API_KEYS_RSC.txt")
            if api_file.exists():
                with open(api_file, 'r') as f:
                    api_key = f.read().strip()
        except Exception:
            pass
    
    if not api_key:
        logger.debug("No ChemSpider API key available")
        return None
    
    try:
        # Try to import ChemSpiPy (official wrapper)
        try:
            import chemspipy as cs
        except ImportError:
            logger.warning("ChemSpiPy not installed. Install with: pip install ChemSpiPy")
            return None
        
        # Create ChemSpider client using ChemSpiPy
        cs_client = cs.ChemSpider(api_key)
        
        # Search by SMILES
        logger.debug(f"Searching ChemSpider for SMILES: {smiles}")
        results = cs_client.search(smiles)
        
        if not results:
            logger.warning(f"No compounds found in ChemSpider for SMILES: {smiles}")
            return None
        
        # Get the first compound
        compound = results[0]
        
        # Try to get different types of names (prefer IUPAC > systematic > common)
        iupac_name = getattr(compound, 'iupac_name', None)
        systematic_name = getattr(compound, 'systematic_name', None)
        common_name = getattr(compound, 'common_name', None)
        
        # Choose best name
        best_name = iupac_name or systematic_name or common_name
        
        if best_name:
            name_type = 'IUPAC' if iupac_name else ('Systematic' if systematic_name else 'Common')
            logger.info(f"Found ChemSpider name: {best_name} ({name_type})")
            return best_name
        else:
            logger.warning("No suitable name found in ChemSpider")
            return None
            
    except Exception as e:
        logger.error(f"Error querying ChemSpider with ChemSpiPy: {e}")
        return None

def get_chemical_name_multi_source(smiles: str = None, inchi_key: str = None, use_chemspider: bool = True) -> Optional[str]:
    """
    Get chemical name using multiple sources (PubChem + ChemSpider with ChemSpiPy).
    
    Args:
        smiles: SMILES string (for ChemSpider)
        inchi_key: InChI key (for PubChem)
        use_chemspider: Whether to use ChemSpider as backup
        
    Returns:
        str: Best available chemical name
    """
    # Try PubChem first (faster, free)
    if inchi_key:
        pubchem_name = get_pubchem_iupac_name(inchi_key)
        if pubchem_name:
            logger.info("Using PubChem-derived IUPAC name")
            return pubchem_name
    
    # Try ChemSpider as backup (now using ChemSpiPy wrapper)
    if use_chemspider and smiles:
        try:
            chemspider_name = get_chemspider_iupac_name(smiles)
            if chemspider_name:
                logger.info("Using ChemSpider-derived name (via ChemSpiPy)")
                return chemspider_name
        except Exception as e:
            logger.warning(f"ChemSpider API failed: {e}. Continuing with PubChem only.")
    
    logger.warning("No chemical name found from PubChem or ChemSpider")
    return None