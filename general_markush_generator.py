"""
Enhanced Markush Structure Generator with Proper R-Group Visualization

This module provides a general solution for generating proper Markush structures
with visible R-group attachment points for any input CSV file.

Key improvements:
1. Proper R-group attachment point detection and labeling
2. Multiple fallback methods for structure visualization
3. Template-based approach for complex SMARTS patterns
4. General solution working with any compound dataset
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw, rdFMCS, rdDepictor
from rdkit.Chem.rdRGroupDecomposition import RGroupDecomposition
import os
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Any
import logging

logger = logging.getLogger(__name__)

class GeneralMarkushGenerator:
    """
    General Markush structure generator that works with any compound dataset.
    
    Solves the "?" atom problem by using multiple visualization strategies:
    1. Template-based rendering with R-group highlights
    2. Core structure with manual R-group labeling
    3. Representative compound approach
    4. Simplified SMARTS conversion
    """
    
    def __init__(self):
        self.valid_mols = []
        self.mcs_result = None
        self.mcs_mol = None
        self.r_group_positions = []
        self.attachment_points = []
        
    def analyze_compound_set(self, smiles_list: List[str]) -> Dict[str, Any]:
        """
        Analyze a set of compounds to identify the common core and R-group positions.
        
        Args:
            smiles_list: List of SMILES strings
            
        Returns:
            Analysis results with MCS and R-group information
        """
        logger.info(f"Analyzing {len(smiles_list)} compounds for Markush structure")
        
        # Convert SMILES to molecules with proper sanitization
        self.valid_mols = []
        for i, smiles in enumerate(smiles_list):
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # Sanitize the molecule to ensure proper valence calculation
                    try:
                        Chem.SanitizeMol(mol)
                        # Ensure explicit valence is calculated
                        for atom in mol.GetAtoms():
                            atom.CalcExplicitValence()
                        self.valid_mols.append(mol)
                    except Exception as e:
                        logger.warning(f"Sanitization failed for SMILES at position {i}: {smiles} - {e}")
                        # Try without full sanitization
                        try:
                            mol = Chem.MolFromSmiles(smiles, sanitize=False)
                            if mol:
                                Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_PROPERTIES)
                                self.valid_mols.append(mol)
                        except:
                            logger.warning(f"Complete sanitization failed for SMILES at position {i}: {smiles}")
                else:
                    logger.warning(f"Invalid SMILES at position {i}: {smiles}")
            except Exception as e:
                logger.warning(f"Error processing SMILES at position {i}: {smiles} - {e}")
        
        logger.info(f"Successfully parsed {len(self.valid_mols)} valid molecules")
        
        if len(self.valid_mols) < 2:
            logger.error("Need at least 2 valid molecules for MCS analysis")
            return {"error": "Insufficient valid molecules"}
        
        # Find maximum common substructure with multiple approaches
        self.mcs_result = self._find_robust_mcs(self.valid_mols)
        
        if self.mcs_result and self.mcs_result.numAtoms > 0:
            try:
                self.mcs_mol = Chem.MolFromSmarts(self.mcs_result.smartsString)
                if self.mcs_mol:
                    # Ensure the MCS molecule is properly sanitized
                    try:
                        Chem.SanitizeMol(self.mcs_mol)
                    except:
                        # If sanitization fails, try alternative approach
                        logger.warning("MCS molecule sanitization failed, using alternative approach")
                logger.info(f"MCS found: {self.mcs_result.numAtoms} atoms, {self.mcs_result.numBonds} bonds")
            except Exception as e:
                logger.warning(f"Error creating MCS molecule: {e}")
                self.mcs_mol = None
        else:
            logger.warning("No suitable MCS found")
            return {"error": "No common substructure found"}
        
        # Identify R-group attachment points with error handling
        try:
            self.attachment_points = self._identify_attachment_points()
        except Exception as e:
            logger.warning(f"Error identifying attachment points: {e}")
            self.attachment_points = []
        
        return {
            "mcs_atoms": self.mcs_result.numAtoms,
            "mcs_bonds": self.mcs_result.numBonds,
            "mcs_smarts": self.mcs_result.smartsString,
            "attachment_points": len(self.attachment_points),
            "valid_molecules": len(self.valid_mols)
        }
    
    def _find_robust_mcs(self, mols: List[Chem.Mol]) -> Optional[Any]:
        """Find MCS using multiple fallback approaches."""
        
        approaches = [
            # Approach 1: Simple MCS
            lambda: rdFMCS.FindMCS(mols),
            
            # Approach 2: With threshold
            lambda: rdFMCS.FindMCS(mols, threshold=0.8, timeout=30),
            
            # Approach 3: Ring-focused
            lambda: rdFMCS.FindMCS(mols, ringMatchesRingOnly=True, timeout=30),
            
            # Approach 4: Subset if too many molecules
            lambda: rdFMCS.FindMCS(mols[:50], timeout=60) if len(mols) > 50 else None
        ]
        
        for i, approach in enumerate(approaches, 1):
            try:
                logger.info(f"Trying MCS approach {i}...")
                result = approach()
                if result and result.numAtoms >= 5:  # Minimum meaningful core size
                    logger.info(f"MCS approach {i} succeeded: {result.numAtoms} atoms")
                    return result
                else:
                    logger.info(f"MCS approach {i} gave small result: {result.numAtoms if result else 0} atoms")
            except Exception as e:
                logger.warning(f"MCS approach {i} failed: {e}")
        
        logger.error("All MCS approaches failed")
        return None
    
    def _identify_attachment_points(self) -> List[Dict[str, Any]]:
        """Identify potential R-group attachment points in the core structure."""
        
        if not self.mcs_mol:
            return []
        
        attachment_points = []
        
        try:
            # Strategy 1: Find atoms with free valence or low connectivity
            for i, atom in enumerate(self.mcs_mol.GetAtoms()):
                try:
                    degree = atom.GetDegree()
                    symbol = atom.GetSymbol()
                    
                    # Use degree-based analysis instead of valence (which causes issues)
                    is_attachment = False
                    
                    if symbol == 'C' and degree <= 3:
                        is_attachment = True
                    elif symbol == 'N' and degree <= 2:
                        is_attachment = True
                    elif symbol == 'O' and degree <= 1:
                        is_attachment = True
                    
                    if is_attachment:
                        attachment_points.append({
                            'atom_idx': i,
                            'atom_symbol': symbol,
                            'degree': degree,
                            'r_group_label': f'R{len(attachment_points) + 1}'
                        })
                except Exception as e:
                    logger.warning(f"Error processing atom {i}: {e}")
                    continue
        
        except Exception as e:
            logger.warning(f"Error in attachment point identification: {e}")
            # Fallback: create generic attachment points
            try:
                num_atoms = self.mcs_mol.GetNumAtoms()
                for i in range(min(4, num_atoms)):  # Maximum 4 R-groups
                    attachment_points.append({
                        'atom_idx': i,
                        'atom_symbol': 'C',
                        'degree': 2,
                        'r_group_label': f'R{i + 1}'
                    })
            except:
                pass
        
        logger.info(f"Identified {len(attachment_points)} potential attachment points")
        return attachment_points
    
    def generate_markush_structure_methods(self, output_prefix: str) -> Dict[str, str]:
        """
        Generate Markush structure using multiple methods to avoid '?' symbols.
        
        Args:
            output_prefix: Prefix for output file names
            
        Returns:
            Dictionary of generated image paths
        """
        
        if not self.valid_mols:
            logger.error("No valid molecules available for structure generation")
            return {}
        
        output_dir = Path("./New_version")
        output_dir.mkdir(exist_ok=True)
        
        generated_images = {}
        
        # Method 1: Template-based with highlights (PRIORITY - preserves ring structures)
        template_path = output_dir / f"{output_prefix}_template_markush.png"
        if self._generate_template_markush(template_path):
            generated_images['template'] = str(template_path)
            logger.info("✅ Template method generated - preserves ring structures correctly")
        else:
            logger.warning("❌ Template method failed - this is the preferred method for ring systems")
        
        # Method 2: Core structure with R-group labels (if MCS available and rings preserved)
        if self.mcs_mol:
            core_path = output_dir / f"{output_prefix}_core_with_rgroups.png"
            if self._generate_core_with_rgroups(core_path):
                generated_images['core_with_rgroups'] = str(core_path)
                logger.info("✅ Core with R-groups generated")
            else:
                logger.warning("⚠️  Core with R-groups method failed")
        
        # Method 3: Representative compound (fallback only)
        repr_path = output_dir / f"{output_prefix}_representative.png"
        if self._generate_representative_structure(repr_path):
            generated_images['representative'] = str(repr_path)
            logger.info("✅ Representative structure generated (fallback)")
        
        # Method 4: Simplified scaffold (may lose ring information)
        if self.mcs_mol:
            scaffold_path = output_dir / f"{output_prefix}_scaffold.png"
            if self._generate_simplified_scaffold(scaffold_path):
                generated_images['scaffold'] = str(scaffold_path)
                logger.warning("⚠️  Scaffold method may not preserve ring structures properly")
        
        # Method 5: Additional examples only if template method failed
        if 'template' not in generated_images:
            logger.info("Generating additional examples since template method failed")
            for i in range(min(3, len(self.valid_mols))):
                fallback_path = output_dir / f"{output_prefix}_example_{i+1}.png"
                try:
                    mol = self.valid_mols[i]
                    img = Draw.MolToImage(mol, size=(1000, 800))
                    img.save(str(fallback_path))
                    generated_images[f'example_{i+1}'] = str(fallback_path)
                    logger.info(f"Generated example structure {i+1}")
                except Exception as e:
                    logger.warning(f"Failed to generate example {i+1}: {e}")
        
        # Choose the best one as the main Markush structure
        main_markush_path = output_dir / f"{output_prefix}_markush_structure.png"
        best_method = self._select_best_markush_image(generated_images)
        
        if best_method:
            import shutil
            shutil.copy2(generated_images[best_method], main_markush_path)
            generated_images['main'] = str(main_markush_path)
            logger.info(f"Selected {best_method} method for main Markush structure")
        
        logger.info(f"Generated {len(generated_images)} structure images")
        return generated_images
    
    def _generate_template_markush(self, output_path: Path) -> bool:
        """Generate Markush structure using the first compound as template with MCS highlighted."""
        
        try:
            if not self.valid_mols or not self.mcs_mol:
                return False
            
            template_mol = self.valid_mols[0]
            
            # Find MCS match in template
            match = template_mol.GetSubstructMatch(self.mcs_mol)
            
            if match:
                # Create highlight colors - core in blue, rest in light gray
                highlight_colors = {i: (0.8, 0.9, 1.0) for i in match}  # Light blue for core
                
                # Add R-group labels to non-core atoms
                mol_copy = Chem.Mol(template_mol)
                r_group_counter = 1
                
                for atom_idx in range(mol_copy.GetNumAtoms()):
                    if atom_idx not in match:
                        atom = mol_copy.GetAtomWithIdx(atom_idx)
                        # Check if this atom is connected to the core
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetIdx() in match:
                                atom.SetProp('atomLabel', f'R{r_group_counter}')
                                r_group_counter += 1
                                break
                
                img = Draw.MolToImage(
                    mol_copy,
                    size=(1000, 800),
                    highlightAtoms=match,
                    highlightColors=highlight_colors,
                    kekulize=True
                )
                
                img.save(str(output_path))
                logger.info(f"Template Markush structure saved: {output_path}")
                return True
            else:
                # Fallback: just render the template
                img = Draw.MolToImage(template_mol, size=(1000, 800))
                img.save(str(output_path))
                logger.info(f"Template structure (no highlight) saved: {output_path}")
                return True
                
        except Exception as e:
            logger.error(f"Template Markush generation failed: {e}")
            return False
    
    def _generate_core_with_rgroups(self, output_path: Path) -> bool:
        """Generate core structure with proper R-group labels."""
        
        try:
            if not self.mcs_mol:
                return False
            
            # Try to convert SMARTS to SMILES for better visualization
            try:
                core_smiles = Chem.MolToSmiles(self.mcs_mol)
                if core_smiles:
                    display_mol = Chem.MolFromSmiles(core_smiles)
                    if display_mol:
                        mol_to_render = display_mol
                        logger.info("Using SMILES representation for core structure")
                    else:
                        mol_to_render = self.mcs_mol
                else:
                    mol_to_render = self.mcs_mol
            except:
                mol_to_render = self.mcs_mol
            
            # Add R-group labels based on attachment points
            mol_copy = Chem.Mol(mol_to_render)
            
            for ap in self.attachment_points:
                atom_idx = ap['atom_idx']
                if atom_idx < mol_copy.GetNumAtoms():
                    atom = mol_copy.GetAtomWithIdx(atom_idx)
                    atom.SetProp('atomLabel', ap['r_group_label'])
            
            # If no specific attachment points, label terminal/low-degree atoms
            if not self.attachment_points:
                r_counter = 1
                for atom in mol_copy.GetAtoms():
                    if atom.GetDegree() <= 2 and atom.GetSymbol() in ['C', 'N', 'O']:
                        atom.SetProp('atomLabel', f'R{r_counter}')
                        r_counter += 1
                        if r_counter > 6:  # Limit to R1-R6
                            break
            
            img = Draw.MolToImage(
                mol_copy,
                size=(1000, 800),
                kekulize=True
            )
            
            img.save(str(output_path))
            logger.info(f"Core with R-groups saved: {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"Core with R-groups generation failed: {e}")
            return False
    
    def _generate_representative_structure(self, output_path: Path) -> bool:
        """Generate a representative compound structure."""
        
        try:
            if not self.valid_mols:
                return False
            
            # Choose a representative molecule (middle of the list)
            repr_mol = self.valid_mols[len(self.valid_mols) // 2]
            
            img = Draw.MolToImage(repr_mol, size=(1000, 800))
            img.save(str(output_path))
            
            logger.info(f"Representative structure saved: {output_path}")
            return True
            
        except Exception as e:
            logger.error(f"Representative structure generation failed: {e}")
            return False
    
    def _generate_simplified_scaffold(self, output_path: Path) -> bool:
        """Generate a simplified scaffold structure."""
        
        try:
            if not self.mcs_mol:
                return False
            
            # Create a simplified version by removing complex SMARTS features
            mcs_smarts = Chem.MolToSmarts(self.mcs_mol)
            
            # Simplify the SMARTS
            simplified_smarts = mcs_smarts
            simplifications = [
                ('&!R', ''),  # Remove ring/non-ring specifications
                ('&!@', ''),  # Remove ring bond specifications
                (':&@', ':'), # Simplify aromatic bonds
                ('!R', ''),   # Remove more ring specifications
                ('!@', ''),   # Remove more bond specifications
            ]
            
            for old, new in simplifications:
                simplified_smarts = simplified_smarts.replace(old, new)
            
            # Try to create molecule from simplified SMARTS
            simplified_mol = Chem.MolFromSmarts(simplified_smarts)
            
            if simplified_mol:
                # Add generic R-group labels
                for i, atom in enumerate(simplified_mol.GetAtoms()):
                    if atom.GetDegree() <= 2:
                        atom.SetProp('atomLabel', f'R{(i % 6) + 1}')
                
                img = Draw.MolToImage(simplified_mol, size=(1000, 800))
                img.save(str(output_path))
                
                logger.info(f"Simplified scaffold saved: {output_path}")
                return True
            else:
                logger.warning("Could not create simplified scaffold")
                return False
                
        except Exception as e:
            logger.error(f"Simplified scaffold generation failed: {e}")
            return False
    
    def _select_best_markush_image(self, generated_images: Dict[str, str]) -> Optional[str]:
        """Select the best Markush structure image from generated options."""
        
        # Priority order: template first (shows rings correctly) > core_with_rgroups > scaffold > representative > examples
        # Template method is prioritized because it correctly preserves ring structures
        priority_order = [
            'template', 'core_with_rgroups', 'scaffold', 'representative',
            'example_1', 'example_2', 'example_3'
        ]
        
        for method in priority_order:
            if method in generated_images and Path(generated_images[method]).exists():
                logger.info(f"Selected {method} method (shows proper ring structures)")
                return method
        
        # If nothing from priority list, return any available method
        for method, path in generated_images.items():
            if Path(path).exists():
                logger.warning(f"Fallback to {method} method")
                return method
        
        return None

# Utility functions for integration

def generate_markush_structure_general(smiles_list: List[str], output_prefix: str) -> Dict[str, Any]:
    """
    General function to generate Markush structures from any SMILES list.
    
    Args:
        smiles_list: List of SMILES strings
        output_prefix: Prefix for output files
        
    Returns:
        Dictionary with analysis results and generated image paths
    """
    
    generator = GeneralMarkushGenerator()
    
    # Analyze the compound set
    analysis = generator.analyze_compound_set(smiles_list)
    
    if "error" in analysis:
        logger.error(f"Analysis failed: {analysis['error']}")
        return analysis
    
    # Generate structure images
    images = generator.generate_markush_structure_methods(output_prefix)
    
    return {
        "analysis": analysis,
        "images": images,
        "success": len(images) > 0
    }

def process_csv_to_markush(csv_path: str, smiles_col: str = 'Smiles', output_prefix: Optional[str] = None) -> Dict[str, Any]:
    """
    Process CSV file to generate Markush structures.
    
    Args:
        csv_path: Path to CSV file
        smiles_col: Name of SMILES column
        output_prefix: Output prefix (derived from filename if None)
        
    Returns:
        Results dictionary
    """
    
    try:
        # Load CSV
        df = pd.read_csv(csv_path)
        logger.info(f"Loaded CSV with {len(df)} rows")
        
        # Find SMILES column with flexible matching
        if smiles_col not in df.columns:
            # Try exact matches first
            smiles_variants = ['SMILES', 'smiles', 'Smiles', 'smi', 'SMI']
            for variant in smiles_variants:
                if variant in df.columns:
                    smiles_col = variant
                    break
            else:
                # Try case-insensitive and stripped matches
                for col in df.columns:
                    col_clean = col.strip().lower()
                    if col_clean in ['smiles', 'smi']:
                        smiles_col = col  # Use original column name
                        break
                else:
                    # Look for columns containing 'smiles'
                    for col in df.columns:
                        col_clean = col.strip().lower()
                        if 'smiles' in col_clean or 'smi' in col_clean:
                            smiles_col = col  # Use original column name
                            break
                    else:
                        raise ValueError(f"SMILES column not found. Available columns: {list(df.columns)}")
        
        smiles_list = df[smiles_col].dropna().tolist()
        logger.info(f"Found {len(smiles_list)} SMILES strings")
        
        # Generate output prefix
        if output_prefix is None:
            output_prefix = Path(csv_path).stem
        
        # Generate Markush structures
        results = generate_markush_structure_general(smiles_list, output_prefix)
        
        # Add CSV information
        results['csv_info'] = {
            'file': csv_path,
            'total_rows': len(df),
            'smiles_column': smiles_col,
            'valid_smiles': len(smiles_list)
        }
        
        return results
        
    except Exception as e:
        logger.error(f"CSV processing failed: {e}")
        return {"error": str(e)}

# Integration with existing patent generator

def get_enhanced_markush_images(csv_path: str, output_prefix: str) -> Dict[str, str]:
    """
    Get enhanced Markush structure images for patent generation.
    
    This function replaces the problematic SMARTS-based approach with
    the new general method that avoids '?' symbols.
    
    Args:
        csv_path: Path to CSV file
        output_prefix: Output prefix
        
    Returns:
        Dictionary of image paths
    """
    
    results = process_csv_to_markush(csv_path, output_prefix=output_prefix)
    
    if "error" in results:
        logger.error(f"Markush generation failed: {results['error']}")
        return {}
    
    images = results.get('images', {})
    
    # Return standardized image paths for patent generator
    return {
        'markush_structure': images.get('main', ''),
        'core_structure': images.get('core_with_rgroups', ''),
        'structure_image': images.get('template', ''),
        'representative': images.get('representative', ''),
        'scaffold': images.get('scaffold', '')
    }