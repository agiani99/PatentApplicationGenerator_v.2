"""
Universal Patent Generator Runner

This script uses the enhanced general Markush approach to generate
patent applications from any CSV file with SMILES data.

Usage:
    python universal_patent_runner.py <csv_file>
    python universal_patent_runner.py <csv_file> <output_prefix>

Features:
- Works with any CSV file containing SMILES
- Generates proper Markush structures with R-group labels  
- Avoids the '?' symbol problem
- Creates multiple visualization methods
- Saves everything in ./New_version folder
"""

import sys
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def main():
    """Main runner function."""
    
    print("🧬 Universal Patent Generator")
    print("=" * 40)
    
    # Parse command line arguments
    if len(sys.argv) < 2:
        print("❌ Usage: python universal_patent_runner.py <csv_file> [output_prefix]")
        print("\nExamples:")
        print("  python universal_patent_runner.py CDK5_ACT_IC50.csv")
        print("  python universal_patent_runner.py my_compounds.csv my_patent")
        return
    
    csv_file = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else None
    
    # Validate input file
    if not Path(csv_file).exists():
        print(f"❌ CSV file not found: {csv_file}")
        return
    
    print(f"📁 Input file: {csv_file}")
    print(f"📂 Output folder: ./New_version")
    
    if output_prefix:
        print(f"🏷️  Output prefix: {output_prefix}")
    else:
        output_prefix = Path(csv_file).stem
        print(f"🏷️  Output prefix: {output_prefix} (from filename)")
    
    # Ensure New_version directory exists
    Path("./New_version").mkdir(exist_ok=True)
    
    try:
        # Import and run the enhanced generator
        from enhanced_patent_generator import generate_enhanced_patent_universal
        
        print(f"\n🚀 Starting patent generation...")
        results = generate_enhanced_patent_universal(csv_file, output_prefix)
        
        if results.get('success', False):
            print(f"\n✅ SUCCESS!")
            print(f"📄 Patent application: {results['patent_application']}")
            print(f"📋 Metadata: {results['metadata_file']}")
            print(f"🖼️  Generated images:")
            
            for method, path in results['markush_images'].items():
                if Path(path).exists():
                    print(f"   ✅ {method}: {path}")
                else:
                    print(f"   ❌ {method}: {path} (missing)")
            
            print(f"\n💡 Next steps:")
            print(f"1. Review the patent document: {results['patent_application']}")
            print(f"2. Check the Markush structure images")
            print(f"3. Verify R-group labeling is correct")
            print(f"4. Make any final adjustments as needed")
            
        else:
            print(f"\n❌ Patent generation failed")
            if 'error' in results:
                print(f"Error: {results['error']}")
        
    except ImportError as e:
        print(f"❌ Import error: {e}")
        print("💡 Make sure all required files are in the New_version folder")
    except Exception as e:
        print(f"❌ Error during patent generation: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()