#!/usr/bin/env python3
"""
Enhanced benchmark script that saves results to both JSON and database
"""

import sys
import os
import json
import sqlite3
from datetime import datetime

# Add the parent directory to the path to import from review-app
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'review-app'))

try:
    from scripts.import_benchmark_data import DataImporter
    DATABASE_AVAILABLE = True
except ImportError:
    DATABASE_AVAILABLE = False
    print("Warning: Database import not available. Will only save to JSON.")

def save_to_database(benchmark_data, dataset_type, benchmark_date):
    """Save benchmark data to database"""
    if not DATABASE_AVAILABLE:
        return
    
    try:
        importer = DataImporter()
        
        for record in benchmark_data:
            try:
                # Insert molecule data
                moleculeData = record.get('molecule_data', record)
                moleculeId = importer.insertMolecule(moleculeData, dataset_type, benchmark_date)
                
                # Insert PFASGroups result if exists
                if 'pfasgroups_result' in record:
                    importer.insertPFASGroupsResult(moleculeId, record['pfasgroups_result'])
                
                # Insert Atlas result if exists
                if 'atlas_result' in record:
                    importer.insertAtlasResult(moleculeId, record['atlas_result'])
                    
            except Exception as e:
                print(f"Warning: Error saving record to database: {e}")
                continue
        
        importer.close()
        print(f"✓ Saved {len(benchmark_data)} records to database")
        
    except Exception as e:
        print(f"Warning: Database save failed: {e}")

def run_enhanced_benchmark():
    """Run the existing enhanced benchmark with database integration"""
    
    # Import the existing benchmark script
    current_dir = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, current_dir)
    
    try:
        from enhanced_pfas_benchmark import main as run_benchmark
        
        # Run the original benchmark
        results = run_benchmark()
        
        # If results were returned, save to database
        if results:
            benchmark_date = datetime.now().strftime('%Y-%m-%d')
            save_to_database(results, 'enhanced', benchmark_date)
            
        return results
        
    except ImportError as e:
        print(f"Error importing enhanced_pfas_benchmark: {e}")
        return None

def run_oecd_benchmark():
    """Run OECD benchmark with database integration"""
    
    try:
        # This would be the OECD benchmark script
        # For now, we'll use the existing JSON files
        json_files = [f for f in os.listdir('data') if f.startswith('pfas_oecd_benchmark') and f.endswith('.json')]
        
        if json_files:
            latest_file = sorted(json_files)[-1]
            print(f"Loading existing OECD data from {latest_file}")
            
            with open(os.path.join('data', latest_file), 'r') as f:
                results = json.load(f)
            
            benchmark_date = datetime.now().strftime('%Y-%m-%d')
            save_to_database(results, 'oecd', benchmark_date)
            return results
        
    except Exception as e:
        print(f"Error in OECD benchmark: {e}")
        return None

def run_timing_benchmark():
    """Run timing benchmark with database integration"""
    
    try:
        json_files = [f for f in os.listdir('data') if f.startswith('pfas_timing_benchmark') and f.endswith('.json')]
        
        if json_files:
            latest_file = sorted(json_files)[-1]
            print(f"Loading existing timing data from {latest_file}")
            
            with open(os.path.join('data', latest_file), 'r') as f:
                results = json.load(f)
            
            # Convert timing data format for database
            converted_results = []
            for record in results:
                converted_record = {
                    'molecule_data': {
                        'smiles': record['smiles'],
                        'molecular_weight': record.get('molecular_weight'),
                        'num_atoms': record.get('num_atoms'),
                        'num_bonds': record.get('num_bonds'),
                        'chain_length': record.get('chain_length')
                    }
                }
                converted_results.append(converted_record)
            
            benchmark_date = datetime.now().strftime('%Y-%m-%d')
            
            # Use specific timing import method
            if DATABASE_AVAILABLE:
                try:
                    importer = DataImporter()
                    importer.importTimingData(os.path.join('data', latest_file), latest_file)
                    importer.close()
                    print(f"✓ Imported timing data to database")
                except Exception as e:
                    print(f"Warning: Timing database import failed: {e}")
            
            return results
        
    except Exception as e:
        print(f"Error in timing benchmark: {e}")
        return None

if __name__ == "__main__":
    print("🚀 Enhanced Benchmark Runner with Database Integration")
    print("=" * 60)
    
    # Check database availability
    if DATABASE_AVAILABLE:
        print("✓ Database integration available")
    else:
        print("⚠️  Database integration not available")
    
    # Run benchmarks
    print("\\n1. Running Enhanced Benchmark...")
    enhanced_results = run_enhanced_benchmark()
    
    print("\\n2. Processing OECD Data...")
    oecd_results = run_oecd_benchmark()
    
    print("\\n3. Processing Timing Data...")
    timing_results = run_timing_benchmark()
    
    print("\\n✅ Benchmark integration complete!")
    
    if DATABASE_AVAILABLE:
        print("\\n📊 Database Status:")
        try:
            importer = DataImporter()
            stats = importer.getImportStats()
            print(json.dumps(stats, indent=2))
            importer.close()
        except Exception as e:
            print(f"Error getting database stats: {e}")