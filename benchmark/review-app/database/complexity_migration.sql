-- ============================================================================
-- Graph Complexity Metrics Migration
-- ============================================================================
-- This migration adds graph complexity metric columns to the molecules table
-- for enhanced timing analysis and performance prediction
-- 
-- Run with: sqlite3 pfas_benchmark.db < complexity_migration.sql
-- ============================================================================

-- Add complexity metric columns
ALTER TABLE molecules ADD COLUMN complexity_diameter REAL;
ALTER TABLE molecules ADD COLUMN complexity_radius REAL;
ALTER TABLE molecules ADD COLUMN complexity_avg_eccentricity REAL;
ALTER TABLE molecules ADD COLUMN complexity_max_eccentricity REAL;
ALTER TABLE molecules ADD COLUMN complexity_avg_degree REAL;
ALTER TABLE molecules ADD COLUMN complexity_density REAL;
ALTER TABLE molecules ADD COLUMN complexity_num_cycles INTEGER;
ALTER TABLE molecules ADD COLUMN complexity_avg_betweenness REAL;
ALTER TABLE molecules ADD COLUMN complexity_max_betweenness REAL;
ALTER TABLE molecules ADD COLUMN complexity_avg_clustering REAL;
ALTER TABLE molecules ADD COLUMN complexity_score REAL;

-- Create indexes for commonly queried complexity metrics
CREATE INDEX IF NOT EXISTS idx_molecules_complexity_score ON molecules(complexity_score);
CREATE INDEX IF NOT EXISTS idx_molecules_complexity_diameter ON molecules(complexity_diameter);
CREATE INDEX IF NOT EXISTS idx_molecules_num_cycles ON molecules(complexity_num_cycles);

-- Verify the migration
SELECT COUNT(*) as total_molecules FROM molecules;

-- Show sample of columns (will be NULL until data is re-imported)
SELECT 
    id, 
    smiles, 
    num_atoms, 
    complexity_score, 
    complexity_diameter 
FROM molecules 
LIMIT 5;

-- Success message
SELECT 'Migration completed successfully! Complexity columns added to molecules table.' as status;
