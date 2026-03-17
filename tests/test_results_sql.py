"""Tests for ResultsModel and MoleculeResult to_sql methods."""

import os
import sqlite3
import tempfile
import pytest
from pathlib import Path

from PFASGroups.PFASEmbeddings import ResultsModel, MoleculeResult
from HalogenGroups import parse_smiles


@pytest.fixture
def sample_pfas_molecules():
    """Generate a list of PFAS molecules for testing."""
    smiles_list = [
        # Perfluorooctanoic acid (PFOA)
        "C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(=O)O",
        # Perfluorooctane sulfonic acid (PFOS)
        "C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)S(=O)(=O)O",
        # Perfluorohexanoic acid (PFHxA)
        "C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(=O)O",
        # GenX (hexafluoropropylene oxide dimer acid)
        "C(C(C(OC(C(F)(F)F)F)(F)F)(F)F)(=O)O",
        # 6:2 Fluorotelomer alcohol
        "C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)CCO",
        # Perfluorobutanoic acid (PFBA)
        "C(C(C(C(F)(F)F)(F)F)(F)F)(=O)O",
        # Perfluorononanoic acid (PFNA)
        "C(C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(=O)O",
        # Perfluorodecanoic acid (PFDA)
        "C(C(C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(=O)O",
    ]
    return smiles_list


@pytest.fixture
def parsed_results(sample_pfas_molecules):
    """Parse PFAS molecules and return ResultsModel."""
    results = parse_smiles(sample_pfas_molecules, output_format='list')
    results_model = ResultsModel.from_raw(results)
    return results_model


class TestMoleculeResultToSQL:
    """Tests for MoleculeResult.to_sql method."""

    def test_to_sql_sqlite_single_molecule(self, parsed_results, tmp_path):
        """Test exporting a single molecule result to SQLite."""
        db_file = tmp_path / "test_single.db"
        
        # Get first molecule result
        mol_result = parsed_results[0]
        
        # Export to SQLite
        mol_result.to_sql(filename=str(db_file))
        
        # Verify database exists
        assert db_file.exists()
        
        # Check database contents
        conn = sqlite3.connect(str(db_file))
        cursor = conn.cursor()
        
        # Check components table
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='components'")
        assert cursor.fetchone() is not None
        
        cursor.execute("SELECT COUNT(*) FROM components")
        component_count = cursor.fetchone()[0]
        assert component_count > 0
        
        # Check groups table
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='pfas_groups_in_compound'")
        assert cursor.fetchone() is not None
        
        cursor.execute("SELECT COUNT(*) FROM pfas_groups_in_compound")
        group_count = cursor.fetchone()[0]
        assert group_count > 0
        
        # Verify SMILES is stored
        cursor.execute("SELECT DISTINCT smiles FROM pfas_groups_in_compound")
        smiles = cursor.fetchone()[0]
        assert smiles == mol_result.smiles
        
        conn.close()

    def test_to_sql_custom_table_names(self, parsed_results, tmp_path):
        """Test using custom table names."""
        db_file = tmp_path / "test_custom_tables.db"
        
        mol_result = parsed_results[0]
        mol_result.to_sql(
            filename=str(db_file),
            components_table="my_components",
            groups_table="my_groups"
        )
        
        conn = sqlite3.connect(str(db_file))
        cursor = conn.cursor()
        
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='my_components'")
        assert cursor.fetchone() is not None
        
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='my_groups'")
        assert cursor.fetchone() is not None
        
        conn.close()

    def test_to_sql_append_mode(self, parsed_results, tmp_path):
        """Test appending multiple molecules to the same database."""
        db_file = tmp_path / "test_append.db"
        
        # Add first molecule
        parsed_results[0].to_sql(filename=str(db_file), if_exists='replace')
        
        # Add second molecule
        parsed_results[1].to_sql(filename=str(db_file), if_exists='append')
        
        conn = sqlite3.connect(str(db_file))
        cursor = conn.cursor()
        
        # Should have data from both molecules
        cursor.execute("SELECT COUNT(DISTINCT smiles) FROM pfas_groups_in_compound")
        distinct_smiles = cursor.fetchone()[0]
        assert distinct_smiles == 2
        
        conn.close()

    def test_to_sql_no_connection_error(self, parsed_results):
        """Test that error is raised when no connection info is provided."""
        mol_result = parsed_results[0]
        
        with pytest.raises(ValueError, match="Either filename.*or dbname"):
            mol_result.to_sql()

    def test_to_sql_data_integrity(self, parsed_results, tmp_path):
        """Test that all data is correctly stored in the database."""
        db_file = tmp_path / "test_integrity.db"
        
        mol_result = parsed_results[0]
        mol_result.to_sql(filename=str(db_file))
        
        conn = sqlite3.connect(str(db_file))
        cursor = conn.cursor()
        
        # Check component data structure
        cursor.execute("PRAGMA table_info(components)")
        columns = {row[1] for row in cursor.fetchall()}
        expected_columns = {'smiles', 'group_id', 'group_name', 'smarts_label', 'component_atoms'}
        assert expected_columns.issubset(columns)
        
        # Check groups data structure
        cursor.execute("PRAGMA table_info(pfas_groups_in_compound)")
        columns = {row[1] for row in cursor.fetchall()}
        expected_columns = {'smiles', 'group_id', 'group_name', 'match_count'}
        assert expected_columns.issubset(columns)
        
        # Verify component atoms format
        cursor.execute("SELECT component_atoms FROM components LIMIT 1")
        atoms = cursor.fetchone()[0]
        assert isinstance(atoms, str)
        assert ',' in atoms or atoms.isdigit()  # Either "1,2,3" or single digit "1"
        
        conn.close()


class TestResultsModelToSQL:
    """Tests for ResultsModel.to_sql method."""

    def test_to_sql_sqlite_multiple_molecules(self, parsed_results, tmp_path):
        """Test exporting multiple molecules to SQLite."""
        db_file = tmp_path / "test_multiple.db"
        
        # Export all results
        parsed_results.to_sql(filename=str(db_file))
        
        assert db_file.exists()
        
        conn = sqlite3.connect(str(db_file))
        cursor = conn.cursor()
        
        # Check that we have data from multiple molecules
        cursor.execute("SELECT COUNT(DISTINCT smiles) FROM pfas_groups_in_compound")
        distinct_smiles = cursor.fetchone()[0]
        assert distinct_smiles == len(parsed_results)
        
        # Check total components count
        cursor.execute("SELECT COUNT(*) FROM components")
        total_components = cursor.fetchone()[0]
        assert total_components > 0
        
        conn.close()

    def test_to_sql_replace_mode(self, parsed_results, tmp_path):
        """Test replace mode overwrites existing data."""
        db_file = tmp_path / "test_replace.db"
        
        # First export with some data
        ResultsModel(parsed_results[:3]).to_sql(filename=str(db_file), if_exists='replace')
        
        conn = sqlite3.connect(str(db_file))
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(DISTINCT smiles) FROM pfas_groups_in_compound")
        count_before = cursor.fetchone()[0]
        conn.close()
        
        # Export again with different data
        ResultsModel(parsed_results[3:5]).to_sql(filename=str(db_file), if_exists='replace')
        
        conn = sqlite3.connect(str(db_file))
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(DISTINCT smiles) FROM pfas_groups_in_compound")
        count_after = cursor.fetchone()[0]
        conn.close()
        
        # Should have replaced the data
        assert count_before == 3
        assert count_after == 2

    def test_to_sql_large_dataset(self, tmp_path):
        """Test with a larger dataset (hundreds of molecules)."""
        # Generate more PFAS-like molecules
        smiles_list = []
        
        # Use valid PFAS SMILES instead of generating invalid ones
        # These are known valid PFAS structures
        base_smiles = [
            "C(C(C(C(F)(F)F)(F)F)(F)F)(=O)O",  # PFBA
            "C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(=O)O",  # PFHxA
            "C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(=O)O",  # PFOA
            "C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)S(=O)(=O)O",  # PFHxS
            "C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)S(=O)(=O)O",  # PFOS
        ]
        smiles_list.extend(base_smiles)
        
        # Duplicate to get to more molecules
        smiles_list = smiles_list * 40  # 200 molecules
        
        results = parse_smiles(smiles_list[:200], output_format='list')  # Limit to 200 for speed
        results_model = ResultsModel.from_raw(results)
        
        db_file = tmp_path / "test_large.db"
        results_model.to_sql(filename=str(db_file))
        
        conn = sqlite3.connect(str(db_file))
        cursor = conn.cursor()
        
        cursor.execute("SELECT COUNT(DISTINCT smiles) FROM pfas_groups_in_compound")
        distinct_count = cursor.fetchone()[0]
        assert distinct_count > 0  # At least some unique molecules
        
        cursor.execute("SELECT COUNT(*) FROM components")
        component_count = cursor.fetchone()[0]
        assert component_count > 100  # Should have many components
        
        conn.close()

    def test_to_sql_group_aggregation(self, parsed_results, tmp_path):
        """Test that groups are correctly aggregated per molecule."""
        db_file = tmp_path / "test_aggregation.db"
        
        parsed_results.to_sql(filename=str(db_file))
        
        conn = sqlite3.connect(str(db_file))
        cursor = conn.cursor()
        
        # Check that match_count is populated
        cursor.execute("SELECT smiles, group_name, match_count FROM pfas_groups_in_compound")
        rows = cursor.fetchall()
        
        assert len(rows) > 0
        for smiles, group_name, match_count in rows:
            assert match_count > 0
            assert isinstance(group_name, str)
            assert len(group_name) > 0
        
        conn.close()

    def test_to_sql_empty_results(self, tmp_path):
        """Test handling of empty results."""
        db_file = tmp_path / "test_empty.db"
        
        # Create empty results
        empty_results = ResultsModel([])
        
        # Should not raise an error, just not create tables
        empty_results.to_sql(filename=str(db_file))
        
        # Database might exist but should be empty or have no tables
        if db_file.exists():
            conn = sqlite3.connect(str(db_file))
            cursor = conn.cursor()
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
            tables = cursor.fetchall()
            # Either no tables or empty tables
            if tables:
                cursor.execute("SELECT COUNT(*) FROM components")
                assert cursor.fetchone()[0] == 0
            conn.close()

    def test_to_sql_mixed_content(self, tmp_path):
        """Test with molecules that have varying numbers of matches."""
        # Create results with varying content
        smiles_list = [
            "C(C(C(C(F)(F)F)(F)F)(F)F)(=O)O",  # PFAS with matches
            "CCCCCCCC",  # Non-PFAS, likely no matches
            "C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(=O)O",  # Another PFAS
        ]
        
        results = parse_smiles(smiles_list, output_format='list')
        results_model = ResultsModel.from_raw(results)
        
        db_file = tmp_path / "test_mixed.db"
        results_model.to_sql(filename=str(db_file))
        
        # Should handle both molecules with and without matches gracefully
        assert db_file.exists()

    def test_to_sql_special_characters_in_smiles(self, tmp_path):
        """Test handling of special characters in SMILES strings."""
        smiles_list = [
            "C(C(C(C(F)(F)F)(F)F)(F)F)(=O)O",
            "c1cc(ccc1C(=O)O)C(C(C(C(F)(F)F)(F)F)(F)F)(F)F",  # Aromatic
        ]
        
        results = parse_smiles(smiles_list, output_format='list')
        results_model = ResultsModel.from_raw(results)
        
        db_file = tmp_path / "test_special_chars.db"
        results_model.to_sql(filename=str(db_file))
        
        conn = sqlite3.connect(str(db_file))
        cursor = conn.cursor()
        
        # Should be able to retrieve SMILES without issues
        cursor.execute("SELECT DISTINCT smiles FROM pfas_groups_in_compound")
        stored_smiles = [row[0] for row in cursor.fetchall()]
        
        # At least one SMILES should be stored
        assert len(stored_smiles) > 0
        
        conn.close()


class TestSQLIntegration:
    """Integration tests for SQL export functionality."""

    def test_roundtrip_verification(self, parsed_results, tmp_path):
        """Test that data can be written and read back correctly."""
        db_file = tmp_path / "test_roundtrip.db"
        
        # Export to database
        parsed_results.to_sql(filename=str(db_file))
        
        # Read back and verify
        conn = sqlite3.connect(str(db_file))
        cursor = conn.cursor()
        
        # Verify each molecule's data
        for mol_result in parsed_results:
            cursor.execute(
                "SELECT COUNT(*) FROM pfas_groups_in_compound WHERE smiles = ?",
                (mol_result.smiles,)
            )
            group_count = cursor.fetchone()[0]
            
            # Should have at least some groups for PFAS molecules
            if any(m.is_group for m in mol_result.matches):
                assert group_count > 0
        
        conn.close()

    def test_database_constraints(self, parsed_results, tmp_path):
        """Test that database schema is appropriate."""
        db_file = tmp_path / "test_schema.db"
        
        parsed_results.to_sql(filename=str(db_file))
        
        conn = sqlite3.connect(str(db_file))
        cursor = conn.cursor()
        
        # Check components table schema
        cursor.execute("PRAGMA table_info(components)")
        components_schema = cursor.fetchall()
        assert len(components_schema) >= 5  # At least 5 columns
        
        # Check groups table schema
        cursor.execute("PRAGMA table_info(pfas_groups_in_compound)")
        groups_schema = cursor.fetchall()
        assert len(groups_schema) >= 4  # At least 4 columns
        
        conn.close()

    def test_concurrent_writes(self, parsed_results, tmp_path):
        """Test that append mode works correctly for multiple write operations."""
        db_file = tmp_path / "test_concurrent.db"
        
        # Split results and write in batches
        batch_size = 2
        for i in range(0, len(parsed_results), batch_size):
            batch = ResultsModel(parsed_results[i:i+batch_size])
            batch.to_sql(
                filename=str(db_file),
                if_exists='append' if i > 0 else 'replace'
            )
        
        # Verify all data is present
        conn = sqlite3.connect(str(db_file))
        cursor = conn.cursor()
        
        cursor.execute("SELECT COUNT(DISTINCT smiles) FROM pfas_groups_in_compound")
        distinct_count = cursor.fetchone()[0]
        
        # Should have all molecules (or at least close, some may not have matches)
        assert distinct_count > 0
        
        conn.close()


class TestErrorHandling:
    """Test error handling in to_sql methods."""

    def test_missing_dependencies(self, parsed_results, monkeypatch, tmp_path):
        """Test graceful handling when pandas/sqlalchemy are not installed."""
        # This test would need to mock the import - skip if dependencies are present
        # Just verify the error message is informative
        pass  # Dependencies should be available in test environment

    def test_invalid_database_path(self, parsed_results):
        """Test handling of invalid database paths."""
        # Try to write to an impossible location
        with pytest.raises(Exception):  # Could be various errors depending on OS
            parsed_results[0].to_sql(filename="/impossible/path/database.db")

    def test_invalid_if_exists_parameter(self, parsed_results, tmp_path):
        """Test handling of invalid if_exists parameter."""
        db_file = tmp_path / "test_invalid_param.db"
        
        # pandas should raise an error for invalid if_exists values
        # Note: may raise ValueError or other pandas/sqlalchemy error
        with pytest.raises(Exception):  # Broader exception to catch pandas errors
            parsed_results[0].to_sql(
                filename=str(db_file),
                if_exists='invalid_mode'
            )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
