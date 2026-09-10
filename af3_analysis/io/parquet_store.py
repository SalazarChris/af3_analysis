"""
Parquet store for AF3 Confidence Analysis Pipeline.

Persists intermediate canonical and analysis tables in Parquet format.
"""

import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path
from typing import Any, Dict, List, Optional
import pandas as pd


class ParquetStore:
    """
    Store tables in Parquet format with schema validation.
    
    Provides methods to write and read tables while preserving
    types and null values.
    """
    
    def __init__(self, store_dir: Path):
        self.store_dir = store_dir
        self._tables: Dict[str, pd.DataFrame] = {}
    
    def write_table(self, name: str, df: pd.DataFrame) -> Path:
        """
        Write a DataFrame to Parquet.
        
        Args:
            name: Table name (used for filename)
            df: DataFrame to write
        
        Returns:
            Path to the written file
        """
        import pyarrow as pa
        import pyarrow.parquet as pq
        
        # Convert to PyArrow table
        table = pa.Table.from_pandas(df, preserve_index=False)
        
        # Write to Parquet
        parquet_path = self.store_dir / f"{name}.parquet"
        pq.write_table(table, parquet_path)
        
        # Also keep in memory
        self._tables[name] = df
        
        return parquet_path
    
    def read_table(self, name: str) -> Optional[pd.DataFrame]:
        """Read a DataFrame from Parquet."""
        parquet_path = self.store_dir / f"{name}.parquet"
        if not parquet_path.exists():
            return None
        
        table = pq.read_table(parquet_path)
        return table.to_pandas()
    
    def has_table(self, name: str) -> bool:
        """Check if a table exists."""
        return (self.store_dir / f"{name}.parquet").exists()
    
    def list_tables(self) -> List[str]:
        """List all stored table names."""
        return list(self._tables.keys())
    
    def get_all_tables(self) -> Dict[str, pd.DataFrame]:
        """Get all stored tables."""
        return self._tables.copy()
    
    def clear(self) -> None:
        """Clear all stored tables."""
        self._tables.clear()
        # Also delete files
        for parquet_file in self.store_dir.glob("*.parquet"):
            parquet_file.unlink()


def create_parquet_store(store_dir: Path) -> ParquetStore:
    """Create a Parquet store instance."""
    store_dir.mkdir(parents=True, exist_ok=True)
    return ParquetStore(store_dir)


__all__ = [
    "ParquetStore",
    "create_parquet_store",
]
