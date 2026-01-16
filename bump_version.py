#!/usr/bin/env python3
"""
Bump version in Cargo.toml and Python version.py files.

Usage:
    python bump_version.py 0.1.5
"""

import sys
import re
import pathlib
from typing import Optional


def bump_version_in_file(filepath: pathlib.Path, old_version: str, new_version: str) -> bool:
    """
    Update version string in a file.
    
    Args:
        filepath: Path to file to update
        old_version: Old version string
        new_version: New version string
        
    Returns:
        True if file was updated, False otherwise
    """
    with open(filepath, "r") as f:
        content = f.read()
    
    if old_version not in content:
        return False
    
    new_content = content.replace(old_version, new_version)
    with open(filepath, "w") as f:
        f.write(new_content)
    
    return True


def get_current_version() -> str:
    """Get current version from Cargo.toml."""
    cargo_path = pathlib.Path(__file__).parent / "Cargo.toml"
    with open(cargo_path) as f:
        for line in f:
            if line.startswith("version"):
                match = re.search(r'version\s*=\s*"([^"]+)"', line)
                if match:
                    return match.group(1)
    raise ValueError("Could not find version in Cargo.toml")


def main(new_version: str) -> int:
    """
    Bump version in all version files.
    
    Args:
        new_version: New version to set
        
    Returns:
        0 on success, 1 on failure
    """
    project_root = pathlib.Path(__file__).parent
    
    # Get current version
    old_version = get_current_version()
    print(f"Current version: {old_version}")
    print(f"New version: {new_version}")
    
    # Files to update
    files_to_update = [
        ("Cargo.toml", project_root / "Cargo.toml"),
        ("pyproject.toml", project_root / "pyproject.toml"),
        ("python/mccnado/version.py", project_root / "python" / "mccnado" / "version.py"),
    ]
    
    # Update all files
    updated_files = []
    for name, filepath in files_to_update:
        if filepath.exists():
            if bump_version_in_file(filepath, old_version, new_version):
                print(f"✓ Updated {name}")
                updated_files.append(name)
            else:
                print(f"✗ Version not found in {name}")
        else:
            print(f"✗ File not found: {name}")
    
    if not updated_files:
        print("\nError: No files were updated")
        return 1
    
    print(f"\n✓ Successfully bumped version from {old_version} to {new_version}")
    print(f"✓ Updated {len(updated_files)} file(s)")
    print("\nNext steps:")
    print("  1. Review the changes: git diff")
    print(f"  2. Commit the changes: git commit -am 'Bump version to {new_version}'")
    print(f"  3. Create a tag: git tag v{new_version}")
    print("  4. Push changes and tag: git push && git push --tags")
    
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python bump_version.py <new_version>")
        print("Example: python bump_version.py 0.1.5")
        sys.exit(1)
    
    new_version = sys.argv[1]
    
    # Validate version format
    if not re.match(r'^\d+\.\d+\.\d+', new_version):
        print(f"Error: Invalid version format '{new_version}'")
        print("Expected format: X.Y.Z (e.g., 0.1.5)")
        sys.exit(1)
    
    sys.exit(main(new_version))
