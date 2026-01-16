"""
Version information for MCCNado.

The version is automatically synced from Cargo.toml by maturin during the build.
"""

try:
    # Try to get version from package metadata (populated by maturin from Cargo.toml)
    from importlib.metadata import version as _get_version
    __version__ = _get_version("mccnado")
except Exception:
    # Fallback version if package metadata is not available
    # This should match the version in Cargo.toml
    __version__ = "0.1.6"

__all__ = ["__version__"]
