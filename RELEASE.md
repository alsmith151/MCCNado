# Release Process

This document describes how to create a new release of MCCNado.

## Version Syncing

MCCNado maintains version consistency across:
- `Cargo.toml` (Rust package version)
- `pyproject.toml` (Python package metadata)
- `python/mccnado/version.py` (Fallback version)

The source of truth is **`Cargo.toml`**. When maturin builds the Python package, it automatically extracts the version from `Cargo.toml`.

## Release Steps

### 1. Bump the Version

Use the provided script to update all version files:

```bash
python bump_version.py 0.1.5
```

This script will:
- Update `Cargo.toml`
- Update `pyproject.toml`
- Update `python/mccnado/version.py`
- Provide instructions for the next steps

### 2. Review Changes

```bash
git diff
```

Verify that all version files have been updated correctly.

### 3. Commit Version Bump

```bash
git commit -am "Bump version to 0.1.5"
```

### 4. Create Git Tag

```bash
git tag v0.1.5
```

**Important**: The tag must be in the format `v<VERSION>` where `<VERSION>` matches the version in `Cargo.toml`.

### 5. Push Changes and Tag

```bash
git push
git push --tags
```

## Automated Validation

When you push a tag (e.g., `v0.1.5`), the GitHub Actions workflow `.github/workflows/validate-release.yml` will:

1. **Validate** that the tag version matches `Cargo.toml` version
2. **Build** the package with maturin
3. **Run tests** to ensure everything works
4. **Verify** version consistency across Python and Rust

If the tag version doesn't match `Cargo.toml`, the workflow will fail with a clear error message.

## Example Release Workflow

```bash
# Bump version
python bump_version.py 0.1.5

# Review changes
git diff

# Commit
git commit -am "Bump version to 0.1.5"

# Create tag
git tag v0.1.5

# Push
git push && git push --tags

# Wait for GitHub Actions validation...
```

## Troubleshooting

### "Version mismatch" error in GitHub Actions

If you see a version mismatch error:

1. Check the tag name: `git tag -l`
2. Check Cargo.toml version: `grep '^version' Cargo.toml`
3. Ensure they match (tag should be `v<VERSION>` where `<VERSION>` is from Cargo.toml)

### To fix a tag

If you created a tag with the wrong version:

```bash
# Delete the tag locally
git tag -d v0.1.4

# Delete the tag remotely
git push origin --delete v0.1.4

# Create the correct tag
git tag v0.1.5
git push --tags
```

## Version Numbering

MCCNado follows [Semantic Versioning](https://semver.org/):

- **MAJOR.MINOR.PATCH** (e.g., 0.1.4)
- Increment MAJOR for incompatible changes
- Increment MINOR for backward-compatible features
- Increment PATCH for bug fixes
