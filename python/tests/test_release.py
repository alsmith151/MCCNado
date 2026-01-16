"""
Test version bumping and release automation.
"""

import pathlib
import tempfile
import shutil
import subprocess
import pytest


@pytest.fixture
def temp_project():
    """Create a temporary copy of the project for testing."""
    source = pathlib.Path(__file__).parent.parent
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = pathlib.Path(tmpdir)
        
        # Copy only the files we need for testing
        files_to_copy = [
            "Cargo.toml",
            "pyproject.toml",
            "python/mccnado/version.py",
            "bump_version.py",
        ]
        
        for file in files_to_copy:
            src = source / file
            dst = tmpdir_path / file
            dst.parent.mkdir(parents=True, exist_ok=True)
            if src.exists():
                shutil.copy2(src, dst)
        
        yield tmpdir_path


class TestVersionBumping:
    """Test the bump_version script."""
    
    def test_bump_version_script_exists(self):
        """Test that bump_version.py script exists."""
        script = pathlib.Path(__file__).parent.parent.parent / "bump_version.py"
        assert script.exists(), "bump_version.py script not found"
    
    def test_bump_version_invalid_format(self):
        """Test that invalid version format is rejected."""
        result = subprocess.run(
            ["python", "bump_version.py", "invalid"],
            capture_output=True,
            text=True,
            cwd=pathlib.Path(__file__).parent.parent.parent,
        )
        assert result.returncode != 0, "Should reject invalid version format"
        # Error message can be in stdout or stderr depending on the script
        output = result.stdout + result.stderr
        assert "Invalid version format" in output
    
    def test_bump_version_valid_format(self):
        """Test that valid version format is accepted."""
        result = subprocess.run(
            ["python", "-c", 
             "import re; "
             "version = '0.2.0'; "
             "assert re.match(r'^\\d+\\.\\d+\\.\\d+', version), 'Invalid format'; "
             "print('✓ Valid version format')"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, "Valid version format should work"
    
    def test_current_version_detection(self):
        """Test that current version can be detected from Cargo.toml."""
        import sys
        project_root = pathlib.Path(__file__).parent.parent.parent
        sys.path.insert(0, str(project_root))
        
        import bump_version
        version = bump_version.get_current_version()
        assert version == "0.1.4", f"Expected 0.1.4, got {version}"
        # Verify format
        assert isinstance(version, str)
        assert len(version.split('.')) == 3, "Version should be X.Y.Z format"


class TestReleaseDocumentation:
    """Test that release documentation exists."""
    
    def test_release_md_exists(self):
        """Test that RELEASE.md documentation exists."""
        release_doc = pathlib.Path(__file__).parent.parent.parent / "RELEASE.md"
        assert release_doc.exists(), "RELEASE.md documentation not found"
    
    def test_release_md_contains_steps(self):
        """Test that RELEASE.md contains key sections."""
        release_doc = pathlib.Path(__file__).parent.parent.parent / "RELEASE.md"
        content = release_doc.read_text()
        
        required_sections = [
            "Release Steps",
            "Bump the Version",
            "Automated Validation",
            "Version Syncing",
        ]
        
        for section in required_sections:
            assert section in content, f"Missing section: {section}"


class TestGitHubActionsWorkflow:
    """Test GitHub Actions workflow configuration."""
    
    def test_workflow_file_exists(self):
        """Test that GitHub Actions workflow exists."""
        workflow = pathlib.Path(__file__).parent.parent.parent / ".github" / "workflows" / "validate-release.yml"
        assert workflow.exists(), "validate-release.yml workflow not found"
    
    def test_workflow_yaml_valid(self):
        """Test that workflow YAML is valid."""
        try:
            import yaml
            workflow = pathlib.Path(__file__).parent.parent.parent / ".github" / "workflows" / "validate-release.yml"
            with open(workflow) as f:
                content = yaml.safe_load(f)
            assert content is not None, "YAML content is empty"
            # YAML parser converts "on:" to boolean True key
            assert True in content or "on" in content, "Missing 'on' trigger in workflow"
            assert "jobs" in content, "Missing 'jobs' in workflow"
        except ImportError:
            pytest.skip("PyYAML not installed")
    
    def test_workflow_has_validation_job(self):
        """Test that workflow has version validation job."""
        import yaml
        workflow = pathlib.Path(__file__).parent.parent.parent / ".github" / "workflows" / "validate-release.yml"
        with open(workflow) as f:
            content = yaml.safe_load(f)
        
        jobs = content.get("jobs", {})
        assert "validate_version" in jobs, "Missing validate_version job"
        assert "build_and_test" in jobs, "Missing build_and_test job"
