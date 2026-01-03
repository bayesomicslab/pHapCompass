import os
import subprocess
import sys
from pathlib import Path
from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.command.develop import develop
from setuptools.command.install import install


def compile_extracthairs():
    """Compile extractHAIRS from source during installation."""
    print("\n" + "="*70)
    print("Compiling extractHAIRS binary...")
    print("="*70)
    
    # Get paths
    repo_root = Path(__file__).parent.resolve()
    extract_poly_dir = repo_root / "third_party" / "extract_poly"
    bin_dir = repo_root / "src" / "phapcompass" / "bin"
    
    if not extract_poly_dir.exists():
        print(f"WARNING: {extract_poly_dir} not found. Skipping extractHAIRS compilation.")
        print("Users will need to install extractHAIRS separately or use --frag-path.")
        return False
    
    # Create bin directory
    bin_dir.mkdir(parents=True, exist_ok=True)
    
    # Compile
    try:
        original_dir = os.getcwd()
        os.chdir(extract_poly_dir)
        
        # Clean previous builds
        subprocess.run(["make", "clean"], check=False, capture_output=True)
        
        # Compile
        print("Running make in third_party/extract_poly...")
        result = subprocess.run(["make"], check=True, capture_output=True, text=True)
        print(result.stdout)
        
        os.chdir(original_dir)
        
        # Find the compiled binary
        build_dir = extract_poly_dir / "build"
        binary_source = build_dir / "extractHAIRS"
        binary_dest = bin_dir / "extractHAIRS"
        
        if not binary_source.exists():
            print(f"ERROR: Compilation succeeded but binary not found at {binary_source}")
            return False
        
        # Copy to package bin directory
        import shutil
        shutil.copy2(binary_source, binary_dest)
        binary_dest.chmod(0o755)
        
        print(f"✓ extractHAIRS compiled and copied to {binary_dest}")
        print(f"  Binary size: {binary_dest.stat().st_size / 1024:.1f} KB")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"\nERROR: Failed to compile extractHAIRS")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        print("\nWARNING: Installation will continue, but extractHAIRS will not be available.")
        print("Users must either:")
        print("  1. Install extractHAIRS separately and add to PATH")
        print("  2. Use --frag-path with precomputed fragment files")
        return False
    except Exception as e:
        print(f"WARNING: Could not compile extractHAIRS: {e}")
        print("Installation will continue without it.")
        return False


class CustomBuildPy(build_py):
    """Custom build_py that compiles extractHAIRS."""
    def run(self):
        compile_extracthairs()
        super().run()


class CustomDevelop(develop):
    """Custom develop that compiles extractHAIRS."""
    def run(self):
        compile_extracthairs()
        super().run()


class CustomInstall(install):
    """Custom install that compiles extractHAIRS."""
    def run(self):
        compile_extracthairs()
        super().run()


if __name__ == "__main__":
    setup(
        cmdclass={
            'build_py': CustomBuildPy,
            'develop': CustomDevelop,
            'install': CustomInstall,
        },
    )


