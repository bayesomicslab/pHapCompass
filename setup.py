import os
import subprocess
import sys
from pathlib import Path
from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.command.develop import develop
from setuptools.command.install import install


def check_build_dependencies():
    """Check if required build tools are available."""
    missing = []
    
    # Check for gcc or clang
    has_compiler = False
    for compiler in ['gcc', 'clang', 'cc']:
        try:
            subprocess.run([compiler, '--version'], capture_output=True, check=True)
            has_compiler = True
            break
        except (subprocess.CalledProcessError, FileNotFoundError):
            continue
    
    if not has_compiler:
        missing.append('C compiler (gcc, clang, or cc)')
    
    # Check for make
    try:
        subprocess.run(['make', '--version'], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        missing.append('make')
    
    # Check for git (needed for submodules)
    try:
        subprocess.run(['git', '--version'], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        missing.append('git')
    
    return missing


def compile_extracthairs():
    """Compile extractHAIRS from source during installation."""
    print("\n" + "="*70)
    print("Compiling extractHAIRS binary...")
    print("="*70)
    
    # Check build dependencies first
    missing_deps = check_build_dependencies()
    if missing_deps:
        print("\n❌ ERROR: Missing required build tools:")
        for dep in missing_deps:
            print(f"  - {dep}")
        print("\nPlease install the missing dependencies:")
        print("\n  On Ubuntu/Debian:")
        print("    sudo apt-get update")
        print("    sudo apt-get install build-essential git")
        print("\n  On macOS:")
        print("    xcode-select --install")
        print("    # or install via Homebrew: brew install gcc git")
        print("\n  On Fedora/RHEL:")
        print("    sudo dnf install gcc make git")
        print("\nThen retry the installation.")
        print("="*70)
        return False
    
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
        
        # Initialize git submodules (htslib and samtools)
        print("\nInitializing git submodules (htslib and samtools)...")
        try:
            # Check if we're in a git repository
            subprocess.run(['git', 'rev-parse', '--git-dir'], 
                         check=True, capture_output=True, cwd=repo_root)
            
            # Initialize submodules
            result = subprocess.run(['git', 'submodule', 'update', '--init', '--recursive'],
                                  check=True, capture_output=True, text=True, cwd=repo_root)
            print("✓ Submodules initialized successfully")
        except subprocess.CalledProcessError as e:
            print("⚠ Warning: Could not initialize git submodules")
            print("  This is expected if you installed via pip from GitHub")
            print("  The Makefile will attempt to download them automatically")
        
        # Clean previous builds
        subprocess.run(["make", "clean"], check=False, capture_output=True)
        
        # Compile
        print("\nRunning make in third_party/extract_poly...")
        print("This may take a few minutes as it builds htslib and samtools...")
        result = subprocess.run(["make"], check=True, capture_output=True, text=True)
        
        # Show compilation output
        if result.stdout:
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
        
        print(f"\n✓ extractHAIRS compiled successfully!")
        print(f"  Binary location: {binary_dest}")
        print(f"  Binary size: {binary_dest.stat().st_size / 1024:.1f} KB")
        print("="*70)
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"\n❌ ERROR: Failed to compile extractHAIRS")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        print("\n" + "="*70)
        print("INSTALLATION WILL CONTINUE, but extractHAIRS will not be available.")
        print("\nTo use pHapCompass, you must either:")
        print("  1. Install extractHAIRS separately:")
        print("     cd third_party/extract_poly")
        print("     make")
        print("     # Add build/extractHAIRS to your PATH")
        print("\n  2. Use --frag-path with precomputed fragment files")
        print("="*70)
        return False
    except Exception as e:
        print(f"⚠ WARNING: Could not compile extractHAIRS: {e}")
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
