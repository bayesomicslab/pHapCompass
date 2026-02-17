# pHapCompass Installation Troubleshooting Guide

This guide addresses common installation issues with pHapCompass, particularly related to compiling the extractHAIRS dependency.

## Understanding the extractHAIRS Dependency

pHapCompass uses a modified version of extractHAIRS (from HAPCUT2) to extract fragment information from BAM files. This tool must be compiled from C source code, which requires:

1. A C compiler (gcc, clang, or cc)
2. The make build tool
3. Git (for downloading required submodules: htslib and samtools)

## Common Installation Issues

### Issue 1: "gcc: command not found" or "make: command not found"

**Cause:** Required build tools are not installed on your system.

**Solution:**

**On Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install build-essential git
```

**On macOS:**
```bash
# Install Xcode Command Line Tools
xcode-select --install

# Verify installation
gcc --version
make --version
```

**On Fedora/RHEL/CentOS:**
```bash
sudo dnf install gcc gcc-c++ make git
```

**On Arch Linux:**
```bash
sudo pacman -S base-devel git
```

### Issue 2: Git submodule errors

**Symptoms:**
- "fatal: No url found for submodule path"
- "Makefile:XX: sub/htslib/Makefile: No such file or directory"

**Cause:** The git submodules (htslib and samtools) were not downloaded.

**Solution:**

**If you cloned the repository:**
```bash
cd pHapCompass
git submodule update --init --recursive
```

**If you installed via pip:**
```bash
# Uninstall first
pip uninstall phapcompass

# Clone properly with submodules
git clone --recursive https://github.com/bayesomicslab/pHapCompass.git
cd pHapCompass
pip install -e .
```

### Issue 3: Compilation fails with htslib/samtools errors

**Symptoms:**
- Errors mentioning "htslib.h" or "sam.h" not found
- Linker errors with "-lhts" or "-lbam"

**Cause:** The htslib/samtools libraries failed to compile.

**Solution:**

```bash
cd pHapCompass/third_party/extract_poly

# Check if submodules exist
ls -la sub/

# If they don't exist, initialize them
git submodule update --init --recursive

# Clean and rebuild everything
make clean
make

# If successful, the binary should exist at:
ls -lh build/extractHAIRS
```

### Issue 4: extractHAIRS binary not found at runtime

**Symptoms:**
```
FileNotFoundError: extractHAIRS binary not found
```

**Cause:** The binary was not compiled or not copied to the correct location.

**Solution:**

**Check if binary exists:**
```bash
# From pHapCompass root directory
ls -lh src/phapcompass/bin/extractHAIRS
```

**If not found, manually compile and copy:**
```bash
cd third_party/extract_poly
make clean
make

# Verify binary was created
ls -lh build/extractHAIRS

# Copy to package bin directory
mkdir -p ../../src/phapcompass/bin
cp build/extractHAIRS ../../src/phapcompass/bin/
chmod +x ../../src/phapcompass/bin/extractHAIRS
```

### Issue 5: Permission errors during compilation

**Symptoms:**
- "Permission denied" when running make
- Cannot create build directory

**Solution:**

**Check directory permissions:**
```bash
ls -ld third_party/extract_poly
```

**Fix permissions if needed:**
```bash
chmod -R u+w third_party/extract_poly
```

**If installed in system Python:**
```bash
# Try installing in user mode instead
pip install --user -e .
```

### Issue 6: Platform-specific compilation issues

**macOS: "ld: library not found for -lz"**

**Solution:**
```bash
# Install zlib via Homebrew
brew install zlib

# Export library path
export LDFLAGS="-L/usr/local/opt/zlib/lib"
export CPPFLAGS="-I/usr/local/opt/zlib/include"

# Retry compilation
cd third_party/extract_poly
make clean
make
```

**macOS: "clang: error: unsupported option '-fcommon'"**

**Solution:**
This is an older flag. Edit `third_party/extract_poly/Makefile` and change line 5:
```makefile
# FROM:
CC=gcc -g -O3 -Wall -D_GNU_SOURCE -fcommon

# TO:
CC=gcc -g -O3 -Wall -D_GNU_SOURCE
```

**Linux: Missing development libraries**

**Solution:**
```bash
# Ubuntu/Debian
sudo apt-get install zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev

# Fedora/RHEL
sudo dnf install zlib-devel bzip2-devel xz-devel libcurl-devel
```

## Alternative: Using Pre-computed Fragment Files

If compilation continues to fail and you have fragment files (`.frag`), you can skip extractHAIRS:

```bash
# Install without extractHAIRS
pip install git+https://github.com/bayesomicslab/pHapCompass.git

# Use with fragment files
phapcompass --data-type short \
  --frag-path precomputed.frag \
  --vcf-path sample.vcf.gz \
  --result-path output.vcf.gz
```

## Alternative: Manual extractHAIRS Installation

You can also install the original HAPCUT2 separately:

```bash
# Clone HAPCUT2
git clone --recursive https://github.com/vibansal/HapCUT2.git
cd HapCUT2
make

# Add to PATH
export PATH=$PATH:$(pwd)/build

# Now pHapCompass will use the system extractHAIRS
```

## Verifying Successful Installation

After resolving issues, verify everything works:

```bash
# Test pHapCompass installation
phapcompass --help

# Test extractHAIRS specifically
which extractHAIRS
# OR
python -c "from phapcompass.utils import find_extracthairs; print(find_extracthairs())"

# Run a test with example data
phapcompass --data-type short \
  --bam-path test_data/short_data_example/0.bam \
  --vcf-path test_data/ref_example/Chr1_unphased.vcf \
  --result-path test_output.vcf.gz
```

## Still Having Issues?

If none of these solutions work:

1. **Collect diagnostic information:**
```bash
# System information
uname -a
gcc --version
make --version
git --version
python --version

# Package information
pip show phapcompass
ls -lh src/phapcompass/bin/
```

2. **Open a GitHub issue** with:
   - Your operating system and version
   - The complete error message
   - Output from the diagnostic commands above
   - Steps you've already tried

3. **Ask on the GitHub Discussions** page for community support
