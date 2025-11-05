# Travel Ratio Computation Pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4383281.svg)](http://doi.org/10.5281/zenodo.4383281)

A refactored and improved R pipeline for computing Travel Ratio (TR) from ChIP-seq signal tracks.

## Overview

**Travel Ratio** is a metric that measures transcriptional elongation efficiency by computing the ratio of gene body signal to TSR (Transcription Start Region) signal in ChIP-seq data. Higher TR values indicate more efficient elongation.

### Version 2.0 Improvements

This refactored version includes major improvements over the original:

✅ **Modular code architecture** - Reusable functions with clear responsibilities
✅ **Configuration file support** - YAML-based parameter management
✅ **Comprehensive error handling** - Validates inputs and provides informative messages
✅ **Native R implementation** - Replaced Python and bash code with R functions
✅ **Better documentation** - Roxygen2-style function docs and detailed comments
✅ **Enhanced visualizations** - Multiple plot types with publication-quality styling
✅ **Logging and progress tracking** - Real-time progress messages and log files
✅ **Reproducibility** - Session info, parameter logging, and structured outputs
✅ **Flexible execution** - Both interactive (Rmd) and command-line (R script) modes

---

## Installation

### Prerequisites

**R Packages** (automatically installed when running the script):
- `yaml` - Configuration file parsing
- `dplyr`, `tidyr` - Data manipulation
- `ggplot2` - Visualization
- `readr`, `stringr` - Data I/O and string operations
- `optparse` - Command-line argument parsing
- `logger` - Logging functionality

**External Tools** (must be installed separately):
- [bedtools](https://bedtools.readthedocs.io/) - BED file manipulation
- [bigWigAverageOverBed](http://hgdownload.soe.ucsc.edu/admin/exe/) - Signal computation over genomic regions

### Installing External Tools

**On Ubuntu/Debian:**
```bash
# Install bedtools
sudo apt-get install bedtools

# Install bigWigAverageOverBed
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed
chmod +x bigWigAverageOverBed
sudo mv bigWigAverageOverBed /usr/local/bin/
```

**On macOS (with Homebrew):**
```bash
# Install bedtools
brew install bedtools

# Install UCSC tools (includes bigWigAverageOverBed)
brew install --cask ucsc-genome-browser
```

**Using conda:**
```bash
conda install -c bioconda bedtools ucsc-bigwigaverageoverbed
```

---

## Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/yourusername/travelRatio.git
cd travelRatio
```

### 2. Configure Your Analysis

Edit `config.yaml` to specify your input files and parameters:

```yaml
input:
  homer_peak_file: "path/to/your/genes.bed"
  selection_bigwig: "path/to/control.bigwig"
  sample_bigwig: "path/to/sample.bigwig"

output:
  output_dir: "output"
  prefix: "my_sample"
```

### 3. Run the Pipeline

**Option A: Command-line execution (recommended for batch processing)**

```bash
Rscript compute_travel_ratio.R --config config.yaml
```

**Option B: Interactive R Markdown (recommended for exploration)**

Open `computeTR_refactored.Rmd` in RStudio and knit the document, or run:

```r
rmarkdown::render("computeTR_refactored.Rmd")
```

---

## Usage

### Command-Line Options

```bash
Rscript compute_travel_ratio.R [options]

Options:
  -c, --config CONFIG     Path to configuration YAML file [default: config.yaml]
  -o, --output DIR        Output directory (overrides config file)
  -v, --verbose           Print detailed progress messages
  -h, --help              Show this help message
```

**Examples:**

```bash
# Use default config.yaml
Rscript compute_travel_ratio.R

# Specify custom config file
Rscript compute_travel_ratio.R --config my_config.yaml

# Override output directory
Rscript compute_travel_ratio.R --config config.yaml --output results_2025

# Verbose mode for debugging
Rscript compute_travel_ratio.R --verbose
```

### Configuration File

The `config.yaml` file contains all configurable parameters. Key sections:

#### Input Files
```yaml
input:
  homer_peak_file: "refGene_hg19_genebody.new.bed"
  selection_bigwig: "control.signal.bigwig"
  sample_bigwig: "sample.fc.signal.bw"
```

#### Gene Filtering
```yaml
gene_filters:
  min_gene_length: 1000      # Minimum gene length (bp)
  distance_cutoff: 1000      # Minimum distance between genes (bp)
```

#### TSR Parameters
```yaml
tsr_parameters:
  upstream_extension: 50     # bp upstream of TSS
  downstream_extension: 300  # bp downstream of TSS
```

#### Transcript Selection
```yaml
transcript_selection:
  signal_cutoff: 0.25        # Minimum mean signal for selection
  selection_method: "highest"
```

#### Gene Body Extension
```yaml
genebody_parameters:
  auto_calculate_extension: TRUE  # Auto-calculate based on distance_cutoff
  manual_extension: 500           # Manual value (if auto is FALSE)
```

See `config.yaml` for all available options.

---

## Output Files

The pipeline creates the following directory structure:

```
output/
├── bed_files/
│   ├── travel_ratio_genes_gt1000bp.bed
│   ├── travel_ratio_unique_genes.bed
│   ├── travel_ratio_nonoverlapping_genes.bed
│   ├── travel_ratio_tsr_final.bed
│   └── travel_ratio_genebody_final.bed
├── signals/
│   ├── travel_ratio_selection_signal.txt
│   ├── travel_ratio_unique_transcripts.txt
│   ├── travel_ratio_tsr_signal.txt
│   └── travel_ratio_genebody_signal.txt
├── plots/
│   ├── travel_ratio_selection_signal_cdf.png
│   ├── travel_ratio_travel_ratio_cdf.png
│   ├── travel_ratio_travel_ratio_hist.png
│   └── travel_ratio_tsr_vs_gb_scatter.png
├── logs/
│   └── travel_ratio.log
└── travel_ratio_travel_ratio.txt  # MAIN OUTPUT
```

### Main Output File

`travel_ratio_travel_ratio.txt` contains two columns:
1. Gene/transcript name
2. Travel Ratio value (Gene Body signal / TSR signal)

**Format:**
```
transcript1!gene1!info    2.456
transcript2!gene2!info    1.234
transcript3!gene3!info    3.789
...
```

---

## Pipeline Steps

The pipeline performs the following steps:

1. **Gene Filtering** - Convert gene annotations to BED format and filter by length (≥1kb)
2. **TSR Generation** - Define TSR regions around transcription start sites
3. **Signal Computation** - Compute ChIP-seq signal over genes for transcript selection
4. **Transcript Selection** - For each gene, select the transcript with highest signal
5. **Overlap Removal** - Remove genes that overlap or are too close to each other
6. **Final Region Definition** - Generate final TSR and gene body BED files
7. **Signal Quantification** - Compute signals over TSR and gene body regions
8. **TR Calculation** - Calculate Travel Ratio (GB signal / TSR signal)
9. **Visualization** - Generate diagnostic plots and summary statistics

---

## Biological Interpretation

### What is Travel Ratio?

Travel Ratio measures the efficiency of transcriptional elongation:

- **High TR (e.g., > 2)**: Strong gene body signal relative to TSR, indicating efficient elongation
- **Low TR (e.g., < 1)**: Weak gene body signal relative to TSR, indicating paused or stalled polymerase
- **TR ≈ 1**: Similar signal in TSR and gene body

### Applications

Travel Ratio is useful for studying:
- RNA Polymerase II elongation dynamics
- Effects of elongation factors or inhibitors
- Differences between paused and actively elongating genes
- Transcriptional regulation during development or stress responses

### References

If you use this pipeline, please cite:

```
Xu, X. M., & Li, J. (2020). travelRatio: A R notebook to compute travel ratio
using ChIP-seq signal track (Version 1.0). Zenodo.
https://doi.org/10.5281/zenodo.4383281
```

---

## Troubleshooting

### Common Issues

**1. "Tool not found" error**

Ensure bedtools and bigWigAverageOverBed are installed and in your PATH:
```bash
which bedtools
which bigWigAverageOverBed
```

If not in PATH, specify full paths in `config.yaml`:
```yaml
tools:
  bedtools: "/path/to/bedtools"
  bigwig_average_over_bed: "/path/to/bigWigAverageOverBed"
```

**2. "Input file not found" error**

Check that all input file paths in `config.yaml` are correct and files exist.

**3. "No genes pass filtering" warning**

Try adjusting these parameters in `config.yaml`:
- Lower `min_gene_length` (e.g., 500 instead of 1000)
- Lower `signal_cutoff` (e.g., 0.1 instead of 0.25)
- Check that your bigWig files have actual signal data

**4. Memory issues with large datasets**

For very large genomes or many genes:
- Process chromosomes separately
- Increase system memory
- Use a compute cluster

**5. R package installation fails**

Try installing packages manually:
```r
install.packages(c("yaml", "dplyr", "tidyr", "ggplot2",
                   "readr", "stringr", "optparse", "logger"))
```

---

## Advanced Usage

### Running on a Compute Cluster

Example SLURM submission script:

```bash
#!/bin/bash
#SBATCH --job-name=travel_ratio
#SBATCH --output=tr_%j.log
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

module load R/4.3
module load bedtools

Rscript compute_travel_ratio.R --config config.yaml --verbose
```

### Processing Multiple Samples

Create a wrapper script to process multiple samples:

```bash
#!/bin/bash

for sample in sample1 sample2 sample3; do
  echo "Processing ${sample}..."

  # Update config with sample-specific paths
  sed "s/SAMPLE_NAME/${sample}/g" config_template.yaml > config_${sample}.yaml

  # Run pipeline
  Rscript compute_travel_ratio.R \
    --config config_${sample}.yaml \
    --output output_${sample}
done
```

### Custom Analysis

To use individual functions in your own R scripts:

```r
# Source the functions
source("compute_travel_ratio.R")

# Load config
config <- yaml::read_yaml("config.yaml")

# Use individual functions
bed_filtered <- convert_and_filter_peaks(
  "genes.bed",
  "output/filtered.bed",
  min_length = 1000
)

# Compute signals
compute_signal_over_bed(
  "sample.bw",
  bed_filtered,
  "output/signals.txt"
)
```

---

## Comparison with Original Version

| Feature | Original | Refactored v2.0 |
|---------|----------|-----------------|
| Parameterization | Hardcoded | YAML config file |
| Language | R + bash + Python | Pure R |
| Error handling | None | Comprehensive |
| Code organization | Linear script | Modular functions |
| Documentation | Minimal | Extensive |
| Logging | None | Full logging system |
| Plots | Basic | Publication-quality |
| Command-line support | No | Yes |
| Reproducibility | Limited | Full |

---

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/improvement`)
3. Make your changes
4. Add tests if applicable
5. Commit your changes (`git commit -am 'Add new feature'`)
6. Push to the branch (`git push origin feature/improvement`)
7. Create a Pull Request

---

## License

This project maintains the same license as the original work. Please check with the original authors for licensing details.

---

## Contact

**Original Authors:**
- Xiaojun Max Xu
- Jinting Li

**For questions about the refactored version:**
- Open an issue on GitHub
- Check the documentation in the code

---

## Changelog

### Version 2.0 (2025)
- Complete refactor with modular architecture
- Added YAML configuration support
- Replaced Python code with native R implementation
- Added comprehensive error handling and validation
- Enhanced visualizations and logging
- Added command-line interface
- Improved documentation

### Version 1.0 (2020)
- Original implementation
- Basic Travel Ratio computation
- Initial release on Zenodo
