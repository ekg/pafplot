# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

pafplot is a Rust-based command-line tool for generating dotplot visualizations from PAF (Pairwise Alignment Format) files. It rasterizes sequence alignments into PNG images to visualize homology relationships between genomic sequences.

## Key Commands

### Build and Installation
```bash
# Build debug version
cargo build

# Build optimized release version
cargo build --release

# Install to Cargo's bin directory
cargo install --force --path .
```

### Development Commands
```bash
# Run with arguments during development
cargo run -- <paf_file>

# Format code according to Rust style
cargo fmt

# Run Rust linter (if clippy is installed)
cargo clippy

# Run tests (note: no tests currently exist)
cargo test
```

### Usage Example
```bash
# Basic usage - generates output.paf.html
pafplot input.paf

# With custom output and size
pafplot -o output.html -s 2000 input.paf

# Dark theme
pafplot -d input.paf

# Generate PNG output (HTML is default)
pafplot -p input.paf

# BED file overlays
pafplot input.paf --bed regions.bed

# Multiple BED files
pafplot input.paf --bed file1.bed --bed file2.bed

# BEDPE paired region overlays
pafplot input.paf --bedpe pairs.bedpe

# Combine options - HTML viewer with dark theme, custom size, and BED overlays
pafplot -d -s 1500 -o custom.html input.paf --bed regions.bed
```

## Code Architecture

The entire implementation is contained in a single file: `src/main.rs`. Key components:

1. **PAF Parsing**: The `Paf` struct represents alignment records. PAF lines must include CIGAR strings in the `cg:Z:` tag.

2. **BED/BEDPE Support**: BED regions (genomic intervals) and BEDPE regions (paired intervals) can be overlaid on the plot as semi-transparent annotations.

3. **Coordinate Mapping**: Uses `boomphf` for minimal perfect hashing to map sequence names to coordinates efficiently.

4. **Rasterization**: 
   - Creates a 2D pixel buffer based on total sequence lengths
   - Uses Xiaolin Wu's antialiased line algorithm for smooth line drawing
   - Draws grey lines to delimit sequence boundaries

5. **Main Flow**:
   - Parse command-line arguments with `clap`
   - Read and parse PAF file
   - Parse BED/BEDPE files if provided
   - Build coordinate mappings for queries (y-axis) and targets (x-axis)
   - Create raster image and draw alignment matches
   - Encode and save as PNG

## Interactive HTML Viewer

The default output is an interactive HTML file with these features:

- **Pure Canvas Rendering**: Direct alignment rendering without PNG embedding
- **Coordinate Consistency**: Uses identical projection math as PNG generation
- **Level-of-Detail**: Summary view for >100k total view size, detailed view below
- **Interactive Zoom & Pan**: 
  - Scroll wheel to zoom toward pointer location (prevents page scrolling)
  - Left-click drag to pan
  - Right-click drag to zoom-to-box
- **Smart Performance**: Only renders visible elements for smooth interaction
- **Sequence Labels**: Target sequences below plot, query sequences on right (fixed rotation)
- **Real-time Tooltips**: Follow mouse with precise coordinate information
- **Color-Coded Alignments**: Forward/reverse alignments shown in different colors
- **Dark/Light Themes**: Automatically matches the PNG's color scheme
- **BED/BEDPE Overlays**: Semi-transparent genomic region annotations that render behind alignments

The HTML viewer uses vanilla JavaScript and requires no external dependencies.

## Important Notes

- Alignments must include CIGAR strings (produced by `minimap2 -c`, `wfmash`, or `lastz --format=paf:wfmash`)
- Sequences are ordered by length: queries on y-axis, targets on x-axis
- The `-r, --range` zooming parameter has a known bug
- All functionality is in a single source file with no external modules
- HTML output filename automatically derived from input filename (`.paf` → `.paf.html`)
- BED regions render as boundary lines at start/end positions
- BEDPE regions render as 2D rectangles (when implemented)

## Commit Identity

When making git commits for this project, use the identity:
Erik Garrison <erik.garrison@gmail.com>