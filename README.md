# pafplot

An interactive HTML dotplot viewer and rasterizer for sequence alignments with genomic region overlay support

## overview

In the process of generating alignments between whole genomes, we often need to understand the base-level alignment between particular sequences.
`pafplot` creates interactive HTML visualizations and static PNG images from PAF (Pairwise Alignment Format) files.
It renders each alignment match as a line, providing a high-level view of the structure of the alignments and the homology relationships between sequences.

Additionally, `pafplot` supports overlaying genomic regions from BED and BEDPE files as semi-transparent annotations, allowing you to visualize how specific genomic features relate to alignment patterns.

## installation

`pafplot` is built with rust, and so we install using cargo:

```
git clone https://github.com/ekg/pafplot
cd pafplot
cargo install --force --path .
```

## usage

Generate alignments with the cigar string attached to the `cg:Z:` tag.
These can be made by several aligners, including `minimap2 -c`, `wfmash`, or `lastz --format=paf:wfmash`.

### Default HTML output

By default, `pafplot` generates an interactive HTML viewer:

```bash
pafplot aln.paf  # Creates aln.paf.html
```

### Interactive HTML viewer features

The HTML viewer provides a rich, interactive experience with:

- **Zoom and pan**: 
  - Scroll wheel to zoom toward mouse position
  - Left-click and drag to pan around the plot
  - Right-click and drag to zoom into a specific region
- **Sequence labels**: Target sequences displayed below the plot, query sequences on the right
- **Coordinate tooltips**: Real-time coordinate display following the mouse cursor
- **Color-coded alignments**: Forward matches in blue/cyan, reverse matches in red/magenta
- **Performance optimization**: Level-of-detail rendering for smooth interaction with large datasets
- **Dark/light themes**: Automatically matches the theme specified with `-d, --dark`
- **Self-contained**: No external dependencies, works offline

### Command-line options

```bash
# Generate HTML only (default)
pafplot aln.paf

# Generate PNG only
pafplot -p aln.paf

# Generate both HTML and PNG
pafplot -p -h aln.paf

# Specify custom output filename
pafplot -o output.html aln.paf

# Dark theme
pafplot -d aln.paf

# Custom size (affects both PNG and HTML canvas)
pafplot -s 2000 aln.paf

# BED file overlays
pafplot --bed regions.bed aln.paf

# Multiple BED files
pafplot --bed file1.bed --bed file2.bed aln.paf

# BEDPE paired region overlays
pafplot --bedpe pairs.bedpe aln.paf

# Combine all options
pafplot -p -h -d -s 1500 --bed regions.bed --bedpe pairs.bedpe -o custom.html aln.paf
```

### BED and BEDPE overlay support

`pafplot` can overlay genomic regions from BED and BEDPE files as semi-transparent annotations:

#### BED format support
- **Standard BED columns**: `chr start end [name] [score] [strand] [thickStart] [thickEnd] [itemRgb]`
- **Minimum required**: `chr start end`
- **Color support**: Uses `itemRgb` column if present, otherwise generates colors based on region name
- **Rendering**: BED regions appear as filled rectangles spanning the genomic coordinates
- **Multiple files**: Use multiple `--bed` flags to overlay regions from different files

#### BEDPE format support  
- **Standard BEDPE columns**: `chr1 start1 end1 chr2 start2 end2 [name] [score] [strand1] [strand2] [itemRgb]`
- **Minimum required**: `chr1 start1 end1 chr2 start2 end2`
- **Rendering**: BEDPE regions appear as 2D rectangles connecting paired genomic coordinates
- **Use case**: Perfect for visualizing structural variants, Hi-C contacts, or other paired genomic features

#### Example BED/BEDPE usage
```bash
# Overlay gene annotations
pafplot --bed genes.bed genome_alignment.paf

# Overlay repeat regions and structural variants
pafplot --bed repeats.bed --bedpe structural_variants.bedpe alignment.paf

# Complex visualization with custom styling
pafplot -d -s 2000 --bed exons.bed --bed introns.bed --bedpe inversions.bedpe -o detailed.html genome.paf
```

### Output specifications

- **Sequence ordering**: Queries are ordered by length on the y-axis, targets on the x-axis
- **Visual style**: Default is black on white, with grey lines delimiting sequence boundaries
- **Dark theme**: Use `-d, --dark` for white on black rendering

### Complete options reference

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| INPUT | | Input PAF file (required) | - |
| --output | -o | Output filename | input.paf.html |
| --png | -p | Generate PNG image output only | false (HTML default) |
| --html | -h | Generate HTML viewer (use with -p for both) | true when no -p |
| --dark | -d | Use dark theme (white on black) | false |
| --size | -s | Major axis size in pixels | 1000 |
| --range | -r | Zoom to specific range (buggy) | - |
| --bed | | BED file(s) for region overlays (multiple allowed) | none |
| --bedpe | | BEDPE file(s) for paired region overlays (multiple allowed) | none |
| --help | | Display help information | - |

## considerations / bugs

- The `-r, --range` parameter for command-line zooming has a known bug
- The HTML viewer provides full interactive zooming and sequence labels, addressing limitations of the static PNG output
- When focusing on specific sequence pairs, pre-filtering the PAF file to those sequences may provide clearer visualization
- BED/BEDPE overlays render behind alignment data to avoid obscuring the alignments
- Color parsing from BED `itemRgb` fields is planned but not yet implemented
- BEDPE rectangle rendering is planned but not yet implemented (currently only BED regions work)

## author

Erik Garrison

## license

MIT
