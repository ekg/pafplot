# pafplot

An interactive HTML dotplot viewer and rasterizer for sequence alignments

## overview

In the process of generating alignments between whole genomes, we often need to understand the base-level alignment between particular sequences.
`pafplot` creates interactive HTML visualizations and static PNG images from PAF (Pairwise Alignment Format) files.
It renders each alignment match as a line, providing a high-level view of the structure of the alignments and the homology relationships between sequences.

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

# Generate both HTML and PNG
pafplot -p aln.paf

# Specify custom output filename
pafplot -o output.html aln.paf

# Dark theme
pafplot -d aln.paf

# Custom size (affects both PNG and HTML canvas)
pafplot -s 2000 aln.paf

# Combine options
pafplot -p -d -s 1500 -o custom.html aln.paf
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
| --png | -p | Generate PNG image output | false (HTML only) |
| --dark | -d | Use dark theme (white on black) | false |
| --size | -s | Major axis size in pixels | 1000 |
| --range | -r | Zoom to specific range (buggy) | - |
| --html | -h | Generate HTML viewer | true (default) |
| --help | | Display help information | - |

## considerations / bugs

- The `-r, --range` parameter for command-line zooming has a known bug
- The HTML viewer provides full interactive zooming and sequence labels, addressing limitations of the static PNG output
- When focusing on specific sequence pairs, pre-filtering the PAF file to those sequences may provide clearer visualization

## author

Erik Garrison

## license

MIT
