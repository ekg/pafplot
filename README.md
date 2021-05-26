# pafplot

A base-level sequence alignment rasterizer / dotplot generator

## overview

In the process of generating alignments between whole genomes, we often need to understand the base-level alignment between particular sequences.
`pafplot` allows us to do so by rasterizing the matches alignment set.
It draws a line on a raster image to represent each match found in a set of alignments.
The resulting image provides a high-level view of the structure of the alignments, and in consequence the homology relationships between the sequences in consideration.

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
With alignments in `aln.paf`, we would render a plot in `aln.paf.png` using this call:

```
pafplot aln.paf
```

The scale of the plot can be changed with the `-s, --size` parameter.
We supply the number of pixels in the longest axis of the rendering.
It's also possible to specify a different output raster image with `-p, --png`.

Queries are ordered by length on the y-axis, while targets are ordered by length on the x-axis.
The plot is white on black, with grey lines delinating sequence boundaries.
These colors may be inverted for a dark theme by adding the `-d, --dark` parameter.

## considerations / bugs

`pafplot` should support zooming functionality, but currently the `-r, --range` parameter has a bug.
The raster image would benefit from sequence name labels, but these are currently not implemented.
As the order is given by the input sequence lengths, it is possible to manually determine the query and target identities.
But, when considering a single pair of sequences, it is often easiest to filter the input alignments down to those between the pair of sequences in question.

## author

Erik Garrison

## license

MIT
