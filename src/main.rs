use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::path::Path;
//use std::cmp;

use boomphf::*;

use rgb::*;
extern crate line_drawing;
use line_drawing::XiaolinWu;

use itertools::Itertools;

extern crate clap;
use clap::{App, Arg}; //, SubCommand};

extern crate base64;

fn for_each_line_in_file(paf_filename: &str, mut callback: impl FnMut(&str)) {
    let file = File::open(paf_filename).unwrap();
    let reader: Box<dyn BufRead> = if paf_filename.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    for line in reader.lines() {
        callback(&line.unwrap());
    }
}

#[derive(Debug, Clone)]
struct AlignedSeq {
    // name of the given sequence
    name: String,
    // its length
    length: usize,
    // its rank among other seqs in the query or target set
    rank: usize,
    // its start offset in the global all-to-all alignment matrix
    offset: usize,
}

impl AlignedSeq {
    fn new() -> Self {
        AlignedSeq {
            name: String::new(),
            length: 0,
            rank: 0,
            offset: 0,
        }
    }
}

#[derive(Debug, Clone)]
struct BedRegion {
    chr: String,
    start: usize,
    end: usize,
    name: Option<String>,
    color: Option<String>,
}

#[derive(Debug, Clone)]
struct BedpeRegion {
    chr1: String,
    start1: usize,
    end1: usize,
    chr2: String,
    start2: usize,
    end2: usize,
    name: Option<String>,
    color: Option<String>,
}

struct PafFile {
    // our input file
    filename: String,
    // each target name in order of first appearance
    targets: Vec<AlignedSeq>,
    // each query name in order of first appearance
    queries: Vec<AlignedSeq>,
    // maps from sequence name to internal id
    target_mphf: Mphf<String>,
    // maps from sequence name to internal id
    query_mphf: Mphf<String>,
    // target axis length
    target_length: f64,
    // query axis length
    query_length: f64,
    // optional PAF tag used to color alignment lines (e.g. "bi", "gi", "dv")
    identity_tag: Option<String>,
}

fn paf_query(line: &str) -> String {
    line.split('\t').next().unwrap().into()
}

fn paf_query_length(line: &str) -> usize {
    line.split('\t').nth(1).unwrap().parse::<usize>().unwrap()
}

fn paf_query_begin(line: &str) -> usize {
    line.split('\t').nth(2).unwrap().parse::<usize>().unwrap()
}

fn paf_query_end(line: &str) -> usize {
    line.split('\t').nth(3).unwrap().parse::<usize>().unwrap()
}

fn paf_query_is_rev(line: &str) -> bool {
    line.split('\t').nth(4).unwrap() == "-"
}

fn paf_target(line: &str) -> String {
    line.split('\t').nth(5).unwrap().into()
}

fn paf_target_length(line: &str) -> usize {
    line.split('\t').nth(6).unwrap().parse::<usize>().unwrap()
}

fn paf_target_begin(line: &str) -> usize {
    line.split('\t').nth(7).unwrap().parse::<usize>().unwrap()
}

fn paf_target_end(line: &str) -> usize {
    line.split('\t').nth(8).unwrap().parse::<usize>().unwrap()
}

/// BLAST identity from the PAF fixed fields: matches / alignment block length.
fn paf_blast_identity(line: &str) -> f64 {
    let matches: f64 = line
        .split('\t')
        .nth(9)
        .and_then(|s| s.parse().ok())
        .unwrap_or(0.0);
    let aln_len: f64 = line
        .split('\t')
        .nth(10)
        .and_then(|s| s.parse().ok())
        .unwrap_or(1.0);
    if aln_len > 0.0 {
        matches / aln_len
    } else {
        0.0
    }
}

/// Value of a numeric PAF tag (e.g. bi:f:, gi:f:, dv:f:) as f64, if present.
fn paf_tag_f64(line: &str, tag: &str) -> Option<f64> {
    for field in line.split('\t').skip(11) {
        let mut it = field.splitn(3, ':');
        if it.next() == Some(tag) {
            let ty = it.next().unwrap_or("");
            let val = it.next().unwrap_or("");
            if ty == "i" || ty == "f" || ty == "d" || ty == "u" {
                return val.parse::<f64>().ok();
            }
        }
    }
    None
}

/// Identity-style score for a PAF line: the given tag when present, else
/// BLAST identity from the fixed fields.
fn paf_line_identity(line: &str, tag: &str) -> f64 {
    paf_tag_f64(line, tag).unwrap_or_else(|| paf_blast_identity(line))
}

/// Turbo colormap (Google's polynomial approximation), t in 0..1.
fn turbo(t: f64) -> RGB8 {
    let t = t.clamp(0.0, 1.0);
    let poly = |c: [f64; 6]| -> f64 {
        (c[0]
            + c[1] * t
            + c[2] * t.powi(2)
            + c[3] * t.powi(3)
            + c[4] * t.powi(4)
            + c[5] * t.powi(5))
        .clamp(0.0, 1.0)
    };
    let r = poly([
        0.13572138,
        4.61539230,
        -42.66032258,
        132.13108234,
        -152.94239396,
        59.28637943,
    ]);
    let g = poly([
        0.09140261,
        2.19418839,
        4.84296658,
        -14.18503333,
        4.27729857,
        2.82956604,
    ]);
    let b = poly([
        0.10667430,
        12.64194608,
        -60.58204836,
        110.36276771,
        -89.90310912,
        27.31224958,
    ]);
    RGB8 {
        r: (r * 255.0).round() as u8,
        g: (g * 255.0).round() as u8,
        b: (b * 255.0).round() as u8,
    }
}

impl PafFile {
    fn new(filename: &str) -> Self {
        let mut query_names: Vec<String> = Vec::new();
        let mut target_names: Vec<String> = Vec::new();
        for_each_line_in_file(filename, |l: &str| {
            query_names.push(paf_query(l));
            target_names.push(paf_target(l));
        });
        query_names = query_names.into_iter().sorted().dedup().collect();
        target_names = target_names.into_iter().sorted().dedup().collect();
        let query_mphf = Mphf::new(1.7, &query_names);
        let target_mphf = Mphf::new(1.7, &target_names);
        let mut seen_queries = vec![false; query_names.len()];
        let mut seen_targets = vec![false; target_names.len()];
        let mut queries: Vec<AlignedSeq> = vec![AlignedSeq::new(); query_names.len()];
        let mut targets: Vec<AlignedSeq> = vec![AlignedSeq::new(); target_names.len()];
        for_each_line_in_file(filename, |l: &str| {
            let query_name: String = paf_query(l);
            let query_id = query_mphf.hash(&query_name) as usize;
            if !seen_queries[query_id] {
                seen_queries[query_id] = true;
                let query = &mut queries[query_id];
                query.name = query_name;
                query.length = paf_query_length(l);
            }
            let target_name: String = paf_target(l);
            let target_id = target_mphf.hash(&target_name) as usize;
            if !seen_targets[target_id] {
                seen_targets[target_id] = true;
                let target = &mut targets[target_id];
                target.name = target_name;
                target.length = paf_target_length(l);
            }
        });
        let mut targets_sort = targets.clone();
        targets_sort.sort_by(|a, b| b.length.partial_cmp(&a.length).unwrap());
        let mut target_idx: usize = 0;
        let mut target_offset: usize = 0;
        targets_sort.iter().for_each(|t| {
            let target_id = target_mphf.hash(&t.name) as usize;
            let target = &mut targets[target_id];
            target.rank = target_idx;
            target_idx += 1;
            target.offset = target_offset;
            target_offset += target.length;
        });
        let mut queries_sort = queries.clone();
        queries_sort.sort_by(|a, b| b.length.partial_cmp(&a.length).unwrap());
        let mut query_idx: usize = 0;
        let mut query_offset: usize = 0;
        queries_sort.iter().for_each(|q| {
            let query_id = query_mphf.hash(&q.name) as usize;
            let query = &mut queries[query_id];
            query.rank = query_idx;
            query_idx += 1;
            query.offset = query_offset;
            query_offset += query.length;
        });
        PafFile {
            filename: filename.to_string(),
            targets,
            queries,
            target_mphf,
            query_mphf,
            target_length: target_offset as f64,
            query_length: query_offset as f64,
            identity_tag: None,
        }
    }
    fn query_range(self: &PafFile, name: &str, start: usize, end: usize) -> (usize, usize) {
        //println!("query_range {} {}", start, end);
        let query_id = self.query_mphf.hash(&name.into()) as usize;
        //let length = self.query_length(query_id);
        let gstart = self.global_query_start(query_id);
        //println!("global query start {}", gstart);
        //println!("query length {}", length);
        let final_start = gstart + start;
        let final_end = gstart + end;
        (final_start, final_end)
    }
    fn target_range(self: &PafFile, name: &str, start: usize, end: usize) -> (usize, usize) {
        //println!("target_range {} {}", start, end);
        let target_id = self.target_mphf.hash(&name.into()) as usize;
        //let length = self.target_length(target_id);
        let gstart = self.global_target_start(target_id);
        //println!("global target start {}", gstart);
        let final_start = gstart + start;
        let final_end = gstart + end;
        (final_start, final_end)
    }
    fn global_query_start(self: &PafFile, idx: usize) -> usize {
        self.queries[idx].offset
    }
    #[allow(dead_code)]
    fn query_length(self: &PafFile, idx: usize) -> usize {
        self.queries[idx].length
    }
    fn global_target_start(self: &PafFile, idx: usize) -> usize {
        self.targets[idx].offset
    }
    #[allow(dead_code)]
    fn target_length(self: &PafFile, idx: usize) -> usize {
        self.targets[idx].length
    }
    fn global_start(self: &PafFile, line: &str, query_rev: bool) -> (usize, usize) {
        let query_id = self.query_mphf.hash(&paf_query(line)) as usize;
        let target_id = self.target_mphf.hash(&paf_target(line)) as usize;
        (
            self.global_target_start(target_id) + paf_target_begin(line),
            if query_rev {
                self.global_query_start(query_id) + paf_query_end(line)
            } else {
                self.global_query_start(query_id) + paf_query_begin(line)
            },
        )
    }
    fn for_each_match<F>(self: &PafFile, line: &str, mut func: F)
    where
        F: FnMut(char, f64, usize, bool, usize, usize),
    {
        let query_rev = paf_query_is_rev(line);
        let identity = match self.identity_tag.as_ref() {
            Some(tag) => paf_line_identity(line, tag.as_str()),
            None => 1.0,
        };
        let (x, y) = self.global_start(line, query_rev);
        let mut target_pos = x;
        let mut query_pos = y;
        // find and walk the cigar string
        //println!("{}", line);
        let mut cigars = line
            .split('\t')
            .filter(|s| s.starts_with("cg:Z:"))
            .map(|s| s.strip_prefix("cg:Z:").unwrap())
            .collect::<Vec<&str>>();

        // Compute the fake CIGAR, in case the real one is missing in the line
        let query_len = paf_query_end(line) - paf_query_begin(line);
        let target_len = paf_target_end(line) - paf_target_begin(line);
        let fake_cigar = format!(
            "{}M",
            if query_len < target_len {
                query_len
            } else {
                target_len
            }
        );
        if cigars.is_empty() {
            cigars.push(fake_cigar.as_str());
        }

        for cigar in cigars {
            //println!("{}", cigar);
            let mut first: usize = 0;
            for (i, b) in cigar.bytes().enumerate() {
                let c = b as char;
                //println!("{} {}", i, b as char);
                match c {
                    'M' | '=' | 'X' => {
                        let n = cigar[first..i].parse::<usize>().unwrap();
                        func(c, identity, target_pos, query_rev, query_pos, n);
                        query_pos = if query_rev {
                            query_pos.saturating_sub(n)
                        } else {
                            query_pos + n
                        };
                        target_pos += n;
                        first = i + 1;
                    }
                    'D' => {
                        let n = cigar[first..i].parse::<usize>().unwrap();
                        target_pos += n;
                        first = i + 1;
                    }
                    'I' => {
                        let n = cigar[first..i].parse::<usize>().unwrap();
                        query_pos = if query_rev {
                            query_pos.saturating_sub(n)
                        } else {
                            query_pos + n
                        };
                        first = i + 1;
                    }
                    _ => {}
                }
            }
        }
    }
    fn for_each_match_in_file<F>(self: &PafFile, mut func: F)
    where
        F: FnMut(char, f64, usize, bool, usize, usize),
    {
        for_each_line_in_file(&self.filename, |line: &str| {
            /*
            let (x, y) = self.global_start(line);
            println!(
                "{} {} {} {} {} {}",
                paf_query(line),
                paf_query_begin(line),
                paf_target(line),
                paf_target_begin(line),
                x,
                y
            );
             */
            self.for_each_match(line, &mut func);
            //println!();
        });
    }
    fn get_axes(self: &PafFile, major_axis: usize) -> (usize, usize) {
        //let max_length = cmp::max(self.query_length, self.target_length) as f64;
        //println!("query = {} target = {}", self.query_length, self.target_length);
        if self.target_length > self.query_length {
            let ratio = self.query_length / self.target_length;
            (major_axis, (major_axis as f64 * ratio) as usize)
        } else {
            let ratio = self.target_length / self.query_length;
            ((major_axis as f64 * ratio) as usize, major_axis)
        }
    }
    fn get_axes_zoom(
        self: &PafFile,
        major_axis: usize,
        zoom: ((usize, usize), (usize, usize)),
    ) -> (usize, usize) {
        let t_length = zoom.0 .1 - zoom.0 .0;
        let q_length = zoom.1 .1 - zoom.1 .0;
        if t_length > q_length {
            let ratio = q_length as f64 / t_length as f64;
            (major_axis, (major_axis as f64 * ratio) as usize)
        } else {
            let ratio = t_length as f64 / q_length as f64;
            ((major_axis as f64 * ratio) as usize, major_axis)
        }
    }
    fn project_xy(self: &PafFile, x: usize, y: usize, axes: (usize, usize)) -> (f64, f64) {
        //println!("axes {} {}", axes.0, axes.1);
        (
            axes.0 as f64 * (x as f64 / self.target_length),
            axes.1 as f64 * (y as f64 / self.query_length),
        )
    }
    fn project_xy_zoom(
        self: &PafFile,
        x: usize,
        y: usize,
        axes: (usize, usize),
        zoom: ((usize, usize), (usize, usize)),
    ) -> (f64, f64) {
        //println!("axes {} {}", axes.0, axes.1);
        (
            axes.0 as f64 * (((x as f64) - (zoom.0 .0 as f64)) / ((zoom.0 .1 - zoom.0 .0) as f64)),
            axes.1 as f64 * (((y as f64) - (zoom.1 .0 as f64)) / ((zoom.1 .1 - zoom.1 .0) as f64)),
        )
    }
    /*
        fn get_pixel(self: &PafFile, x: i64, y: i64, axes: (usize, usize)) -> usize {
            let (q, t) = axes;
            x * (self.query_length / q as f64).round() as usize * q + y * (self.target_length / t as f64).round() as usize
        }
    */
}

fn parse_bed_file(filename: &str) -> Vec<BedRegion> {
    let mut regions = Vec::new();
    for_each_line_in_file(filename, |line| {
        if line.starts_with('#') || line.is_empty() {
            return;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 3 {
            let chr = parts[0].to_string();
            let start = parts[1].parse::<usize>().unwrap_or(0);
            let end = parts[2].parse::<usize>().unwrap_or(0);
            let name = if parts.len() > 3 && !parts[3].is_empty() {
                Some(parts[3].to_string())
            } else {
                None
            };
            let color = if parts.len() > 4 && !parts[4].is_empty() {
                Some(parts[4].to_string())
            } else {
                None
            };
            regions.push(BedRegion {
                chr,
                start,
                end,
                name,
                color,
            });
        }
    });
    regions
}

fn parse_bedpe_file(filename: &str) -> Vec<BedpeRegion> {
    let mut regions = Vec::new();
    for_each_line_in_file(filename, |line| {
        if line.starts_with('#') || line.is_empty() {
            return;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 6 {
            let chr1 = parts[0].to_string();
            let start1 = parts[1].parse::<usize>().unwrap_or(0);
            let end1 = parts[2].parse::<usize>().unwrap_or(0);
            let chr2 = parts[3].to_string();
            let start2 = parts[4].parse::<usize>().unwrap_or(0);
            let end2 = parts[5].parse::<usize>().unwrap_or(0);
            let name = if parts.len() > 6 && !parts[6].is_empty() {
                Some(parts[6].to_string())
            } else {
                None
            };
            let color = if parts.len() > 7 && !parts[7].is_empty() {
                Some(parts[7].to_string())
            } else {
                None
            };
            regions.push(BedpeRegion {
                chr1,
                start1,
                end1,
                chr2,
                start2,
                end2,
                name,
                color,
            });
        }
    });
    regions
}

fn main() {
    let matches = App::new("pafplot")
        .version("0.1.0")
        .author("Erik Garrison <erik.garrison@gmail.com>")
        .about("Generate dot plots from pairwise DNA alignments in PAF format (PNG; optional interactive HTML and labeled PDF with mouseover tooltips)")
        .arg(
            Arg::with_name("INPUT")
                .required(true)
                .takes_value(true)
                .index(1)
                .help("input PAF file"),
        )
        .arg(
            Arg::with_name("html")
                .takes_value(false)
                .long("html")
                .help("Generate an interactive HTML viewer in addition to the default PNG output."),
        )
        .arg(
            Arg::with_name("pdf")
                .takes_value(false)
                .long("pdf")
                .help("Emit a PDF embedding the dot plot, with row/column labels and mouseover tooltips."),
        )
        .arg(
            Arg::with_name("dark")
                .takes_value(false)
                .short("d")
                .long("dark")
                .help("Render using dark theme (white on black)."),
        )
        .arg(
            Arg::with_name("size")
                .takes_value(true)
                .short("s")
                .long("size")
                .help("The major axis of the plot, in pixels. [default: 1000]"),
        )
        .arg(
            Arg::with_name("range")
                .takes_value(true)
                .short("r")
                .long("range")
                .help("Plot the given 2D range in target and query rather than the full matrix seqA:10-200,seqB:300-400"),
        )
        .arg(
            Arg::with_name("html")
                .takes_value(false)
                .short("h")
                .long("html")
                .help("Generate an interactive HTML viewer (default output format)."),
        )
        .arg(
            Arg::with_name("color-by")
                .takes_value(true)
                .long("color-by")
                .help("Color alignment lines by this numeric PAF tag (e.g. bi, gi, dv, de). \
                       Lines missing the tag fall back to BLAST identity (matches/alen)."),
        )
        .arg(
            Arg::with_name("color-by-identity")
                .takes_value(false)
                .long("color-by-identity")
                .help("Color alignment lines by identity (equivalent to --color-by bi). \
                       Uses the bi: tag when present, else matches/alen from the fixed fields."),
        )
        .arg(
            Arg::with_name("color-scale")
                .takes_value(true)
                .long("color-scale")
                .help("Identity range mapping onto the line colormap, min:max in [0,1]. [default: 0.8:1.0]"),
        )
        .arg(
            Arg::with_name("invert-color-scale")
                .takes_value(false)
                .long("invert-color-scale")
                .help("Flip the color scale (for divergence-style tags such as dv or de)."),
        )
        .arg(
            Arg::with_name("output")
                .takes_value(true)
                .short("o")
                .long("output")
                .help("Output filename (defaults to input.paf.png)."),
        )
        .arg(
            Arg::with_name("bed")
                .takes_value(true)
                .multiple(true)
                .long("bed")
                .help("BED file(s) for marking regions as semi-transparent overlays"),
        )
        .arg(
            Arg::with_name("bedpe")
                .takes_value(true)
                .multiple(true)
                .long("bedpe")
                .help("BEDPE file(s) for marking 2D ranges as semi-transparent overlays"),
        )
        .arg(
            Arg::with_name("highlight-fill")
                .takes_value(false)
                .long("highlight-fill")
                .help("Fill highlight boxes with a translucent tint in addition to the outline \
                       (default draws outline-only boxes so the identity coloring stays visible)."),
        )
        .arg(
            Arg::with_name("pdf-label-size")
                .takes_value(true)
                .long("pdf-label-size")
                .help("PDF label font size in points [default: 9]. Increase when the figure will be placed smaller than its page width."),
        )
        .arg(
            Arg::with_name("highlight")
                .takes_value(true)
                .multiple(true)
                .long("highlight")
                .help("Highlight file(s): TSV of 2D regions drawn on the PDF as labeled boxes \
                       (qname\tqstart\tqend\ttname\ttstart\ttend\t[#RRGGBB]\t[label])"),
        )
        .get_matches();

    let filename = matches.value_of("INPUT").unwrap();
    let mut paf = PafFile::new(filename);

    // line coloring by identity (or any numeric tag)
    let color_tag: Option<String> = if matches.is_present("color-by") {
        Some(matches.value_of("color-by").unwrap().to_string())
    } else if matches.is_present("color-by-identity") {
        Some("bi".to_string())
    } else {
        None
    };
    paf.identity_tag = color_tag.clone();
    let (id_min, id_max) = match matches.value_of("color-scale") {
        Some(spec) => {
            let parts: Vec<&str> = spec.split(':').collect();
            if parts.len() != 2 {
                panic!("[pafplot] invalid --color-scale {}, expected min:max", spec);
            }
            (
                parts[0].parse::<f64>().unwrap(),
                parts[1].parse::<f64>().unwrap(),
            )
        }
        None => (0.8, 1.0),
    };
    let invert_scale = matches.is_present("invert-color-scale");

    let major_axis = matches
        .value_of("size")
        .unwrap_or("1000")
        .parse::<usize>()
        .unwrap();

    let default_output = format!("{filename}.png");

    let output_filename = matches.value_of("output").unwrap_or(&default_output);

    let dark = matches.is_present("dark");

    // Parse BED and BEDPE files
    let bed_regions: Vec<BedRegion> = matches
        .values_of("bed")
        .map(|files| files.flat_map(parse_bed_file).collect())
        .unwrap_or_default();

    let bedpe_regions: Vec<BedpeRegion> = matches
        .values_of("bedpe")
        .map(|files| files.flat_map(parse_bedpe_file).collect())
        .unwrap_or_default();

    let highlight_regions: Vec<HighlightRegion> = matches
        .values_of("highlight")
        .map(|files| files.flat_map(parse_highlight_file).collect())
        .unwrap_or_default();

    let using_zoom = matches.is_present("range");
    let (target_range, query_range): ((usize, usize), (usize, usize)) = if using_zoom {
        let splitv = matches
            .value_of("range")
            .unwrap()
            .split(',')
            .collect::<Vec<&str>>();
        if splitv.len() != 2 {
            panic!(
                "[pafplot::main] invalid range specification {}",
                matches.value_of("range").unwrap()
            );
        }

        let target_name = splitv[0].split(':').next().unwrap();
        let target_start = splitv[0]
            .split(':')
            .nth(1)
            .unwrap()
            .split('-')
            .next()
            .unwrap()
            .parse::<usize>()
            .unwrap();
        let target_end = splitv[0]
            .split(':')
            .nth(1)
            .unwrap()
            .split('-')
            .nth(1)
            .unwrap()
            .parse::<usize>()
            .unwrap();

        let query_name = splitv[1].split(':').next().unwrap();
        let query_start = splitv[1]
            .split(':')
            .nth(1)
            .unwrap()
            .split('-')
            .next()
            .unwrap()
            .parse::<usize>()
            .unwrap();
        let query_end = splitv[1]
            .split(':')
            .nth(1)
            .unwrap()
            .split('-')
            .nth(1)
            .unwrap()
            .parse::<usize>()
            .unwrap();

        (
            paf.target_range(target_name, target_start, target_end),
            paf.query_range(query_name, query_start, query_end),
        )
    } else {
        (
            (0, paf.target_length as usize),
            (0, paf.query_length as usize),
        )
    };

    // colors we use
    let white = RGB8 {
        r: 255_u8,
        g: 255,
        b: 255,
    };
    let black = white.map(|ch| 255 - ch);

    /*
    println!(
        "getting axes query=({} {}) target=({} {})",
        query_range.0, query_range.1, target_range.0, target_range.1
    );
     */
    let axes = if using_zoom {
        paf.get_axes_zoom(major_axis, (target_range, query_range))
    } else {
        paf.get_axes(major_axis)
    };
    //println!("axes = {} {}, target_length = {}, query_length = {}", axes.0, axes.1, paf.target_length, paf.query_length);
    let mut raw = vec![0u8; axes.0 * axes.1 * 3];
    let pixels = raw.as_rgb_mut();

    if dark {
        for i in pixels.iter_mut() {
            *i = black;
        }
    } else {
        for i in pixels.iter_mut() {
            *i = white;
        }
    }
    let get_color = if dark {
        |val: f64| RGB8 {
            r: ((255.0 * val).round() as u8),
            g: ((255.0 * val).round() as u8),
            b: ((255.0 * val).round() as u8),
        }
    } else {
        |val: f64| RGB8 {
            r: ((255.0 * (1.0 - val)).round() as u8),
            g: ((255.0 * (1.0 - val)).round() as u8),
            b: ((255.0 * (1.0 - val)).round() as u8),
        }
    };

    let get_coords = |x: usize, y: usize| {
        if using_zoom {
            paf.project_xy_zoom(x, y, axes, (target_range, query_range))
        } else {
            paf.project_xy(x, y, axes)
        }
    };

    // draw border and grid
    // Draw border first
    let border_width = 0.5;
    // Top border
    let start = get_coords(0, 0);
    let end = get_coords(paf.target_length as usize, 0);
    for ((i, j), val) in XiaolinWu::<f64, i64>::new(start, end) {
        if i >= 0 && i < (axes.0 as i64) && j >= 0 && j < (axes.1 as i64) {
            let i: usize = (i as usize) + (((axes.1 - 1) - j as usize) * axes.0);
            pixels[i] = get_color(val * border_width);
        }
    }
    // Right border
    let start = get_coords(paf.target_length as usize, 0);
    let end = get_coords(paf.target_length as usize, paf.query_length as usize);
    for ((i, j), val) in XiaolinWu::<f64, i64>::new(start, end) {
        if i >= 0 && i < (axes.0 as i64) && j >= 0 && j < (axes.1 as i64) {
            let i: usize = (i as usize) + (((axes.1 - 1) - j as usize) * axes.0);
            pixels[i] = get_color(val * border_width);
        }
    }
    // Bottom border
    let start = get_coords(0, paf.query_length as usize);
    let end = get_coords(paf.target_length as usize, paf.query_length as usize);
    for ((i, j), val) in XiaolinWu::<f64, i64>::new(start, end) {
        if i >= 0 && i < (axes.0 as i64) && j >= 0 && j < (axes.1 as i64) {
            let i: usize = (i as usize) + (((axes.1 - 1) - j as usize) * axes.0);
            pixels[i] = get_color(val * border_width);
        }
    }
    // Left border
    let start = get_coords(0, 0);
    let end = get_coords(0, paf.query_length as usize);
    for ((i, j), val) in XiaolinWu::<f64, i64>::new(start, end) {
        if i >= 0 && i < (axes.0 as i64) && j >= 0 && j < (axes.1 as i64) {
            let i: usize = (i as usize) + (((axes.1 - 1) - j as usize) * axes.0);
            pixels[i] = get_color(val * border_width);
        }
    }

    // draw our grid
    paf.targets.iter().for_each(|target| {
        if target.offset > 0 {
            let start = get_coords(target.offset, 0);
            let end = get_coords(target.offset, paf.query_length as usize);
            for ((i, j), val) in XiaolinWu::<f64, i64>::new(start, end) {
                if i >= 0 && i < (axes.0 as i64) && j >= 0 && j < (axes.1 as i64) {
                    let i: usize = (i as usize) + (((axes.1 - 1) - j as usize) * axes.0);
                    pixels[i] = get_color(val * 0.2);
                }
            }
        }
    });
    paf.queries.iter().for_each(|query| {
        if query.offset > 0 {
            let start = get_coords(0, query.offset);
            let end = get_coords(paf.target_length as usize, query.offset);
            for ((i, j), val) in XiaolinWu::<f64, i64>::new(start, end) {
                if i >= 0 && i < (axes.0 as i64) && j >= 0 && j < (axes.1 as i64) {
                    let i: usize = (i as usize) + (((axes.1 - 1) - j as usize) * axes.0);
                    pixels[i] = get_color(val * 0.2);
                }
            }
        }
    });

    // for each match, we draw a line on our raster using Xiaolin Wu's antialiased line algorithm
    // lines are colored by identity (turbo colormap) when --color-by/--color-by-identity is set
    let fg_line = if dark { white } else { black };
    let blend = |bg: RGB8, lc: RGB8, a: f64| -> RGB8 {
        RGB8 {
            r: (bg.r as f64 + (lc.r as f64 - bg.r as f64) * a).round() as u8,
            g: (bg.g as f64 + (lc.g as f64 - bg.g as f64) * a).round() as u8,
            b: (bg.b as f64 + (lc.b as f64 - bg.b as f64) * a).round() as u8,
        }
    };
    let line_color_for = |identity: f64| -> RGB8 {
        match color_tag {
            Some(_) => {
                let mut t = if id_max > id_min {
                    (identity - id_min) / (id_max - id_min)
                } else {
                    0.0
                };
                if invert_scale {
                    t = 1.0 - t;
                }
                turbo(t)
            }
            None => fg_line,
        }
    };
    let draw_match = |_c, identity: f64, x: usize, rev: bool, y: usize, len: usize| {
        //println!("draw_match {} {} {} {} {}", _c, x, rev, y, len);
        let lc = line_color_for(identity);
        let start = get_coords(x, y);
        let end = get_coords(x + len, if rev { y.saturating_sub(len) } else { y + len });

        /*
        println!(
            "start and end ({} {}) ({} {})",
            start.0, start.1, end.0, end.1
        );
         */
        for ((i, j), val) in XiaolinWu::<f64, i64>::new(start, end) {
            //println!("checking pixel {} {} {}", i, j, val);
            if i >= 0 && i < (axes.0 as i64) && j >= 0 && j < (axes.1 as i64) {
                //println!("drawing pixel {} {} {}", i, j, val);
                let i: usize = (i as usize) + (((axes.1 - 1) - j as usize) * axes.0);
                pixels[i] = if color_tag.is_some() {
                    blend(if dark { black } else { white }, lc, val)
                } else {
                    get_color(val)
                };
            }
        }
    };
    paf.for_each_match_in_file(draw_match);

    // per-(query, target) pair alignment statistics, used for PDF cell tooltips
    let mut pair_stats: PairStats = std::collections::HashMap::new();
    for_each_line_in_file(filename, |l: &str| {
        let q = paf_query(l);
        let t = paf_target(l);
        let matches: f64 = l.split('\t').nth(9).unwrap_or("0").parse().unwrap_or(0.0);
        let aln_len: f64 = l.split('\t').nth(10).unwrap_or("1").parse().unwrap_or(1.0);
        let entry = pair_stats.entry((q, t)).or_insert((0usize, 0u64, 0f64));
        entry.0 += 1;
        entry.1 += (paf_query_end(l) - paf_query_begin(l)) as u64;
        if aln_len > 0.0 {
            let id = matches / aln_len;
            if id > entry.2 {
                entry.2 = id;
            }
        }
    });

    // Generate PNG (default output)
    let png_filename = if output_filename.ends_with(".png") {
        output_filename.to_string()
    } else if output_filename.ends_with(".html") {
        output_filename.replace(".html", ".png")
    } else if output_filename.ends_with(".pdf") {
        output_filename.replace(".pdf", ".png")
    } else {
        format!("{output_filename}.png")
    };

    let path = &Path::new(&png_filename);
    // encode_file takes the path to the image, a u8 array,
    // the width, the height, the color mode, and the bit depth
    if let Err(e) = lodepng::encode_file(path, &raw, axes.0, axes.1, lodepng::ColorType::RGB, 8) {
        panic!("failed to write png: {:?}", e);
    }
    eprintln!("[pafplot] PNG image written to: {}", png_filename);

    // Generate PDF viewer if requested
    if matches.is_present("pdf") || output_filename.ends_with(".pdf") {
        let pdf_filename = if output_filename.ends_with(".pdf") {
            output_filename.to_string()
        } else if output_filename.ends_with(".png") {
            output_filename.replace(".png", ".pdf")
        } else if output_filename.ends_with(".html") {
            output_filename.replace(".html", ".pdf")
        } else {
            format!("{output_filename}.pdf")
        };
        generate_pdf_viewer(PdfViewerConfig {
            paf: &paf,
            axes,
            raw: &raw,
            output_filename: &pdf_filename,
            using_zoom,
            ranges: (target_range, query_range),
            pair_stats: &pair_stats,
            highlights: &highlight_regions,
            legend: if color_tag.is_some() {
                Some((id_min, id_max, invert_scale))
            } else {
                None
            },
            highlight_fill: matches.is_present("highlight-fill"),
            label_size: matches
                .value_of("pdf-label-size")
                .and_then(|v| v.parse::<f64>().ok())
                .filter(|v| *v > 0.0)
                .unwrap_or(9.0),
        });
    }

    // Generate HTML viewer if requested
    if matches.is_present("html") {
        let html_filename = if output_filename.ends_with(".html") {
            output_filename.to_string()
        } else if output_filename.ends_with(".png") {
            output_filename.replace(".png", ".html")
        } else {
            format!("{output_filename}.html")
        };
        generate_html_viewer(HtmlViewerConfig {
            paf: &paf,
            _axes: axes,
            output_filename: &html_filename,
            dark,
            using_zoom,
            ranges: (target_range, query_range),
            bed_regions: &bed_regions,
            bedpe_regions: &bedpe_regions,
        });
    }
}

fn collect_alignment_data(paf: &PafFile) -> (String, String, String) {
    let mut alignments = Vec::new();
    let mut detailed_alignments = Vec::new();
    let mut summary_grid = std::collections::HashMap::new();
    let grid_size = 1000; // Grid resolution for summary view

    // Collect both simple alignments and detailed CIGAR data
    for_each_line_in_file(&paf.filename, |line| {
        let query_rev = paf_query_is_rev(line);
        let (x, y) = paf.global_start(line, query_rev);

        // Store the full alignment for line drawing
        let query_len = paf_query_end(line) - paf_query_begin(line);
        let target_len = paf_target_end(line) - paf_target_begin(line);

        alignments.push(format!(
            r#"{{"x":{x},"y":{y},"queryLen":{query_len},"targetLen":{target_len},"rev":{query_rev}}}"#
        ));

        // Update summary grid
        let grid_x = (x as f64 / paf.target_length * grid_size as f64) as usize;
        let grid_y = (y as f64 / paf.query_length * grid_size as f64) as usize;
        let entry = summary_grid.entry((grid_x, grid_y)).or_insert((0, 0, 0));
        entry.0 += 1;
        entry.1 += query_len.max(target_len); // Use the larger segment length for summary
        if query_rev {
            entry.2 += 1;
        }

        // Collect detailed CIGAR operations for base-level view
        let mut cigar_ops = Vec::new();

        paf.for_each_match(line, |c, _identity, target_pos, rev, query_pos, len| {
            let op_type = match c {
                'M' | '=' => "match",
                'X' => "mismatch",
                'I' => "insertion",
                'D' => "deletion",
                _ => "unknown",
            };

            cigar_ops.push(format!(
                r#"{{"type":"{op_type}","x":{target_pos},"y":{query_pos},"len":{len},"rev":{rev}}}"#
            ));
        });

        if !cigar_ops.is_empty() {
            detailed_alignments.push(format!(
                r#"{{"x":{},"y":{},"queryLen":{},"targetLen":{},"rev":{},"ops":[{}]}}"#,
                x,
                y,
                query_len,
                target_len,
                query_rev,
                cigar_ops.join(",")
            ));
        }
    });

    // Convert summary grid to JSON
    let mut summary_data = Vec::new();
    for ((gx, gy), (count, total_len, rev_count)) in summary_grid {
        let density = count as f64;
        let avg_len = total_len as f64 / count as f64;
        let rev_ratio = rev_count as f64 / count as f64;

        summary_data.push(format!(
            r#"{{"gx":{gx},"gy":{gy},"density":{density},"avgLen":{avg_len},"revRatio":{rev_ratio}}}"#
        ));
    }

    let alignments_json = format!("[{}]", alignments.join(","));
    let summary_json = format!("[{}]", summary_data.join(","));
    let detailed_json = format!("[{}]", detailed_alignments.join(","));

    (alignments_json, summary_json, detailed_json)
}

struct HtmlViewerConfig<'a> {
    paf: &'a PafFile,
    _axes: (usize, usize),
    output_filename: &'a str,
    dark: bool,
    using_zoom: bool,
    ranges: ((usize, usize), (usize, usize)),
    bed_regions: &'a [BedRegion],
    bedpe_regions: &'a [BedpeRegion],
}

fn generate_html_viewer(config: HtmlViewerConfig) {
    // Collect alignment data for canvas rendering
    let (alignments_json, summary_json, detailed_json) = collect_alignment_data(config.paf);

    // Generate sequence metadata as JSON
    let mut targets_json = String::from("[");
    for (i, target) in config.paf.targets.iter().enumerate() {
        if i > 0 {
            targets_json.push(',');
        }
        targets_json.push_str(&format!(
            r#"{{"name":"{}","length":{},"rank":{},"offset":{}}}"#,
            target.name.replace('"', r#"\""#),
            target.length,
            target.rank,
            target.offset
        ));
    }
    targets_json.push(']');

    let mut queries_json = String::from("[");
    for (i, query) in config.paf.queries.iter().enumerate() {
        if i > 0 {
            queries_json.push(',');
        }
        queries_json.push_str(&format!(
            r#"{{"name":"{}","length":{},"rank":{},"offset":{}}}"#,
            query.name.replace('"', r#"\""#),
            query.length,
            query.rank,
            query.offset
        ));
    }
    queries_json.push(']');

    // Convert BED regions to JSON
    let mut bed_json = String::from("[");
    for (i, region) in config.bed_regions.iter().enumerate() {
        if i > 0 {
            bed_json.push(',');
        }
        bed_json.push_str(&format!(
            r#"{{"chr":"{}","start":{},"end":{},"name":{},"color":{}}}"#,
            region.chr,
            region.start,
            region.end,
            region
                .name
                .as_ref()
                .map(|n| format!(r#""{n}""#))
                .unwrap_or_else(|| "null".to_string()),
            region
                .color
                .as_ref()
                .map(|c| format!(r#""{c}""#))
                .unwrap_or_else(|| "null".to_string())
        ));
    }
    bed_json.push(']');

    // Convert BEDPE regions to JSON
    let mut bedpe_json = String::from("[");
    for (i, region) in config.bedpe_regions.iter().enumerate() {
        if i > 0 {
            bedpe_json.push(',');
        }
        bedpe_json.push_str(&format!(
            r#"{{"chr1":"{}","start1":{},"end1":{},"chr2":"{}","start2":{},"end2":{},"name":{},"color":{}}}"#,
            region.chr1,
            region.start1,
            region.end1,
            region.chr2,
            region.start2,
            region.end2,
            region.name.as_ref().map(|n| format!(r#""{n}""#)).unwrap_or_else(|| "null".to_string()),
            region.color.as_ref().map(|c| format!(r#""{c}""#)).unwrap_or_else(|| "null".to_string())
        ));
    }
    bedpe_json.push(']');

    #[allow(clippy::format_in_format_args)]
    let html_content = format!(
        r#"<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PAF Plot Viewer - {}</title>
    <style>
        body {{
            font-family: 'Monaco', 'Menlo', 'Ubuntu Mono', monospace;
            background-color: {};
            color: {};
            margin: 0;
            padding: 0;
            overflow: hidden;
            width: 100vw;
            height: 100vh;
        }}
        #plotCanvas {{
            position: absolute;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            cursor: crosshair;
        }}
        .tooltip {{
            position: absolute;
            background-color: {};
            border: 1px solid {};
            padding: 8px;
            border-radius: 4px;
            font-size: 11px;
            pointer-events: none;
            z-index: 1000;
            white-space: nowrap;
            display: none;
        }}
        .labels {{
            position: fixed;
            font-size: 10px;
            font-family: 'Monaco', 'Menlo', 'Ubuntu Mono', monospace;
            pointer-events: none;
            z-index: 15;
        }}
        .target-labels {{
            bottom: 10px;
            left: 0;
            width: 100%;
            height: 100px;
        }}
        .query-labels {{
            right: 10px;
            top: 0;
            width: 150px;
            height: 100%;
        }}
        .label {{
            position: absolute;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
            background-color: {};
            padding: 2px 4px;
            border-radius: 2px;
            font-size: 9px;
            opacity: 0.9;
        }}
        .target-label {{
            transform-origin: left top;
            transform: rotate(-45deg);
            max-width: 100px;
        }}
        .query-label {{
            text-align: right;
            right: 0;
            max-width: 140px;
        }}
        .info-panel {{
            position: fixed;
            top: 20px;
            right: 20px;
            background-color: {};
            border: 1px solid {};
            padding: 15px;
            border-radius: 5px;
            min-width: 250px;
            font-size: 11px;
            z-index: 20;
            opacity: 0.9;
        }}
        .config-panel {{
            position: fixed;
            top: 20px;
            left: 20px;
            background-color: {};
            border: 1px solid {};
            padding: 15px;
            border-radius: 5px;
            min-width: 200px;
            font-size: 11px;
            z-index: 20;
            opacity: 0.9;
        }}
        .config-panel h3 {{
            margin-top: 0;
            margin-bottom: 10px;
        }}
        .config-control {{
            margin-bottom: 10px;
        }}
        .config-control label {{
            display: inline-block;
            width: 100px;
            font-size: 10px;
        }}
        .config-control input[type="range"] {{
            width: 100px;
        }}
        .config-control span {{
            font-size: 10px;
            margin-left: 5px;
        }}
        #renderingMode {{
            font-weight: bold;
            color: {};
        }}
        .minimap {{
            position: fixed;
            bottom: 20px;
            right: 20px;
            width: 200px;
            height: 150px;
            background-color: {};
            border: 2px solid {};
            z-index: 25;
            overflow: hidden;
        }}
        .minimap-canvas {{
            width: 100%;
            height: 100%;
        }}
        .minimap-viewport {{
            position: absolute;
            border: 2px solid {};
            background-color: {};
            pointer-events: none;
        }}
        .controls {{
            margin-bottom: 20px;
        }}
        .zoom-info {{
            font-size: 11px;
            color: {};
            margin-top: 10px;
        }}
    </style>
</head>
<body>
    <canvas id="plotCanvas"></canvas>
    <div class="tooltip" id="tooltip"></div>
    <div class="labels target-labels" id="targetLabels"></div>
    <div class="labels query-labels" id="queryLabels"></div>
    
    <div class="config-panel">
        <h3>Settings</h3>
        <div class="config-control">
            <label>Line Width:</label>
            <input type="range" id="lineWidthSlider" min="0.1" max="2" step="0.1" value="0.5">
            <span id="lineWidthValue">0.5x</span>
        </div>
        <div class="config-control">
            <label>Grid Opacity:</label>
            <input type="range" id="gridOpacitySlider" min="0" max="1" step="0.1" value="0.3">
            <span id="gridOpacityValue">30%</span>
        </div>
        <div class="config-control">
            <label>Show Details:</label>
            <input type="checkbox" id="showDetailsCheckbox" checked>
        </div>
        <div class="config-control">
            <label>Show Labels:</label>
            <input type="checkbox" id="showLabelsCheckbox" checked>
        </div>
        <div class="config-control">
            <button onclick="resetZoom()">Reset View</button>
            <button id="undoZoomBtn" onclick="undoZoom()" disabled>Undo Zoom</button>
        </div>
        <div class="config-control">
            <label>Rendering Mode:</label>
            <span id="renderingMode">Standard</span>
        </div>
    </div>
        
        <div class="info-panel">
            <h3>Current Position</h3>
            <div id="position">Move mouse over plot</div>
            <h3>Sequence Info</h3>
            <div id="sequenceInfo">Hover over plot for details</div>
            <div class="zoom-info">
                Zoom: <span id="zoomLevel">1.0x</span><br>
                Pan: <span id="panInfo">0, 0</span><br>
                Viewport: <span id="viewportInfo">0 bp</span>
            </div>
        </div>
        
        <div class="minimap" id="minimap">
            <canvas class="minimap-canvas" id="minimapCanvas"></canvas>
            <div class="minimap-viewport" id="minimapViewport"></div>
        </div>
    </div>

    <script>
        const canvas = document.getElementById('plotCanvas');
        const gl = canvas.getContext('webgl');
        const tooltip = document.getElementById('tooltip');
        const position = document.getElementById('position');
        const sequenceInfo = document.getElementById('sequenceInfo');
        const zoomLevel = document.getElementById('zoomLevel');
        const panInfo = document.getElementById('panInfo');
        const minimapCanvas = document.getElementById('minimapCanvas');
        const minimapCtx = minimapCanvas.getContext('2d');
        const minimapViewport = document.getElementById('minimapViewport');
        
        // Configuration controls
        const lineWidthSlider = document.getElementById('lineWidthSlider');
        const lineWidthValue = document.getElementById('lineWidthValue');
        const gridOpacitySlider = document.getElementById('gridOpacitySlider');
        const gridOpacityValue = document.getElementById('gridOpacityValue');
        const showDetailsCheckbox = document.getElementById('showDetailsCheckbox');
        const showLabelsCheckbox = document.getElementById('showLabelsCheckbox');
        const targetLabels = document.getElementById('targetLabels');
        const queryLabels = document.getElementById('queryLabels');
        const renderingMode = document.getElementById('renderingMode');
        const viewportInfo = document.getElementById('viewportInfo');
        
        // Configuration values
        let lineWidthMultiplier = 0.5;
        let gridOpacity = 0.3;
        let showDetails = true;
        let showLabels = true;

        // Plot data
        const alignments = {};
        const summaryData = {};
        const detailedAlignments = {};
        const targets = {};
        const queries = {};
        const bedRegions = {};
        const bedpeRegions = {};
        let canvasWidth = window.innerWidth;
        let canvasHeight = window.innerHeight;
        const targetLength = {};
        const queryLength = {};
        
        // Calculate plot dimensions maintaining 1:1 aspect ratio
        let plotWidth, plotHeight, plotOffsetX, plotOffsetY;
        
        function calculatePlotDimensions() {{
            const aspectRatio = targetLength / queryLength;
            const viewportAspect = canvasWidth / canvasHeight;
            
            if (aspectRatio > viewportAspect) {{
                // Target is wider - fit to width
                plotWidth = canvasWidth * 0.9;  // 90% of viewport width
                plotHeight = plotWidth / aspectRatio;
                plotOffsetX = canvasWidth * 0.05;
                plotOffsetY = (canvasHeight - plotHeight) / 2;
            }} else {{
                // Query is taller - fit to height
                plotHeight = canvasHeight * 0.9;  // 90% of viewport height
                plotWidth = plotHeight * aspectRatio;
                plotOffsetX = (canvasWidth - plotWidth) / 2;
                plotOffsetY = canvasHeight * 0.05;
            }}
        }}
        const usingZoom = {};
        const targetRange = {};
        const queryRange = {};
        const darkMode = {};
        const GRID_SIZE = 1000;
        const BASE_VIEW_THRESHOLD = 100; // Show base view when < 100bp visible

        // Zoom and pan state
        let zoom = 1.0;
        let panX = 0;
        let panY = 0;
        let isDragging = false;
        let dragStartX = 0;
        let dragStartY = 0;
        let dragStartPanX = 0;
        let dragStartPanY = 0;
        
        // Zoom-to-box state
        let isBoxSelecting = false;
        let boxStartX = 0;
        let boxStartY = 0;
        let boxEndX = 0;
        let boxEndY = 0;
        
        // Zoom history for undo
        const zoomHistory = [];
        const maxZoomHistory = 50;
        
        function saveZoomState() {{
            zoomHistory.push({{ zoom: zoom, panX: panX, panY: panY }});
            if (zoomHistory.length > maxZoomHistory) {{
                zoomHistory.shift();
            }}
            updateUndoButton();
        }}
        
        function undoZoom() {{
            if (zoomHistory.length > 0) {{
                const state = zoomHistory.pop();
                zoom = state.zoom;
                panX = state.panX;
                panY = state.panY;
                updateDisplay();
                updateUndoButton();
            }}
        }}
        
        function updateUndoButton() {{
            const undoBtn = document.getElementById('undoZoomBtn');
            if (undoBtn) {{
                undoBtn.disabled = zoomHistory.length === 0;
            }}
        }}

        // Colors
        const backgroundColor = darkMode ? '#1a1a1a' : '#ffffff';
        const lineColor = darkMode ? '#ffffff' : '#000000';
        const borderColor = darkMode ? '#444444' : '#cccccc';
        
        // WebGL setup
        let shaderProgram;
        let positionBuffer;
        let colorBuffer;
        let lineVertices = [];
        let lineColors = [];
        
        // Vertex shader source
        const vsSource = `
            attribute vec2 a_position;
            attribute vec4 a_color;
            uniform vec2 u_resolution;
            uniform vec2 u_pan;
            uniform float u_zoom;
            varying vec4 v_color;
            
            void main() {{
                // Apply zoom and pan
                vec2 position = (a_position * u_zoom + u_pan) / u_resolution * 2.0 - 1.0;
                position.y *= -1.0; // Flip Y axis
                gl_Position = vec4(position, 0.0, 1.0);
                v_color = a_color;
            }}
        `;
        
        // Fragment shader source
        const fsSource = `
            precision mediump float;
            varying vec4 v_color;
            
            void main() {{
                gl_FragColor = v_color;
            }}
        `;
        
        // Initialize WebGL
        function initWebGL() {{
            if (!gl) {{
                alert('WebGL not supported');
                return false;
            }}
            
            // Create shaders
            const vertexShader = createShader(gl.VERTEX_SHADER, vsSource);
            const fragmentShader = createShader(gl.FRAGMENT_SHADER, fsSource);
            
            // Create program
            shaderProgram = gl.createProgram();
            gl.attachShader(shaderProgram, vertexShader);
            gl.attachShader(shaderProgram, fragmentShader);
            gl.linkProgram(shaderProgram);
            
            if (!gl.getProgramParameter(shaderProgram, gl.LINK_STATUS)) {{
                console.error('Unable to initialize shader program:', gl.getProgramInfoLog(shaderProgram));
                return false;
            }}
            
            // Create buffers
            positionBuffer = gl.createBuffer();
            colorBuffer = gl.createBuffer();
            
            // Set clear color - use a slightly different color to distinguish plot area
            const bgColor = hexToRgb(darkMode ? '#0a0a0a' : '#f8f8f8');
            gl.clearColor(bgColor.r, bgColor.g, bgColor.b, 1.0);
            
            // Enable blending for anti-aliasing
            gl.enable(gl.BLEND);
            gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
            
            // Line width will be set dynamically based on zoom
            
            return true;
        }}
        
        function createShader(type, source) {{
            const shader = gl.createShader(type);
            gl.shaderSource(shader, source);
            gl.compileShader(shader);
            
            if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {{
                console.error('Shader compilation error:', gl.getShaderInfoLog(shader));
                gl.deleteShader(shader);
                return null;
            }}
            
            return shader;
        }}
        
        function hexToRgb(hex) {{
            const result = /^#?([a-f\d]{{2}})([a-f\d]{{2}})([a-f\d]{{2}})$/i.exec(hex);
            return result ? {{
                r: parseInt(result[1], 16) / 255,
                g: parseInt(result[2], 16) / 255,
                b: parseInt(result[3], 16) / 255
            }} : {{ r: 0, g: 0, b: 0 }};
        }}

        // Coordinate projection functions - project to plot area
        function projectCoords(x, y) {{
            let projX, projY;
            if (usingZoom) {{
                // Match project_xy_zoom from Rust code
                projX = plotWidth * ((x - targetRange[0]) / (targetRange[1] - targetRange[0]));
                projY = plotHeight * ((y - queryRange[0]) / (queryRange[1] - queryRange[0]));
            }} else {{
                // Match project_xy from Rust code
                projX = plotWidth * (x / targetLength);
                projY = plotHeight * (y / queryLength);
            }}
            // Y-axis flip to match PNG (bottom-left to top-left origin)
            projY = plotHeight - projY;
            // Add plot offset
            projX += plotOffsetX;
            projY += plotOffsetY;
            return {{ x: projX, y: projY }};
        }}

        // Unified viewport calculation function
        function calculateViewportInfo() {{
            const viewLeft = -panX / zoom;
            const viewRight = (canvasWidth - panX) / zoom;
            const viewTop = -panY / zoom;
            const viewBottom = (canvasHeight - panY) / zoom;
            const viewWidth = viewRight - viewLeft;
            const viewHeight = viewBottom - viewTop;
            
            // Calculate viewport width in base pairs - use target axis for consistency
            const viewportWidthBp = usingZoom ? 
                ((viewWidth / plotWidth) * (targetRange[1] - targetRange[0])) :
                ((viewWidth / plotWidth) * targetLength);
            
            return {{
                viewLeft,
                viewRight,
                viewTop,
                viewBottom,
                viewWidth,
                viewHeight,
                viewportWidthBp
            }};
        }}

        function prepareLineData() {{
            lineVertices = [];
            lineColors = [];
            
            // Get consistent viewport info
            const currentViewport = calculateViewportInfo();
            const {{ viewLeft, viewRight, viewTop, viewBottom, viewWidth, viewHeight, viewportWidthBp }} = currentViewport;
            
            // Performance optimization: switch rendering modes based on viewport size
            const VIEWPORT_THRESHOLD = 2000000; // 2Mbp threshold for performance
            const shouldShowAlignmentSegments = viewportWidthBp > VIEWPORT_THRESHOLD;
            const shouldShowDetails = showDetails && !shouldShowAlignmentSegments && viewportWidthBp < 200000; // Show CIGAR details when zoomed in to 200kb or less
            
            // Add grid lines (sequence boundaries) - gray like the border
            const gridColor = hexToRgb(borderColor);
            const endLineColor = hexToRgb(borderColor);
            
            // Draw plot border
            const borderCoords = [
                projectCoords(0, 0),
                projectCoords(targetLength, 0),
                projectCoords(targetLength, queryLength),
                projectCoords(0, queryLength)
            ];
            
            // Top border
            lineVertices.push(borderCoords[0].x, borderCoords[0].y);
            lineVertices.push(borderCoords[1].x, borderCoords[1].y);
            lineColors.push(gridColor.r, gridColor.g, gridColor.b, 1.0);
            lineColors.push(gridColor.r, gridColor.g, gridColor.b, 1.0);
            
            // Right border
            lineVertices.push(borderCoords[1].x, borderCoords[1].y);
            lineVertices.push(borderCoords[2].x, borderCoords[2].y);
            lineColors.push(gridColor.r, gridColor.g, gridColor.b, 1.0);
            lineColors.push(gridColor.r, gridColor.g, gridColor.b, 1.0);
            
            // Bottom border
            lineVertices.push(borderCoords[2].x, borderCoords[2].y);
            lineVertices.push(borderCoords[3].x, borderCoords[3].y);
            lineColors.push(gridColor.r, gridColor.g, gridColor.b, 1.0);
            lineColors.push(gridColor.r, gridColor.g, gridColor.b, 1.0);
            
            // Left border
            lineVertices.push(borderCoords[3].x, borderCoords[3].y);
            lineVertices.push(borderCoords[0].x, borderCoords[0].y);
            lineColors.push(gridColor.r, gridColor.g, gridColor.b, 1.0);
            lineColors.push(gridColor.r, gridColor.g, gridColor.b, 1.0);
            
            // Target boundaries (vertical lines)
            targets.forEach(target => {{
                if (target.offset > 0) {{
                    const coords = projectCoords(target.offset, 0);
                    const endCoords = projectCoords(target.offset, queryLength);
                    
                    lineVertices.push(coords.x, coords.y);
                    lineVertices.push(endCoords.x, endCoords.y);
                    lineColors.push(gridColor.r, gridColor.g, gridColor.b, 1.0);
                    lineColors.push(gridColor.r, gridColor.g, gridColor.b, 1.0);
                }}
                
                // Add sequence end lines
                const endX = target.offset + target.length;
                const endStart = projectCoords(endX, 0);
                const endEnd = projectCoords(endX, queryLength);
                lineVertices.push(endStart.x, endStart.y);
                lineVertices.push(endEnd.x, endEnd.y);
                lineColors.push(endLineColor.r, endLineColor.g, endLineColor.b, 1.0);
                lineColors.push(endLineColor.r, endLineColor.g, endLineColor.b, 1.0);
            }});
            
            // Query boundaries (horizontal lines)
            queries.forEach(query => {{
                if (query.offset > 0) {{
                    const coords = projectCoords(0, query.offset);
                    const endCoords = projectCoords(targetLength, query.offset);
                    
                    lineVertices.push(coords.x, coords.y);
                    lineVertices.push(endCoords.x, endCoords.y);
                    lineColors.push(gridColor.r, gridColor.g, gridColor.b, 1.0);
                    lineColors.push(gridColor.r, gridColor.g, gridColor.b, 1.0);
                }}
                
                // Add sequence end lines
                const endY = query.offset + query.length;
                const endStart = projectCoords(0, endY);
                const endEnd = projectCoords(targetLength, endY);
                lineVertices.push(endStart.x, endStart.y);
                lineVertices.push(endEnd.x, endEnd.y);
                lineColors.push(endLineColor.r, endLineColor.g, endLineColor.b, 1.0);
                lineColors.push(endLineColor.r, endLineColor.g, endLineColor.b, 1.0);
            }});
            
            // Draw BED regions before alignments so alignments render on top
            drawBedRegions();
            
            // Draw alignments - three rendering modes based on viewport size
            const matchColor = hexToRgb(lineColor);
            const mismatchColor = hexToRgb(darkMode ? '#ff4444' : '#cc0000');
            const indelColor = hexToRgb(darkMode ? '#4444ff' : '#0000cc');
            const summaryColor = hexToRgb(lineColor); // Use same color as standard matches
            
            if (shouldShowAlignmentSegments) {{
                // Mode 1: Alignment segments view (viewport > 2Mbp)
                // Random sparsification for performance at extreme zoom levels
                let alignmentsToRender = alignments;
                const targetRenderCount = 100000; // Target number of alignments to render
                
                if (alignments.length > targetRenderCount) {{
                    // Calculate sampling probability
                    const sampleProbability = targetRenderCount / alignments.length;
                    
                    // Deterministic sparsification based on alignment hash
                    alignmentsToRender = alignments.filter(alignment => {{
                        // Create a simple hash from alignment features
                        const hashStr = `${{alignment.x}}_${{alignment.y}}_${{alignment.targetLen}}_${{alignment.queryLen}}_${{alignment.rev}}`;
                        let hash = 0;
                        for (let i = 0; i < hashStr.length; i++) {{
                            const char = hashStr.charCodeAt(i);
                            hash = ((hash << 5) - hash) + char;
                            hash = hash & hash; // Convert to 32-bit integer
                        }}
                        // Use hash to deterministically decide if this alignment should be rendered
                        const normalizedHash = (Math.abs(hash) % 1000000) / 1000000;
                        return normalizedHash < sampleProbability;
                    }});
                }}
                
                alignmentsToRender.forEach(alignment => {{
                    
                    const startCoords = projectCoords(alignment.x, alignment.y);
                    const endX = alignment.x + alignment.targetLen;
                    const endY = alignment.y + (alignment.rev ? -alignment.queryLen : alignment.queryLen);
                    const endCoords = projectCoords(endX, endY);
                    
                    // Enhanced frustum culling with margin for edge cases
                    const minX = Math.min(startCoords.x, endCoords.x);
                    const maxX = Math.max(startCoords.x, endCoords.x);
                    const minY = Math.min(startCoords.y, endCoords.y);
                    const maxY = Math.max(startCoords.y, endCoords.y);
                    
                    // Add margin to prevent edge case culling issues
                    const margin = 100; // pixels
                    if (maxX >= viewLeft - margin && minX <= viewRight + margin && 
                        maxY >= viewTop - margin && minY <= viewBottom + margin) {{
                        // Use consistent high opacity for visibility at all zoom levels
                        const alpha = 0.9;
                        lineVertices.push(startCoords.x, startCoords.y);
                        lineVertices.push(endCoords.x, endCoords.y);
                        lineColors.push(summaryColor.r, summaryColor.g, summaryColor.b, alpha);
                        lineColors.push(summaryColor.r, summaryColor.g, summaryColor.b, alpha);
                    }}
                }});
            }} else if (shouldShowDetails && detailedAlignments.length > 0) {{
                // Mode 2: Detailed CIGAR operations (viewport < 200kbp and details enabled)
                // Draw gray boxes for each alignment segment first
                const segmentColor = hexToRgb(darkMode ? '#333333' : '#e0e0e0');
                
                detailedAlignments.forEach(alignment => {{
                    // First check if the entire alignment is within or intersects the viewport
                    const segmentStartCoords = projectCoords(alignment.x, alignment.y);
                    const segmentEndX = alignment.x + alignment.targetLen;
                    const segmentEndY = alignment.y + (alignment.rev ? -alignment.queryLen : alignment.queryLen);
                    const segmentEndCoords = projectCoords(segmentEndX, segmentEndY);
                    
                    // Alignment bounding box in screen coordinates
                    const alignmentMinX = Math.min(segmentStartCoords.x, segmentEndCoords.x);
                    const alignmentMaxX = Math.max(segmentStartCoords.x, segmentEndCoords.x);
                    const alignmentMinY = Math.min(segmentStartCoords.y, segmentEndCoords.y);
                    const alignmentMaxY = Math.max(segmentStartCoords.y, segmentEndCoords.y);
                    
                    // Check if alignment intersects viewport (with generous margins for partial visibility)
                    const margin = 50; // Allow rendering of partially visible alignments
                    const alignmentVisible = (
                        alignmentMaxX >= viewLeft - margin && 
                        alignmentMinX <= viewRight + margin && 
                        alignmentMaxY >= viewTop - margin && 
                        alignmentMinY <= viewBottom + margin
                    );
                    
                    if (!alignmentVisible) {{
                        return; // Skip this entire alignment if it's not visible
                    }}
                    
                    // Draw segment boundary box
                    const boxCorners = [
                        projectCoords(alignment.x, alignment.y),
                        projectCoords(alignment.x + alignment.targetLen, alignment.y),
                        projectCoords(alignment.x + alignment.targetLen, alignment.y + (alignment.rev ? -alignment.queryLen : alignment.queryLen)),
                        projectCoords(alignment.x, alignment.y + (alignment.rev ? -alignment.queryLen : alignment.queryLen))
                    ];
                    
                    // Top edge
                    lineVertices.push(boxCorners[0].x, boxCorners[0].y);
                    lineVertices.push(boxCorners[1].x, boxCorners[1].y);
                    lineColors.push(segmentColor.r, segmentColor.g, segmentColor.b, 0.6);
                    lineColors.push(segmentColor.r, segmentColor.g, segmentColor.b, 0.6);
                    
                    // Right edge
                    lineVertices.push(boxCorners[1].x, boxCorners[1].y);
                    lineVertices.push(boxCorners[2].x, boxCorners[2].y);
                    lineColors.push(segmentColor.r, segmentColor.g, segmentColor.b, 0.6);
                    lineColors.push(segmentColor.r, segmentColor.g, segmentColor.b, 0.6);
                    
                    // Bottom edge
                    lineVertices.push(boxCorners[2].x, boxCorners[2].y);
                    lineVertices.push(boxCorners[3].x, boxCorners[3].y);
                    lineColors.push(segmentColor.r, segmentColor.g, segmentColor.b, 0.6);
                    lineColors.push(segmentColor.r, segmentColor.g, segmentColor.b, 0.6);
                    
                    // Left edge
                    lineVertices.push(boxCorners[3].x, boxCorners[3].y);
                    lineVertices.push(boxCorners[0].x, boxCorners[0].y);
                    lineColors.push(segmentColor.r, segmentColor.g, segmentColor.b, 0.6);
                    lineColors.push(segmentColor.r, segmentColor.g, segmentColor.b, 0.6);
                    
                    // Now draw detailed CIGAR operations within the box - with individual viewport culling
                    alignment.ops.forEach(op => {{
                        // For reverse alignments, we need to transform both X and Y coordinates
                        // The stored op.y values decrement from alignment.y (query end)
                        // But to match PNG rendering, we need Y to increment from query start
                        const opY = op.rev ? (alignment.y - alignment.queryLen + (alignment.y - op.y)) : op.y;
                        // For reverse alignments, X needs to go backwards from the end
                        const opX = op.rev ? (alignment.x + alignment.targetLen - (op.x - alignment.x)) : op.x;
                        const startCoords = projectCoords(opX, opY);
                        
                        if (op.type === 'match' || op.type === 'unknown') {{
                            const endX = opX + (op.rev ? -op.len : op.len);
                            const endY = opY + op.len;
                            const endCoords = projectCoords(endX, endY);
                            
                            // Check if this match operation intersects the viewport
                            const opMinX = Math.min(startCoords.x, endCoords.x);
                            const opMaxX = Math.max(startCoords.x, endCoords.x);
                            const opMinY = Math.min(startCoords.y, endCoords.y);
                            const opMaxY = Math.max(startCoords.y, endCoords.y);
                            
                            if (opMaxX >= viewLeft && opMinX <= viewRight && opMaxY >= viewTop && opMinY <= viewBottom) {{
                                lineVertices.push(startCoords.x, startCoords.y);
                                lineVertices.push(endCoords.x, endCoords.y);
                                lineColors.push(matchColor.r, matchColor.g, matchColor.b, 0.9);
                                lineColors.push(matchColor.r, matchColor.g, matchColor.b, 0.9);
                            }}
                        }} else if (op.type === 'mismatch') {{
                            // Check if mismatch region intersects viewport before rendering individual segments
                            const mismatchEndX = opX + (op.rev ? -op.len : op.len);
                            const mismatchEndY = opY + op.len;
                            const mismatchEndCoords = projectCoords(mismatchEndX, mismatchEndY);
                            
                            const mismatchMinX = Math.min(startCoords.x, mismatchEndCoords.x);
                            const mismatchMaxX = Math.max(startCoords.x, mismatchEndCoords.x);
                            const mismatchMinY = Math.min(startCoords.y, mismatchEndCoords.y);
                            const mismatchMaxY = Math.max(startCoords.y, mismatchEndCoords.y);
                            
                            if (mismatchMaxX >= viewLeft && mismatchMinX <= viewRight && mismatchMaxY >= viewTop && mismatchMinY <= viewBottom) {{
                                // Draw mismatches as individual segments, but limit count for performance
                                const maxSegments = Math.min(op.len, 100); // Limit to prevent too many segments
                                for (let i = 0; i < maxSegments; i++) {{
                                    const mx = opX + (op.rev ? -i : i);
                                    const my = opY + i;
                                    const mStartCoords = projectCoords(mx, my);
                                    const mEndCoords = projectCoords(mx + (op.rev ? -1 : 1), my + 1);
                                    
                                    // Individual segment viewport check for very fine-grained culling
                                    const segMinX = Math.min(mStartCoords.x, mEndCoords.x);
                                    const segMaxX = Math.max(mStartCoords.x, mEndCoords.x);
                                    const segMinY = Math.min(mStartCoords.y, mEndCoords.y);
                                    const segMaxY = Math.max(mStartCoords.y, mEndCoords.y);
                                    
                                    if (segMaxX >= viewLeft && segMinX <= viewRight && segMaxY >= viewTop && segMinY <= viewBottom) {{
                                        lineVertices.push(mStartCoords.x, mStartCoords.y);
                                        lineVertices.push(mEndCoords.x, mEndCoords.y);
                                        lineColors.push(mismatchColor.r, mismatchColor.g, mismatchColor.b, 1.0);
                                        lineColors.push(mismatchColor.r, mismatchColor.g, mismatchColor.b, 1.0);
                                    }}
                                }}
                            }}
                        }} else if (op.type === 'insertion') {{
                            const endCoords = projectCoords(opX, opY + op.len);
                            
                            // Check if insertion line intersects viewport
                            const insMinX = Math.min(startCoords.x, endCoords.x);
                            const insMaxX = Math.max(startCoords.x, endCoords.x);
                            const insMinY = Math.min(startCoords.y, endCoords.y);
                            const insMaxY = Math.max(startCoords.y, endCoords.y);
                            
                            if (insMaxX >= viewLeft && insMinX <= viewRight && insMaxY >= viewTop && insMinY <= viewBottom) {{
                                lineVertices.push(startCoords.x, startCoords.y);
                                lineVertices.push(endCoords.x, endCoords.y);
                                lineColors.push(indelColor.r, indelColor.g, indelColor.b, 0.9);
                                lineColors.push(indelColor.r, indelColor.g, indelColor.b, 0.9);
                            }}
                        }} else if (op.type === 'deletion') {{
                            const endCoords = projectCoords(opX + (op.rev ? -op.len : op.len), opY);
                            
                            // Check if deletion line intersects viewport
                            const delMinX = Math.min(startCoords.x, endCoords.x);
                            const delMaxX = Math.max(startCoords.x, endCoords.x);
                            const delMinY = Math.min(startCoords.y, endCoords.y);
                            const delMaxY = Math.max(startCoords.y, endCoords.y);
                            
                            if (delMaxX >= viewLeft && delMinX <= viewRight && delMaxY >= viewTop && delMinY <= viewBottom) {{
                                lineVertices.push(startCoords.x, startCoords.y);
                                lineVertices.push(endCoords.x, endCoords.y);
                                lineColors.push(indelColor.r, indelColor.g, indelColor.b, 0.9);
                                lineColors.push(indelColor.r, indelColor.g, indelColor.b, 0.9);
                            }}
                        }}
                    }});
                }});
            }} else {{
                // Mode 3: Standard alignment lines (200kbp - 2Mbp viewport)
                alignments.forEach(alignment => {{
                    const startCoords = projectCoords(alignment.x, alignment.y);
                    const endX = alignment.x + alignment.targetLen;
                    const endY = alignment.y + (alignment.rev ? -alignment.queryLen : alignment.queryLen);
                    const endCoords = projectCoords(endX, endY);
                    
                    // Frustum culling with margin
                    const minX = Math.min(startCoords.x, endCoords.x);
                    const maxX = Math.max(startCoords.x, endCoords.x);
                    const minY = Math.min(startCoords.y, endCoords.y);
                    const maxY = Math.max(startCoords.y, endCoords.y);
                    
                    const margin = 100; // pixels
                    if (maxX >= viewLeft - margin && minX <= viewRight + margin && 
                        maxY >= viewTop - margin && minY <= viewBottom + margin) {{
                        lineVertices.push(startCoords.x, startCoords.y);
                        lineVertices.push(endCoords.x, endCoords.y);
                        lineColors.push(matchColor.r, matchColor.g, matchColor.b, 0.9);
                        lineColors.push(matchColor.r, matchColor.g, matchColor.b, 0.9);
                    }}
                }});
            }}
            
        }}

        function drawBedRegions() {{
            // BED rendering - draw filled rectangles
            bedRegions.forEach(region => {{
                let targetSeq = targets.find(t => t.name === region.chr);
                let querySeq = queries.find(q => q.name === region.chr);
                
                // Use red color for BED regions (we can add color parsing later)
                const color = {{ r: 1, g: 0, b: 0, a: 0.4 }};
                
                if (targetSeq) {{
                    // Draw filled vertical rectangle for target sequence
                    const startX = targetSeq.offset + region.start;
                    const endX = targetSeq.offset + region.end;
                    
                    // Calculate step size based on both zoom and base pair coverage
                    // At high zoom, use base pair resolution; at low zoom, use pixel resolution
                    const currentViewport = calculateViewportInfo();
                    const bpPerPixel = currentViewport.viewportWidthBp / canvasWidth;
                    const stepSize = Math.min(Math.max(1, Math.floor(bpPerPixel)), Math.max(1, (endX - startX) / 1000));
                    
                    for (let x = startX; x <= endX; x += stepSize) {{
                        const top = projectCoords(x, 0);
                        const bottom = projectCoords(x, queryLength);
                        lineVertices.push(top.x, top.y);
                        lineVertices.push(bottom.x, bottom.y);
                        lineColors.push(color.r, color.g, color.b, color.a);
                        lineColors.push(color.r, color.g, color.b, color.a);
                    }}
                }}
                
                if (querySeq) {{
                    // Draw filled horizontal rectangle for query sequence
                    const startY = querySeq.offset + region.start;
                    const endY = querySeq.offset + region.end;
                    
                    // Calculate step size based on both zoom and base pair coverage
                    // At high zoom, use base pair resolution; at low zoom, use pixel resolution
                    const currentViewport = calculateViewportInfo();
                    const bpPerPixel = currentViewport.viewportWidthBp / canvasWidth;
                    const stepSize = Math.min(Math.max(1, Math.floor(bpPerPixel)), Math.max(1, (endY - startY) / 1000));
                    
                    for (let y = startY; y <= endY; y += stepSize) {{
                        const left = projectCoords(0, y);
                        const right = projectCoords(targetLength, y);
                        lineVertices.push(left.x, left.y);
                        lineVertices.push(right.x, right.y);
                        lineColors.push(color.r, color.g, color.b, color.a);
                        lineColors.push(color.r, color.g, color.b, color.a);
                    }}
                }}
            }});
        }}
        
        function drawPlot() {{
            if (!gl) return;
            
            // Handle canvas size
            canvas.width = canvasWidth;
            canvas.height = canvasHeight;
            
            gl.viewport(0, 0, canvasWidth, canvasHeight);
            gl.clear(gl.COLOR_BUFFER_BIT);
            
            if (lineVertices.length === 0) return;
            
            // Calculate dynamic line width based on zoom and user setting
            const baseLineWidth = 1.0;
            const zoomFactor = Math.max(1, Math.log2(zoom + 1));
            const lineWidth = Math.min(10, baseLineWidth * zoomFactor * lineWidthMultiplier);
            gl.lineWidth(lineWidth);
            
            // Use shader program
            gl.useProgram(shaderProgram);
            
            // Set uniforms
            const resolutionLoc = gl.getUniformLocation(shaderProgram, 'u_resolution');
            const panLoc = gl.getUniformLocation(shaderProgram, 'u_pan');
            const zoomLoc = gl.getUniformLocation(shaderProgram, 'u_zoom');
            
            gl.uniform2f(resolutionLoc, canvasWidth, canvasHeight);
            gl.uniform2f(panLoc, panX, panY);
            gl.uniform1f(zoomLoc, zoom);
            
            // Bind position buffer
            const positionLoc = gl.getAttribLocation(shaderProgram, 'a_position');
            gl.bindBuffer(gl.ARRAY_BUFFER, positionBuffer);
            gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(lineVertices), gl.DYNAMIC_DRAW);
            gl.enableVertexAttribArray(positionLoc);
            gl.vertexAttribPointer(positionLoc, 2, gl.FLOAT, false, 0, 0);
            
            // Bind color buffer
            const colorLoc = gl.getAttribLocation(shaderProgram, 'a_color');
            gl.bindBuffer(gl.ARRAY_BUFFER, colorBuffer);
            gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(lineColors), gl.DYNAMIC_DRAW);
            gl.enableVertexAttribArray(colorLoc);
            gl.vertexAttribPointer(colorLoc, 4, gl.FLOAT, false, 0, 0);
            
            // Draw lines
            gl.drawArrays(gl.LINES, 0, lineVertices.length / 2);
            
            // Draw selection box if active
            if (isBoxSelecting) {{
                drawSelectionBox();
            }}
        }}
        
        function drawSelectionBox() {{
            if (!isBoxSelecting) return;
            
            // Draw selection box using 2D canvas overlay
            let canvas2d = document.getElementById('overlayCanvas');
            if (!canvas2d) {{
                // Create overlay canvas if it doesn't exist
                const overlay = document.createElement('canvas');
                overlay.id = 'overlayCanvas';
                overlay.style.position = 'absolute';
                overlay.style.top = '0';
                overlay.style.left = '0';
                overlay.style.width = '100%';
                overlay.style.height = '100%';
                overlay.style.pointerEvents = 'none';
                overlay.style.zIndex = '10';
                overlay.width = canvasWidth;
                overlay.height = canvasHeight;
                document.body.appendChild(overlay);
                canvas2d = overlay;
            }}
            
            const ctx2d = canvas2d.getContext('2d');
            ctx2d.clearRect(0, 0, canvasWidth, canvasHeight);
            
            // Draw the selection box
            ctx2d.strokeStyle = darkMode ? '#00ff00' : '#0088ff';
            ctx2d.lineWidth = 2;
            ctx2d.setLineDash([5, 5]);
            ctx2d.strokeRect(
                Math.min(boxStartX, boxEndX),
                Math.min(boxStartY, boxEndY),
                Math.abs(boxEndX - boxStartX),
                Math.abs(boxEndY - boxStartY)
            );
            
            // Fill with semi-transparent color
            ctx2d.fillStyle = darkMode ? 'rgba(0, 255, 0, 0.1)' : 'rgba(0, 136, 255, 0.1)';
            ctx2d.fillRect(
                Math.min(boxStartX, boxEndX),
                Math.min(boxStartY, boxEndY),
                Math.abs(boxEndX - boxStartX),
                Math.abs(boxEndY - boxStartY)
            );
        }}

        function pixelToGenomicCoords(pixelX, pixelY) {{
            // Convert screen coordinates to plot coordinates
            const plotX = (pixelX - panX) / zoom;
            const plotY = (pixelY - panY) / zoom;
            
            // Remove plot offset to get relative position in plot
            const relX = plotX - plotOffsetX;
            const relY = plotY - plotOffsetY;
            
            // Convert to genomic coordinates
            let genomicX, genomicY;
            
            if (usingZoom) {{
                genomicX = targetRange[0] + (relX / plotWidth) * (targetRange[1] - targetRange[0]);
                genomicY = queryRange[0] + ((plotHeight - relY) / plotHeight) * (queryRange[1] - queryRange[0]);
            }} else {{
                genomicX = (relX / plotWidth) * targetLength;
                genomicY = ((plotHeight - relY) / plotHeight) * queryLength;
            }}
            
            return {{
                target: Math.max(0, Math.min(targetLength, Math.round(genomicX))),
                query: Math.max(0, Math.min(queryLength, Math.round(genomicY)))
            }};
        }}

        function findSequenceAtCoord(coord, sequences, totalLength) {{
            for (const seq of sequences) {{
                if (coord >= seq.offset && coord < seq.offset + seq.length) {{
                    return {{
                        sequence: seq,
                        relativePos: coord - seq.offset
                    }};
                }}
            }}
            return null;
        }}


        function drawMinimap() {{
            // Set minimap size
            minimapCanvas.width = 200;
            minimapCanvas.height = 150;
            
            // Clear and draw background
            minimapCtx.fillStyle = backgroundColor;
            minimapCtx.fillRect(0, 0, 200, 150);
            
            // Draw simplified plot
            minimapCtx.strokeStyle = borderColor;
            minimapCtx.lineWidth = 0.5;
            
            // Draw a simplified version of alignments
            if (alignments.length < 5000) {{
                minimapCtx.globalAlpha = 0.3;
                minimapCtx.strokeStyle = lineColor;
                alignments.forEach(alignment => {{
                    const startCoords = projectCoords(alignment.x, alignment.y);
                    const endX = alignment.x + (alignment.rev ? -alignment.targetLen : alignment.targetLen);
                    const endY = alignment.y + alignment.queryLen;
                    const endCoords = projectCoords(endX, endY);
                    
                    minimapCtx.beginPath();
                    minimapCtx.moveTo(startCoords.x * 200 / canvasWidth, startCoords.y * 150 / canvasHeight);
                    minimapCtx.lineTo(endCoords.x * 200 / canvasWidth, endCoords.y * 150 / canvasHeight);
                    minimapCtx.stroke();
                }});
            }} else {{
                // Draw grid for large datasets
                minimapCtx.fillStyle = darkMode ? 'rgba(255,255,255,0.1)' : 'rgba(0,0,0,0.1)';
                for (let i = 0; i < 20; i++) {{
                    for (let j = 0; j < 15; j++) {{
                        if ((i + j) % 2 === 0) {{
                            minimapCtx.fillRect(i * 10, j * 10, 10, 10);
                        }}
                    }}
                }}
            }}
            
            // Update viewport rectangle
            const viewLeft = Math.max(0, -panX) / zoom;
            const viewTop = Math.max(0, -panY) / zoom;
            const viewWidth = Math.min(canvasWidth, canvasWidth / zoom);
            const viewHeight = Math.min(canvasHeight, canvasHeight / zoom);
            
            const vpLeft = (viewLeft / canvasWidth) * 200;
            const vpTop = (viewTop / canvasHeight) * 150;
            const vpWidth = (viewWidth / canvasWidth) * 200;
            const vpHeight = (viewHeight / canvasHeight) * 150;
            
            minimapViewport.style.left = vpLeft + 'px';
            minimapViewport.style.top = vpTop + 'px';
            minimapViewport.style.width = vpWidth + 'px';
            minimapViewport.style.height = vpHeight + 'px';
        }}

        function updateLabels() {{
            if (!showLabels) {{
                targetLabels.style.display = 'none';
                queryLabels.style.display = 'none';
                return;
            }}
            
            targetLabels.style.display = 'block';
            queryLabels.style.display = 'block';
            
            // Clear existing labels
            targetLabels.innerHTML = '';
            queryLabels.innerHTML = '';
            
            // Calculate visible plot bounds
            const viewLeft = Math.max(plotOffsetX, (-panX) / zoom);
            const viewRight = Math.min(plotOffsetX + plotWidth, (canvasWidth - panX) / zoom);
            const viewTop = Math.max(plotOffsetY, (-panY) / zoom);
            const viewBottom = Math.min(plotOffsetY + plotHeight, (canvasHeight - panY) / zoom);
            
            // Target labels (bottom)
            targets.forEach(target => {{
                const startCoords = projectCoords(target.offset, 0);
                const endCoords = projectCoords(target.offset + target.length, 0);
                const startX = startCoords.x * zoom + panX;
                const endX = endCoords.x * zoom + panX;
                
                // Check if sequence is at least partially visible
                const seqVisible = (startX <= canvasWidth && endX >= 0);
                
                if (seqVisible) {{
                    const label = document.createElement('div');
                    label.className = 'label target-label';
                    label.textContent = target.name;
                    label.title = `${{target.name}} (${{target.length.toLocaleString()}} bp)`;
                    
                    // Calculate label position, clamping to viewport edges
                    let labelX = startX;
                    
                    // If sequence extends beyond left edge, stick label to left border
                    if (startX < 0 && endX > 0) {{
                        labelX = 5; // Small margin from edge
                        label.style.borderLeft = '3px solid ' + (darkMode ? '#00ff00' : '#0088ff');
                    }}
                    // If sequence starts beyond right edge but is visible, stick to right
                    else if (startX > canvasWidth - 100 && startX < canvasWidth + 500) {{
                        labelX = canvasWidth - 100;
                        label.style.borderRight = '3px solid ' + (darkMode ? '#00ff00' : '#0088ff');
                    }}
                    
                    label.style.left = labelX + 'px';
                    label.style.bottom = '10px';
                    targetLabels.appendChild(label);
                }}
            }});
            
            // Query labels (right side)
            queries.forEach(query => {{
                const startCoords = projectCoords(0, query.offset + query.length);
                const endCoords = projectCoords(0, query.offset);
                const startY = startCoords.y * zoom + panY;
                const endY = endCoords.y * zoom + panY;
                
                // Check if sequence is at least partially visible
                const seqVisible = (Math.min(startY, endY) <= canvasHeight && Math.max(startY, endY) >= 0);
                
                if (seqVisible) {{
                    const label = document.createElement('div');
                    label.className = 'label query-label';
                    label.textContent = query.name;
                    label.title = `${{query.name}} (${{query.length.toLocaleString()}} bp)`;
                    
                    // Calculate label position, clamping to viewport edges
                    let labelY = startY;
                    
                    // If sequence extends beyond top edge, stick label to top border
                    if (startY < 0 && endY > 0) {{
                        labelY = 5; // Small margin from edge
                        label.style.borderTop = '3px solid ' + (darkMode ? '#00ff00' : '#0088ff');
                    }}
                    // If sequence starts beyond bottom edge but is visible, stick to bottom
                    else if (startY > canvasHeight - 50 && startY < canvasHeight + 200) {{
                        labelY = canvasHeight - 30;
                        label.style.borderBottom = '3px solid ' + (darkMode ? '#00ff00' : '#0088ff');
                    }}
                    // For sequences within viewport, clamp to visible area
                    else if (endY > canvasHeight) {{
                        // If sequence end is off screen, position label at a visible spot
                        labelY = Math.min(startY, canvasHeight - 30);
                    }}
                    
                    label.style.top = labelY + 'px';
                    queryLabels.appendChild(label);
                }}
            }});
        }}

        function updateDisplay() {{
            prepareLineData();  // Rebuild line data with frustum culling
            drawPlot();
            drawMinimap();
            updateLabels();
            
            // Get consistent viewport info
            const currentViewport = calculateViewportInfo();
            const {{ viewportWidthBp }} = currentViewport;
            
            // Update rendering mode indicator
            let mode = 'Standard';
            if (viewportWidthBp > 2000000) {{
                mode = 'Segments';
            }} else if (showDetails && viewportWidthBp < 200000) {{
                mode = 'Details';
            }}
            renderingMode.textContent = mode;
            
            // Update zoom info
            zoomLevel.textContent = zoom.toFixed(2) + 'x';
            panInfo.textContent = `${{Math.round(panX)}}, ${{Math.round(panY)}}`;
            viewportInfo.textContent = viewportWidthBp >= 1000 ? 
                `${{Math.round(viewportWidthBp / 1000)}}k bp` : 
                `${{Math.round(viewportWidthBp)}} bp`;
        }}

        // Mouse event handlers
        canvas.addEventListener('mousemove', function(e) {{
            const rect = canvas.getBoundingClientRect();
            const pixelX = e.clientX - rect.left;
            const pixelY = e.clientY - rect.top;
            
            if (isDragging) {{
                panX = dragStartPanX + (pixelX - dragStartX);
                panY = dragStartPanY + (pixelY - dragStartY);
                updateDisplay();
                return;
            }}
            
            if (isBoxSelecting) {{
                boxEndX = pixelX;
                boxEndY = pixelY;
                updateDisplay();
                return;
            }}
            
            // Update tooltip and info
            const genomicCoords = pixelToGenomicCoords(pixelX, pixelY);
            const targetSeq = findSequenceAtCoord(genomicCoords.target, targets, targetLength);
            const querySeq = findSequenceAtCoord(genomicCoords.query, queries, queryLength);
            
            const targetText = targetSeq ? 
                `${{targetSeq.sequence.name}}:${{targetSeq.relativePos.toLocaleString()}}` : 
                `${{genomicCoords.target.toLocaleString()}}`;
            const queryText = querySeq ? 
                `${{querySeq.sequence.name}}:${{querySeq.relativePos.toLocaleString()}}` : 
                `${{genomicCoords.query.toLocaleString()}}`;
            
            position.innerHTML = `Target: ${{targetText}}<br>Query: ${{queryText}}`;
            
            // Update sequence info
            let info = '';
            if (targetSeq) {{
                info += `<strong>Target:</strong> ${{targetSeq.sequence.name}}<br>`;
                info += `Length: ${{targetSeq.sequence.length.toLocaleString()}} bp<br>`;
                info += `Position: ${{targetSeq.relativePos.toLocaleString()}} / ${{targetSeq.sequence.length.toLocaleString()}}<br><br>`;
            }}
            if (querySeq) {{
                info += `<strong>Query:</strong> ${{querySeq.sequence.name}}<br>`;
                info += `Length: ${{querySeq.sequence.length.toLocaleString()}} bp<br>`;
                info += `Position: ${{querySeq.relativePos.toLocaleString()}} / ${{querySeq.sequence.length.toLocaleString()}}`;
            }}
            sequenceInfo.innerHTML = info || 'No sequence information available';
            
            // Show tooltip
            tooltip.style.display = 'block';
            tooltip.style.left = (rect.left + pixelX + 10) + 'px';
            tooltip.style.top = (rect.top + pixelY - 30) + 'px';
            tooltip.innerHTML = `${{targetText}}<br>${{queryText}}`;
        }});

        canvas.addEventListener('mouseleave', function() {{
            tooltip.style.display = 'none';
        }});

        canvas.addEventListener('mousedown', function(e) {{
            e.preventDefault();
            const rect = canvas.getBoundingClientRect();
            const mouseX = e.clientX - rect.left;
            const mouseY = e.clientY - rect.top;
            
            if (e.button === 0) {{ // Left click - pan
                isDragging = true;
                dragStartX = mouseX;
                dragStartY = mouseY;
                dragStartPanX = panX;
                dragStartPanY = panY;
                canvas.style.cursor = 'grabbing';
            }} else if (e.button === 2) {{ // Right click - zoom box
                isBoxSelecting = true;
                boxStartX = mouseX;
                boxStartY = mouseY;
                boxEndX = mouseX;
                boxEndY = mouseY;
                canvas.style.cursor = 'crosshair';
            }}
        }});

        canvas.addEventListener('mouseup', function(e) {{
            if (e.button === 0 && isDragging) {{ // Left click release
                isDragging = false;
                canvas.style.cursor = 'crosshair';
            }} else if (e.button === 2 && isBoxSelecting) {{ // Right click release
                isBoxSelecting = false;
                canvas.style.cursor = 'crosshair';
                
                // Zoom to selected box
                const boxWidth = Math.abs(boxEndX - boxStartX);
                const boxHeight = Math.abs(boxEndY - boxStartY);
                
                if (boxWidth > 10 && boxHeight > 10) {{
                    zoomToBox(
                        Math.min(boxStartX, boxEndX),
                        Math.min(boxStartY, boxEndY),
                        boxWidth,
                        boxHeight
                    );
                }} else {{
                    // Clear the overlay if box was too small
                    const canvas2d = document.getElementById('overlayCanvas');
                    if (canvas2d) {{
                        const ctx2d = canvas2d.getContext('2d');
                        ctx2d.clearRect(0, 0, canvasWidth, canvasHeight);
                    }}
                }}
                
                updateDisplay();
            }}
        }});

        canvas.addEventListener('contextmenu', function(e) {{
            e.preventDefault(); // Prevent context menu
        }});

        canvas.addEventListener('wheel', function(e) {{
            e.preventDefault();
            e.stopPropagation();
            
            const rect = canvas.getBoundingClientRect();
            const mouseX = e.clientX - rect.left;
            const mouseY = e.clientY - rect.top;
            
            // More responsive zoom factors
            const zoomFactor = e.deltaY > 0 ? 0.9 : 1.1;
            const newZoom = Math.max(0.1, Math.min(100000, zoom * zoomFactor));
            
            if (newZoom !== zoom) {{
                // Save current state before zooming
                saveZoomState();
                
                // Calculate the point under the mouse in canvas space (before zoom)
                const pointX = (mouseX - panX) / zoom;
                const pointY = (mouseY - panY) / zoom;
                
                // Update zoom
                zoom = newZoom;
                
                // Adjust pan so the same point stays under the mouse
                panX = mouseX - pointX * zoom;
                panY = mouseY - pointY * zoom;
                
                updateDisplay();
            }}
        }}, {{ passive: false }});

        function zoomToBox(boxX, boxY, boxWidth, boxHeight) {{
            // Save current state before zooming
            saveZoomState();
            
            // Calculate the center of the box in screen coordinates
            const boxCenterX = boxX + boxWidth / 2;
            const boxCenterY = boxY + boxHeight / 2;
            
            // Get current center point in plot coordinates (before zoom change)
            const centerPlotX = (boxCenterX - panX) / zoom;
            const centerPlotY = (boxCenterY - panY) / zoom;
            
            // Calculate new zoom to fit the box
            const zoomFactorX = canvasWidth / boxWidth;
            const zoomFactorY = canvasHeight / boxHeight;
            const newZoom = Math.min(zoomFactorX, zoomFactorY) * zoom * 0.9; // 0.9 for padding
            
            // Clamp zoom
            zoom = Math.max(0.1, Math.min(100000, newZoom));
            
            // Calculate new pan to keep the box center in the screen center
            panX = canvasWidth / 2 - centerPlotX * zoom;
            panY = canvasHeight / 2 - centerPlotY * zoom;
            
            // Clear the overlay canvas
            const canvas2d = document.getElementById('overlayCanvas');
            if (canvas2d) {{
                const ctx2d = canvas2d.getContext('2d');
                ctx2d.clearRect(0, 0, canvasWidth, canvasHeight);
            }}
            
            updateDisplay();
        }}

        function resetZoom() {{
            zoom = 1.0;
            panX = 0;
            panY = 0;
            updateDisplay();
        }}

        // Configuration event handlers
        lineWidthSlider.addEventListener('input', function(e) {{
            lineWidthMultiplier = parseFloat(e.target.value);
            lineWidthValue.textContent = lineWidthMultiplier.toFixed(1) + 'x';
            updateDisplay();
        }});
        
        gridOpacitySlider.addEventListener('input', function(e) {{
            gridOpacity = parseFloat(e.target.value);
            gridOpacityValue.textContent = Math.round(gridOpacity * 100) + '%';
            updateDisplay();
        }});
        
        showDetailsCheckbox.addEventListener('change', function(e) {{
            showDetails = e.target.checked;
            updateDisplay();
        }});
        
        showLabelsCheckbox.addEventListener('change', function(e) {{
            showLabels = e.target.checked;
            updateDisplay();
        }});
        
        // Handle window resize
        window.addEventListener('resize', function() {{
            canvasWidth = window.innerWidth;
            canvasHeight = window.innerHeight;
            calculatePlotDimensions();
            updateDisplay();
        }});
        
        // Initialize
        if (initWebGL()) {{
            calculatePlotDimensions();
            prepareLineData();
            updateDisplay();
            updateUndoButton();
        }} else {{
            document.body.innerHTML = '<h1>WebGL not supported</h1>';
        }}
    </script>
</body>
</html>"#,
        config.output_filename,
        if config.dark { "#1a1a1a" } else { "#ffffff" },
        if config.dark { "#ffffff" } else { "#000000" },
        if config.dark { "#444444" } else { "#cccccc" },
        if config.dark { "#2a2a2a" } else { "#f5f5f5" },
        if config.dark { "#444444" } else { "#cccccc" },
        if config.dark { "#2a2a2a" } else { "#f5f5f5" }, // config panel background
        if config.dark { "#444444" } else { "#cccccc" }, // config panel border
        if config.dark { "#2a2a2a" } else { "#f0f0f0" }, // minimap background
        if config.dark { "#666666" } else { "#999999" }, // minimap border
        if config.dark { "#00ff00" } else { "#0088ff" }, // viewport border
        if config.dark {
            "rgba(0,255,0,0.2)"
        } else {
            "rgba(0,136,255,0.2)"
        }, // viewport fill
        if config.dark { "#666666" } else { "#999999" },
        if config.dark {
            "rgba(26,26,26,0.8)"
        } else {
            "rgba(255,255,255,0.8)"
        }, // label background
        if config.dark { "#1a1a1a" } else { "#ffffff" }, // page background again
        if config.dark { "#00ff00" } else { "#0088ff" }, // render mode color
        alignments_json,
        summary_json,
        detailed_json,
        targets_json,
        queries_json,
        bed_json,
        bedpe_json,
        if config.using_zoom {
            config.ranges.0 .1 - config.ranges.0 .0
        } else {
            config.paf.target_length as usize
        },
        if config.using_zoom {
            config.ranges.1 .1 - config.ranges.1 .0
        } else {
            config.paf.query_length as usize
        },
        config.using_zoom,
        format!("[{}, {}]", config.ranges.0 .0, config.ranges.0 .1),
        format!("[{}, {}]", config.ranges.1 .0, config.ranges.1 .1),
        config.dark,
    );

    let html_filename = if config.output_filename.ends_with(".html") {
        config.output_filename.to_string()
    } else if config.output_filename.ends_with(".png") {
        config.output_filename.replace(".png", ".html")
    } else {
        format!("{}.html", config.output_filename)
    };
    std::fs::write(&html_filename, html_content).expect("Failed to write HTML file");

    println!("Generated interactive HTML viewer: {html_filename}");
}

/// Per-(query, target) pair statistics used to fill PDF cell tooltips:
/// (alignment count, total aligned bp, best identity 0..1)
type PairStats = std::collections::HashMap<(String, String), (usize, u64, f64)>;

fn find_seq<'a>(seqs: &'a [AlignedSeq], name: &str) -> Option<&'a AlignedSeq> {
    seqs.iter().find(|s| s.name == name)
}

/// A 2D highlight region drawn on the PDF dot plot: query row range x target
/// column range, with an optional color (hex, default orange) and label.
#[derive(Debug, Clone)]
struct HighlightRegion {
    qname: String,
    qstart: usize,
    qend: usize,
    tname: String,
    tstart: usize,
    tend: usize,
    color: (f64, f64, f64),
    label: Option<String>,
}

fn parse_hex_color(s: &str) -> Option<(f64, f64, f64)> {
    let h = s.trim_start_matches('#');
    if h.len() == 6 {
        let r = u8::from_str_radix(&h[0..2], 16).ok()?;
        let g = u8::from_str_radix(&h[2..4], 16).ok()?;
        let b = u8::from_str_radix(&h[4..6], 16).ok()?;
        Some((r as f64 / 255.0, g as f64 / 255.0, b as f64 / 255.0))
    } else {
        None
    }
}

fn parse_highlight_file(filename: &str) -> Vec<HighlightRegion> {
    let mut regions = Vec::new();
    let file = File::open(filename).unwrap();
    let reader: Box<dyn BufRead> = if filename.ends_with(".gz") {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 6 {
            eprintln!("[pafplot::highlight] skipping malformed line: {}", line);
            continue;
        }
        let color = if f.len() > 6 {
            parse_hex_color(f[6]).unwrap_or((1.0, 0.55, 0.0))
        } else {
            (1.0, 0.55, 0.0)
        };
        let label = if f.len() > 7 {
            Some(f[7..].join(" "))
        } else {
            None
        };
        regions.push(HighlightRegion {
            qname: f[0].into(),
            qstart: f[1].parse().unwrap(),
            qend: f[2].parse().unwrap(),
            tname: f[3].into(),
            tstart: f[4].parse().unwrap(),
            tend: f[5].parse().unwrap(),
            color,
            label,
        });
    }
    regions
}

struct PdfViewerConfig<'a> {
    paf: &'a PafFile,
    axes: (usize, usize),
    raw: &'a [u8],
    output_filename: &'a str,
    using_zoom: bool,
    ranges: ((usize, usize), (usize, usize)),
    pair_stats: &'a PairStats,
    highlights: &'a [HighlightRegion],
    // identity-color legend: (min, max, inverted) when line coloring is on
    legend: Option<(f64, f64, bool)>,
    // PDF label font size in points (labels shrink with the page at placement, so a
    // figure placed much smaller than its page needs a larger label size)
    label_size: f64,
    // draw highlight boxes with a translucent fill in addition to the outline
    highlight_fill: bool,
}

/// Human-readable sequence length for PDF labels.
fn fmt_len(bp: usize) -> String {
    if bp >= 1_000_000 {
        format!("{:.1} Mb", bp as f64 / 1_000_000.0)
    } else if bp >= 1_000 {
        format!("{} kb", (bp as f64 / 1_000.0).round() as usize)
    } else {
        format!("{} bp", bp)
    }
}

fn pdf_escape(s: &str) -> String {
    let mut out = String::new();
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '(' => out.push_str("\\("),
            ')' => out.push_str("\\)"),
            _ if (c as u32) < 32 || (c as u32) > 126 => out.push('?'),
            _ => out.push(c),
        }
    }
    out
}

/// Generate a PDF that embeds the dot plot raster with labeled rows (queries)
/// and columns (targets), plus mouseover tooltips on the label strips and on
/// every alignment-populated cell of the matrix.
fn generate_pdf_viewer(config: PdfViewerConfig) {
    use std::io::Write;
    let paf = config.paf;
    let (ax_w, ax_h) = (config.axes.0 as f64, config.axes.1 as f64);
    // scale the raster so the page remains a readable size (max ~850 pt)
    let img_scale = (850.0 / ax_w.max(ax_h)).min(1.0);
    let img_w = ax_w * img_scale;
    let img_h = ax_h * img_scale;
    let font_size = config.label_size;
    // margins sized to the longest labels so nothing clips
    let row_labels: Vec<String> = paf
        .queries
        .iter()
        .map(|q| format!("{} ({})", q.name, fmt_len(q.length)))
        .collect();
    let col_labels: Vec<String> = paf
        .targets
        .iter()
        .map(|t| format!("{} ({})", t.name, fmt_len(t.length)))
        .collect();
    let longest_q = row_labels.iter().map(|s| s.len()).max().unwrap_or(10) as f64;
    let longest_t = col_labels.iter().map(|s| s.len()).max().unwrap_or(10) as f64;
    let margin_l = (longest_q * font_size * 0.52 + 14.0).max(90.0);
    let margin_t = (longest_t * font_size * 0.52 + 14.0).max(90.0);
    let margin_b = if config.legend.is_some() { 40.0 } else { 24.0 };
    const MARGIN_R: f64 = 16.0;
    let page_w = margin_l + img_w + MARGIN_R;
    let page_h = margin_t + img_h + margin_b;
    let img_x = margin_l;
    let img_y = margin_b;

    // zlib-compress the raster for the image XObject
    let mut encoder = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
    encoder
        .write_all(config.raw)
        .expect("[pafplot::pdf] deflate raster");
    let img_data = encoder.finish().expect("[pafplot::pdf] finish deflate");

    // band geometry helpers: (start, end) pixel ranges for a sequence's span
    let tx = |off: usize, len: usize| -> (f64, f64) {
        if paf.target_length <= 0.0 {
            (0.0, 0.0)
        } else {
            (
                img_w * (off as f64 / paf.target_length),
                img_w * ((off + len) as f64 / paf.target_length),
            )
        }
    };
    let qy = |off: usize, len: usize| -> (f64, f64) {
        if paf.query_length <= 0.0 {
            (0.0, 0.0)
        } else {
            (
                img_h * (off as f64 / paf.query_length),
                img_h * ((off + len) as f64 / paf.query_length),
            )
        }
    };
    // band coords -> PDF user space. x: left to right. y: the raster stores
    // query offset 0 in its LAST row (bottom of the PNG), so query offset 0
    // maps to the BOTTOM of the image square, not the top.
    let px2pdf = |px: f64| img_x + px;
    let py2pdf = |py: f64| img_y + py;

    // visible text labels + tooltip annotations
    let mut label_ops = String::new();
    label_ops.push_str("BT\n/F1 7 Tf\n");
    struct Annot {
        rect: [f64; 4],
        tip: String,
    }
    let mut annots: Vec<Annot> = Vec::new();

    if config.using_zoom {
        // single-sequence zoom: one row label, one column label, one cell tooltip
        let t_mid = (config.ranges.0 .0 + config.ranges.0 .1) / 2;
        let q_mid = (config.ranges.1 .0 + config.ranges.1 .1) / 2;
        if let Some(t) = paf
            .targets
            .iter()
            .find(|s| t_mid >= s.offset && t_mid < s.offset + s.length)
        {
            label_ops.push_str(&format!(
                "0 1 -1 0 {:.2} {:.2} Tm ({}) Tj\n",
                img_x + 3.0,
                img_y + img_h + 3.0,
                pdf_escape(&t.name)
            ));
            annots.push(Annot {
                rect: [img_x, img_y + img_h + 3.0, img_x + img_w, page_h],
                tip: format!("target {} ({} bp)", t.name, t.length),
            });
        }
        if let Some(q) = paf
            .queries
            .iter()
            .find(|s| q_mid >= s.offset && q_mid < s.offset + s.length)
        {
            let qlabel = format!("{} ({})", q.name, fmt_len(q.length));
            label_ops.push_str(&format!(
                "1 0 0 1 {:.2} {:.2} Tm ({}) Tj\n",
                (margin_l - 8.0 - qlabel.len() as f64 * font_size * 0.52).max(4.0),
                img_y + img_h / 2.0 - 2.5,
                pdf_escape(&qlabel)
            ));
            annots.push(Annot {
                rect: [0.0, img_y, margin_l - 4.0, img_y + img_h],
                tip: format!("query {} ({} bp)", q.name, q.length),
            });
            if let Some(t) = paf
                .targets
                .iter()
                .find(|s| t_mid >= s.offset && t_mid < s.offset + s.length)
            {
                if let Some(stats) = config.pair_stats.get(&(q.name.clone(), t.name.clone())) {
                    annots.push(Annot {
                        rect: [img_x, img_y, img_x + img_w, img_y + img_h],
                        tip: format!(
                            "{} x {}: {} alignments, {} bp, best identity {:.1}%",
                            q.name,
                            t.name,
                            stats.0,
                            stats.1,
                            stats.2 * 100.0
                        ),
                    });
                }
            }
        }
    } else {
        // full matrix: label every row and column
        for (qi, q) in paf.queries.iter().enumerate() {
            let (py0, py1) = qy(q.offset, q.length);
            let y_mid = py2pdf((py0 + py1) / 2.0) - 2.5;
            let label = &row_labels[qi];
            label_ops.push_str(&format!(
                "1 0 0 1 {:.2} {:.2} Tm ({}) Tj\n",
                (margin_l - 8.0 - label.len() as f64 * font_size * 0.52).max(4.0),
                y_mid,
                pdf_escape(label)
            ));
            annots.push(Annot {
                rect: [
                    0.0,
                    py2pdf(py1).min(py2pdf(py0)),
                    margin_l - 4.0,
                    py2pdf(py1).max(py2pdf(py0)),
                ],
                tip: format!("query {} ({} bp)", q.name, q.length),
            });
        }
        for (ti, t) in paf.targets.iter().enumerate() {
            let (px0, px1) = tx(t.offset, t.length);
            let label = &col_labels[ti];
            label_ops.push_str(&format!(
                "0 1 -1 0 {:.2} {:.2} Tm ({}) Tj\n",
                px2pdf(px0) + 2.0,
                img_y + img_h + 3.0,
                pdf_escape(label)
            ));
            annots.push(Annot {
                rect: [px2pdf(px0), img_y + img_h + 3.0, px2pdf(px1), page_h],
                tip: format!("target {} ({} bp)", t.name, t.length),
            });
        }
        // cell tooltips for every alignment-populated pair
        for ((qname, tname), stats) in config.pair_stats.iter() {
            let (q, t) = match (find_seq(&paf.queries, qname), find_seq(&paf.targets, tname)) {
                (Some(q), Some(t)) => (q, t),
                _ => continue,
            };
            let (px0, px1) = tx(t.offset, t.length);
            let (py0, py1) = qy(q.offset, q.length);
            annots.push(Annot {
                rect: [
                    px2pdf(px0),
                    py2pdf(py1).min(py2pdf(py0)),
                    px2pdf(px1),
                    py2pdf(py1).max(py2pdf(py0)),
                ],
                tip: format!(
                    "{} x {}: {} alignments, {} bp, best identity {:.1}%",
                    qname,
                    tname,
                    stats.0,
                    stats.1,
                    stats.2 * 100.0
                ),
            });
        }
    }
    label_ops.push_str("ET\n");

    // highlight boxes for exchange candidates and other 2D regions
    let mut highlight_ops = String::new();
    for hl in config.highlights.iter() {
        // resolve the query and target; swap if the user gave them in reverse
        let (q, t) = match (
            find_seq(&paf.queries, &hl.qname),
            find_seq(&paf.targets, &hl.tname),
        ) {
            (Some(q), Some(t)) => (q, t),
            _ => match (
                find_seq(&paf.queries, &hl.tname),
                find_seq(&paf.targets, &hl.qname),
            ) {
                (Some(q), Some(t)) => (q, t),
                _ => {
                    eprintln!(
                        "[pafplot::pdf] highlight sequences not found, skipping: {} x {}",
                        hl.qname, hl.tname
                    );
                    continue;
                }
            },
        };
        let q_is_swapped = find_seq(&paf.queries, &hl.qname).is_none();
        let (qg0, qg1) = if q_is_swapped {
            (q.offset + hl.tstart, q.offset + hl.tend)
        } else {
            (q.offset + hl.qstart, q.offset + hl.qend)
        };
        let (tg0, tg1) = if q_is_swapped {
            (t.offset + hl.qstart, t.offset + hl.qend)
        } else {
            (t.offset + hl.tstart, t.offset + hl.tend)
        };
        // project the global 2D range into the image area
        let (x0, x1, y0, y1) = if config.using_zoom {
            let (tr, qr) = config.ranges;
            if tg1 <= tr.0 || tg0 >= tr.1 || qg1 <= qr.0 || qg0 >= qr.1 {
                continue; // outside the zoomed range
            }
            let tl = (tr.1 - tr.0) as f64;
            let ql = (qr.1 - qr.0) as f64;
            (
                img_w * ((tg0.max(tr.0) - tr.0) as f64 / tl),
                img_w * ((tg1.min(tr.1) - tr.0) as f64 / tl),
                img_h * ((qg0.max(qr.0) - qr.0) as f64 / ql),
                img_h * ((qg1.min(qr.1) - qr.0) as f64 / ql),
            )
        } else {
            (
                img_w * (tg0 as f64 / paf.target_length),
                img_w * (tg1 as f64 / paf.target_length),
                img_h * (qg0 as f64 / paf.query_length),
                img_h * (qg1 as f64 / paf.query_length),
            )
        };
        let (bx0, bx1) = (px2pdf(x0), px2pdf(x1));
        let (bya, byb) = (py2pdf(y1), py2pdf(y0));
        let (by0, by1) = if bya <= byb { (bya, byb) } else { (byb, bya) };
        let (cr, cg, cb) = hl.color;
        // expand sub-pixel boxes to a visible minimum, around their center,
        // so small tracts (e.g. a 3 kb patch in a 500 kb window) stay visible
        const MIN_BOX_PT: f64 = 6.0;
        let (mut bx0, mut bx1) = (bx0, bx1);
        let (mut by0, mut by1) = (by0, by1);
        if bx1 - bx0 < MIN_BOX_PT {
            let c = (bx0 + bx1) / 2.0;
            bx0 = c - MIN_BOX_PT / 2.0;
            bx1 = c + MIN_BOX_PT / 2.0;
        }
        if by1 - by0 < MIN_BOX_PT {
            let c = (by0 + by1) / 2.0;
            by0 = c - MIN_BOX_PT / 2.0;
            by1 = c + MIN_BOX_PT / 2.0;
        }
        // outline only by default: a filled box tints the dot plot and obscures
        // the identity coloring; --highlight-fill restores the translucent fill
        if config.highlight_fill {
            highlight_ops.push_str(&format!(
                "/GS0 gs\n{:.3} {:.3} {:.3} rg\n{:.2} {:.2} {:.2} {:.2} re f\n/GS1 gs\n{:.3} {:.3} {:.3} RG\n1 w\n{:.2} {:.2} {:.2} {:.2} re S\n",
                cr, cg, cb, bx0, by0, bx1 - bx0, by1 - by0, cr, cg, cb, bx0, by0, bx1 - bx0, by1 - by0
            ));
        } else {
            highlight_ops.push_str(&format!(
                "{:.3} {:.3} {:.3} RG\n1.5 w\n{:.2} {:.2} {:.2} {:.2} re S\n",
                cr,
                cg,
                cb,
                bx0,
                by0,
                bx1 - bx0,
                by1 - by0
            ));
        }
        // label: right of the box when there is room, otherwise left of it
        let label_text = hl.label.clone().unwrap_or_else(|| {
            format!(
                "{}:{}-{} x {}:{}-{}",
                hl.qname, hl.qstart, hl.qend, hl.tname, hl.tstart, hl.tend
            )
        });
        highlight_ops.push_str("BT\n/F1 6.5 Tf\n");
        let (tx_pos, ty_pos) = if bx1 + 150.0 < page_w {
            (bx1 + 3.0, by1 - 7.0)
        } else {
            (
                (bx0 - 3.0 - label_text.len() as f64 * 6.5 * 0.52).max(2.0),
                by1 - 7.0,
            )
        };
        highlight_ops.push_str(&format!(
            "{:.3} {:.3} {:.3} rg\n1 0 0 1 {:.2} {:.2} Tm ({}) Tj\nET\n",
            cr,
            cg,
            cb,
            tx_pos,
            ty_pos,
            pdf_escape(&label_text)
        ));
        annots.push(Annot {
            rect: [bx0, by0, bx1, by1],
            tip: label_text.clone(),
        });
    }

    // identity-color legend strip
    let mut legend_ops = String::new();
    if let Some((lmin, lmax, linv)) = config.legend {
        let strip_w = 200.0;
        let segs = 100;
        legend_ops.push_str(&format!(
            "BT\n/F1 6.5 Tf\n0 0 0 rg\n1 0 0 1 {:.2} {:.2} Tm (alignment identity) Tj\nET\n",
            img_x,
            margin_b - 6.0
        ));
        for s in 0..segs {
            let u = s as f64 / (segs - 1) as f64;
            let t = if linv { 1.0 - u } else { u };
            let c = turbo(t);
            legend_ops.push_str(&format!(
                "{:.3} {:.3} {:.3} rg\n{:.2} {:.2} {:.3} {:.2} re f\n",
                c.r as f64 / 255.0,
                c.g as f64 / 255.0,
                c.b as f64 / 255.0,
                img_x + u * strip_w,
                margin_b - 18.0,
                strip_w / segs as f64 + 0.2,
                8.0
            ));
        }
        legend_ops.push_str(&format!(
            "BT\n/F1 6.5 Tf\n0 0 0 rg\n1 0 0 1 {:.2} {:.2} Tm ({:.0}%) Tj\nET\n",
            img_x,
            margin_b - 30.0,
            lmin * 100.0
        ));
        legend_ops.push_str(&format!(
            "BT\n/F1 6.5 Tf\n0 0 0 rg\n1 0 0 1 {:.2} {:.2} Tm ({:.0}%) Tj\nET\n",
            img_x + strip_w / 2.0 - 12.0,
            margin_b - 30.0,
            (lmin + lmax) / 2.0 * 100.0
        ));
        legend_ops.push_str(&format!(
            "BT\n/F1 6.5 Tf\n0 0 0 rg\n1 0 0 1 {:.2} {:.2} Tm ({:.0}%) Tj\nET\n",
            img_x + strip_w - 20.0,
            margin_b - 30.0,
            lmax * 100.0
        ));
    }

    // assemble the page content stream
    let content = format!(
        "q\n{:.2} 0 0 {:.2} {:.2} {:.2} cm\n/Im0 Do\nQ\n{}{}{}",
        img_w, img_h, img_x, img_y, label_ops, highlight_ops, legend_ops
    );

    // ---- assemble the PDF file
    let n_annots = annots.len();
    let font_obj = 6;
    let gs_fill_obj = 7; // /GS0: translucent fill for highlight boxes
    let gs_full_obj = 8; // /GS1: opaque restore
    let first_annot = 9;
    let total_objs = if n_annots == 0 {
        gs_full_obj
    } else {
        first_annot + n_annots - 1
    }; // objects 1..=total_objs

    let mut buf: Vec<u8> = Vec::new();
    let mut offsets: Vec<usize> = Vec::new();
    let push_obj = |buf: &mut Vec<u8>, offsets: &mut Vec<usize>, body: Vec<u8>| -> usize {
        offsets.push(buf.len());
        let num = offsets.len();
        buf.extend_from_slice(format!("{} 0 obj\n", num).as_bytes());
        buf.extend_from_slice(&body);
        buf.extend_from_slice(b"\nendobj\n");
        num
    };

    buf.extend_from_slice(b"%PDF-1.7\n%\xE2\xE3\xCF\xD3\n");

    // 1: catalog
    push_obj(
        &mut buf,
        &mut offsets,
        b"<< /Type /Catalog /Pages 2 0 R >>".to_vec(),
    );
    // 2: pages
    push_obj(
        &mut buf,
        &mut offsets,
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>".to_vec(),
    );
    // 3: page
    let mut annot_refs = String::new();
    for i in first_annot..=total_objs {
        annot_refs.push_str(&format!("{} 0 R ", i));
    }
    let page_obj = format!(
        "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {:.2} {:.2}] /Resources << /XObject << /Im0 5 0 R >> /Font << /F1 {} 0 R >> /ExtGState << /GS0 {} 0 R /GS1 {} 0 R >> >> /Contents 4 0 R /Annots [{}]>>",
        page_w, page_h, font_obj, gs_fill_obj, gs_full_obj, annot_refs
    );
    push_obj(&mut buf, &mut offsets, page_obj.into_bytes());
    // 4: content stream
    let content_obj = format!(
        "<< /Length {} >>\nstream\n{}\nendstream",
        content.len(),
        content
    );
    push_obj(&mut buf, &mut offsets, content_obj.into_bytes());
    // 5: image XObject (FlateDecode, PNG predictor)
    let image_obj = format!(
        "<< /Type /XObject /Subtype /Image /Width {} /Height {} /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /FlateDecode /Length {} >>\nstream\n",
        config.axes.0, config.axes.1, img_data.len()
    );
    let mut image_body = image_obj.into_bytes();
    image_body.extend_from_slice(&img_data);
    image_body.extend_from_slice(b"\nendstream");
    push_obj(&mut buf, &mut offsets, image_body);
    // 6: font
    push_obj(
        &mut buf,
        &mut offsets,
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>".to_vec(),
    );
    // 7: translucent fill state; 8: opaque state
    push_obj(
        &mut buf,
        &mut offsets,
        b"<< /Type /ExtGState /ca 0.18 /CA 0.18 >>".to_vec(),
    );
    push_obj(
        &mut buf,
        &mut offsets,
        b"<< /Type /ExtGState /ca 1 /CA 1 >>".to_vec(),
    );
    // 9..: link annotations carrying hover tooltips
    for a in annots.iter() {
        let annot_obj = format!(
            "<< /Type /Annot /Subtype /Link /Rect [{:.2} {:.2} {:.2} {:.2}] /Border [0 0 0] /F 4 /TU ({}) >>",
            a.rect[0], a.rect[1], a.rect[2], a.rect[3],
            pdf_escape(&a.tip)
        );
        push_obj(&mut buf, &mut offsets, annot_obj.into_bytes());
    }

    // xref + trailer
    let xref_off = buf.len();
    buf.extend_from_slice(format!("xref\n0 {}\n", total_objs + 1).as_bytes());
    buf.extend_from_slice(b"0000000000 65535 f \n");
    for off in offsets.iter() {
        buf.extend_from_slice(format!("{:010} 00000 n \n", off).as_bytes());
    }
    buf.extend_from_slice(
        format!(
            "trailer << /Size {} /Root 1 0 R >>\nstartxref\n{}\n%%EOF\n",
            total_objs + 1,
            xref_off
        )
        .as_bytes(),
    );

    let path = &Path::new(config.output_filename);
    if let Err(e) = std::fs::write(path, &buf) {
        panic!("[pafplot::pdf] failed to write pdf: {:?}", e);
    }
    eprintln!(
        "[pafplot] PDF viewer written to: {} ({:.0}x{:.0} pt, {} row labels, {} column labels, {} tooltips)",
        config.output_filename,
        page_w,
        page_h,
        paf.queries.len(),
        paf.targets.len(),
        n_annots
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_paf_parsing_functions() {
        let test_line =
            "query1\t1000\t100\t900\t+\ttarget1\t2000\t200\t1800\t800\t800\t60\tcg:Z:800M";

        assert_eq!(paf_query(test_line), "query1");
        assert_eq!(paf_query_length(test_line), 1000);
        assert_eq!(paf_query_begin(test_line), 100);
        assert_eq!(paf_query_end(test_line), 900);
        assert!(!paf_query_is_rev(test_line));

        assert_eq!(paf_target(test_line), "target1");
        assert_eq!(paf_target_length(test_line), 2000);
        assert_eq!(paf_target_begin(test_line), 200);
        assert_eq!(paf_target_end(test_line), 1800);
    }

    #[test]
    fn test_paf_reverse_strand() {
        let test_line =
            "query1\t1000\t100\t900\t-\ttarget1\t2000\t200\t1800\t800\t800\t60\tcg:Z:800M";
        assert!(paf_query_is_rev(test_line));
    }

    #[test]
    fn test_aligned_seq_creation() {
        let seq = AlignedSeq::new();
        assert_eq!(seq.name, "");
        assert_eq!(seq.length, 0);
        assert_eq!(seq.rank, 0);
        assert_eq!(seq.offset, 0);
    }

    #[test]
    fn test_coordinate_projection() {
        // Create a simple test PAF with known dimensions
        std::fs::write(
            "test_simple.paf",
            "query1\t1000\t0\t1000\t+\ttarget1\t2000\t0\t2000\t1000\t1000\t60\tcg:Z:1000M\n",
        )
        .expect("Failed to write test file");

        let paf = PafFile::new("test_simple.paf");
        let axes = (100, 50); // 100x50 canvas

        // Test projection of coordinates
        let (proj_x, proj_y) = paf.project_xy(1000, 500, axes);
        assert!((proj_x - 50.0).abs() < 0.1); // Should be roughly half of target axis
        assert!((proj_y - 25.0).abs() < 0.1); // Should be roughly half of query axis

        // Clean up
        std::fs::remove_file("test_simple.paf").ok();
    }

    #[test]
    fn test_axis_calculation() {
        std::fs::write(
            "test_axis.paf",
            "query1\t1000\t0\t1000\t+\ttarget1\t2000\t0\t2000\t1000\t1000\t60\tcg:Z:1000M\n",
        )
        .expect("Failed to write test file");

        let paf = PafFile::new("test_axis.paf");
        let axes = paf.get_axes(1000);

        // Target is longer than query, so target should get the major axis
        assert_eq!(axes.0, 1000); // target axis gets major axis
        assert_eq!(axes.1, 500); // query axis is scaled down proportionally

        std::fs::remove_file("test_axis.paf").ok();
    }

    #[test]
    fn test_error_handling_invalid_paf() {
        // Test handling of malformed PAF lines
        let invalid_line = "incomplete\tline";

        // These should not panic but may return defaults or handle gracefully
        // The actual behavior depends on implementation details
        let result = std::panic::catch_unwind(|| paf_query_length(invalid_line));
        assert!(result.is_err()); // Should panic on invalid input
    }
}
