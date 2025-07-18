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
    let reader = BufReader::new(file);
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
        F: FnMut(char, usize, bool, usize, usize),
    {
        let query_rev = paf_query_is_rev(line);
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
        let fake_cigar = format!("{}M", if query_len < target_len { query_len } else { target_len });
        if cigars.is_empty() {
            cigars.push(fake_cigar.as_str());
        }

        for cigar in cigars
        {
            //println!("{}", cigar);
            let mut first: usize = 0;
            for (i, b) in cigar.bytes().enumerate() {
                let c = b as char;
                //println!("{} {}", i, b as char);
                match c {
                    'M' | '=' | 'X' => {
                        let n = cigar[first..i].parse::<usize>().unwrap() as usize;
                        func(c, query_pos, query_rev, target_pos, n);
                        query_pos += if query_rev { 0 - n } else { n };
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
                        query_pos += if query_rev { 0 - n } else { n };
                        first = i + 1;
                    }
                    _ => {}
                }
            }
        }
    }
    fn for_each_match_in_file<F>(self: &PafFile, mut func: F)
    where
        F: FnMut(char, usize, bool, usize, usize),
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
            self.for_each_match(line, |c, x, r, y, d| func(c, x, r, y, d));
            //println!();
        });
    }
    fn get_axes(self: &PafFile, major_axis: usize) -> (usize, usize) {
        //let max_length = cmp::max(self.query_length, self.target_length) as f64;
        //println!("query = {} target = {}", self.query_length, self.target_length);
        if self.target_length > self.query_length {
            let ratio = self.query_length as f64 / self.target_length as f64;
            (major_axis, (major_axis as f64 * ratio) as usize)
        } else {
            let ratio = self.target_length as f64 / self.query_length as f64;
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
            axes.0 as f64 * (x as f64 / self.target_length as f64),
            axes.1 as f64 * (y as f64 / self.query_length as f64),
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
            axes.0 as f64
                * (((x as f64) - (zoom.0 .0 as f64)) / ((zoom.0 .1 - zoom.0 .0) as f64)),
            axes.1 as f64
                * (((y as f64) - (zoom.1 .0 as f64)) / ((zoom.1 .1 - zoom.1 .0) as f64)),
        )
    }
    /*
        fn get_pixel(self: &PafFile, x: i64, y: i64, axes: (usize, usize)) -> usize {
            let (q, t) = axes;
            x * (self.query_length / q as f64).round() as usize * q + y * (self.target_length / t as f64).round() as usize
        }
    */
}

fn main() {
    let matches = App::new("pafplot")
        .version("0.1.0")
        .author("Erik Garrison <erik.garrison@gmail.com>")
        .about("Generate a dotplot from pairwise DNA alignments in PAF format")
        .arg(
            Arg::with_name("INPUT")
                .required(true)
                .takes_value(true)
                .index(1)
                .help("input PAF file"),
        )
        .arg(
            Arg::with_name("png")
                .takes_value(true)
                .short("p")
                .long("png")
                .help("Save the dotplot to this file."),
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
                .help("Generate an interactive HTML viewer with sequence labels and coordinates"),
        )
        .get_matches();

    let filename = matches.value_of("INPUT").unwrap();
    let paf = PafFile::new(filename);

    let major_axis = matches
        .value_of("size")
        .unwrap_or(&"1000")
        .parse::<usize>()
        .unwrap();

    let default_output = format!("{}.png", filename);

    let output_png = matches.value_of("png").unwrap_or(&default_output);

    let dark = matches.is_present("dark");

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

        (paf.target_range(target_name, target_start, target_end),
         paf.query_range(query_name, query_start, query_end))
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
    //println!("axes = {} {}", axes.0, axes.1);
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
    let draw_match = |_c, x: usize, rev: bool, y: usize, len: usize| {
        //println!("draw_match {} {} {} {} {}", _c, x, rev, y, len);
        let start = get_coords(x, y);
        let end = get_coords(x + if rev { 0 - len } else { len }, y + len);
        /*
        println!(
            "start and end ({} {}) ({} {})",
            start.0, start.1, end.0, end.1
        );
         */
        for ((j, i), val) in XiaolinWu::<f64, i64>::new(start, end) {
            //println!("checking pixel {} {} {}", i, j, val);
            if i >= 0 && i < (axes.0 as i64) && j >= 0 && j < (axes.1 as i64) {
                //println!("drawing pixel {} {} {}", i, j, val);
                let i: usize = (i as usize) + (((axes.1 - 1) - j as usize) * axes.0);
                pixels[i] = get_color(val);
            }
        }
    };
    paf.for_each_match_in_file(draw_match);

    let path = &Path::new(output_png);

    // encode_file takes the path to the image, a u8 array,
    // the width, the height, the color mode, and the bit depth
    if let Err(e) = lodepng::encode_file(path, &raw, axes.0, axes.1, lodepng::ColorType::RGB, 8) {
        panic!("failed to write png: {:?}", e);
    }

    // Generate HTML viewer if requested
    if matches.is_present("html") {
        generate_html_viewer(&paf, axes, output_png, dark, using_zoom, (target_range, query_range));
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
        let align_len = std::cmp::min(query_len, target_len);
        
        alignments.push(format!(
            r#"{{"x":{},"y":{},"len":{},"rev":{}}}"#,
            x, y, align_len, query_rev
        ));
        
        // Update summary grid
        let grid_x = (x as f64 / paf.target_length * grid_size as f64) as usize;
        let grid_y = (y as f64 / paf.query_length * grid_size as f64) as usize;
        let entry = summary_grid.entry((grid_x, grid_y)).or_insert((0, 0, 0));
        entry.0 += 1;
        entry.1 += align_len;
        if query_rev { entry.2 += 1; }
        
        // Collect detailed CIGAR operations for base-level view
        let mut cigar_ops = Vec::new();
        
        paf.for_each_match(line, |c, target_pos, rev, query_pos, len| {
            let op_type = match c {
                'M' | '=' => "match",
                'X' => "mismatch",
                'I' => "insertion",
                'D' => "deletion",
                _ => "unknown"
            };
            
            cigar_ops.push(format!(
                r#"{{"type":"{}","x":{},"y":{},"len":{},"rev":{}}}"#,
                op_type, target_pos, query_pos, len, rev
            ));
        });
        
        if !cigar_ops.is_empty() {
            detailed_alignments.push(format!(
                r#"{{"x":{},"y":{},"rev":{},"ops":[{}]}}"#,
                x, y, query_rev, cigar_ops.join(",")
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
            r#"{{"gx":{},"gy":{},"density":{},"avgLen":{},"revRatio":{}}}"#,
            gx, gy, density, avg_len, rev_ratio
        ));
    }
    
    let alignments_json = format!("[{}]", alignments.join(","));
    let summary_json = format!("[{}]", summary_data.join(","));
    let detailed_json = format!("[{}]", detailed_alignments.join(","));
    
    (alignments_json, summary_json, detailed_json)
}

fn generate_html_viewer(
    paf: &PafFile,
    _axes: (usize, usize),
    output_filename: &str,
    dark: bool,
    using_zoom: bool,
    ranges: ((usize, usize), (usize, usize)),
) {
    // Collect alignment data for canvas rendering
    let (alignments_json, summary_json, detailed_json) = collect_alignment_data(paf);

    // Generate sequence metadata as JSON
    let mut targets_json = String::from("[");
    for (i, target) in paf.targets.iter().enumerate() {
        if i > 0 { targets_json.push(','); }
        targets_json.push_str(&format!(
            r#"{{"name":"{}","length":{},"rank":{},"offset":{}}}"#,
            target.name.replace('"', r#"\""#), target.length, target.rank, target.offset
        ));
    }
    targets_json.push(']');

    let mut queries_json = String::from("[");
    for (i, query) in paf.queries.iter().enumerate() {
        if i > 0 { queries_json.push(','); }
        queries_json.push_str(&format!(
            r#"{{"name":"{}","length":{},"rank":{},"offset":{}}}"#,
            query.name.replace('"', r#"\""#), query.length, query.rank, query.offset
        ));
    }
    queries_json.push(']');

    let html_content = format!(r#"<!DOCTYPE html>
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

        function prepareLineData() {{
            lineVertices = [];
            lineColors = [];
            
            // Calculate visible bounds
            const viewLeft = -panX / zoom;
            const viewRight = (canvasWidth - panX) / zoom;
            const viewTop = -panY / zoom;
            const viewBottom = (canvasHeight - panY) / zoom;
            const viewWidth = viewRight - viewLeft;
            const viewHeight = viewBottom - viewTop;
            
            // Calculate viewport width in base pairs
            const viewportWidthBp = Math.max(
                usingZoom ? 
                    ((viewWidth / plotWidth) * (targetRange[1] - targetRange[0])) :
                    ((viewWidth / plotWidth) * targetLength),
                usingZoom ? 
                    ((viewHeight / plotHeight) * (queryRange[1] - queryRange[0])) :
                    ((viewHeight / plotHeight) * queryLength)
            );
            
            // Performance optimization: switch rendering modes based on viewport size
            const VIEWPORT_THRESHOLD = 50000; // 50kbp threshold for performance
            const shouldShowAlignmentSegments = viewportWidthBp > VIEWPORT_THRESHOLD;
            const shouldShowDetails = showDetails && !shouldShowAlignmentSegments && viewportWidthBp < 1000; // Show CIGAR details only when very zoomed in
            
            // Add grid lines (sequence boundaries)
            const gridColor = hexToRgb(borderColor);
            const endLineColor = hexToRgb(darkMode ? '#666666' : '#999999');
            
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
                    lineColors.push(gridColor.r, gridColor.g, gridColor.b, gridOpacity);
                    lineColors.push(gridColor.r, gridColor.g, gridColor.b, gridOpacity);
                }}
                
                // Add sequence end lines
                const endX = target.offset + target.length;
                const endStart = projectCoords(endX, 0);
                const endEnd = projectCoords(endX, queryLength);
                lineVertices.push(endStart.x, endStart.y);
                lineVertices.push(endEnd.x, endEnd.y);
                lineColors.push(endLineColor.r, endLineColor.g, endLineColor.b, gridOpacity * 0.7);
                lineColors.push(endLineColor.r, endLineColor.g, endLineColor.b, gridOpacity * 0.7);
            }});
            
            // Query boundaries (horizontal lines)
            queries.forEach(query => {{
                if (query.offset > 0) {{
                    const coords = projectCoords(0, query.offset);
                    const endCoords = projectCoords(targetLength, query.offset);
                    
                    lineVertices.push(coords.x, coords.y);
                    lineVertices.push(endCoords.x, endCoords.y);
                    lineColors.push(gridColor.r, gridColor.g, gridColor.b, gridOpacity);
                    lineColors.push(gridColor.r, gridColor.g, gridColor.b, gridOpacity);
                }}
                
                // Add sequence end lines
                const endY = query.offset + query.length;
                const endStart = projectCoords(0, endY);
                const endEnd = projectCoords(targetLength, endY);
                lineVertices.push(endStart.x, endStart.y);
                lineVertices.push(endEnd.x, endEnd.y);
                lineColors.push(endLineColor.r, endLineColor.g, endLineColor.b, gridOpacity * 0.7);
                lineColors.push(endLineColor.r, endLineColor.g, endLineColor.b, gridOpacity * 0.7);
            }});
            
            // Draw alignments - three rendering modes based on viewport size
            const matchColor = hexToRgb(lineColor);
            const mismatchColor = hexToRgb(darkMode ? '#ff4444' : '#cc0000');
            const indelColor = hexToRgb(darkMode ? '#4444ff' : '#0000cc');
            const summaryColor = hexToRgb(darkMode ? '#888888' : '#666666');
            
            if (shouldShowAlignmentSegments) {{
                // Mode 1: Alignment segments view (viewport > 50kbp)
                // Draw simplified segments with reduced opacity for performance
                let renderCount = 0;
                const maxRenderCount = 10000; // Limit rendering for extreme performance
                
                alignments.forEach(alignment => {{
                    if (renderCount >= maxRenderCount) return;
                    
                    const startCoords = projectCoords(alignment.x, alignment.y);
                    const endX = alignment.x + (alignment.rev ? -alignment.len : alignment.len);
                    const endY = alignment.y + alignment.len;
                    const endCoords = projectCoords(endX, endY);
                    
                    // Enhanced frustum culling for performance
                    const minX = Math.min(startCoords.x, endCoords.x);
                    const maxX = Math.max(startCoords.x, endCoords.x);
                    const minY = Math.min(startCoords.y, endCoords.y);
                    const maxY = Math.max(startCoords.y, endCoords.y);
                    
                    if (maxX >= viewLeft && minX <= viewRight && maxY >= viewTop && minY <= viewBottom) {{
                        // Use lower opacity for performance and visual clarity at this scale
                        const alpha = Math.min(0.6, 1.0 / Math.sqrt(viewportWidthBp / 50000));
                        lineVertices.push(startCoords.x, startCoords.y);
                        lineVertices.push(endCoords.x, endCoords.y);
                        lineColors.push(summaryColor.r, summaryColor.g, summaryColor.b, alpha);
                        lineColors.push(summaryColor.r, summaryColor.g, summaryColor.b, alpha);
                        renderCount++;
                    }}
                }});
            }} else if (shouldShowDetails && detailedAlignments.length > 0) {{
                // Mode 2: Detailed CIGAR operations (viewport < 1kbp and details enabled)
                detailedAlignments.forEach(alignment => {{
                    alignment.ops.forEach(op => {{
                        const startCoords = projectCoords(op.x, op.y);
                        
                        if (op.type === 'match' || op.type === 'unknown') {{
                            const endX = op.x + op.len;
                            const endY = op.y + (op.rev ? -op.len : op.len);
                            const endCoords = projectCoords(endX, endY);
                            
                            lineVertices.push(startCoords.x, startCoords.y);
                            lineVertices.push(endCoords.x, endCoords.y);
                            lineColors.push(matchColor.r, matchColor.g, matchColor.b, 0.9);
                            lineColors.push(matchColor.r, matchColor.g, matchColor.b, 0.9);
                        }} else if (op.type === 'mismatch') {{
                            // Draw mismatches as individual segments
                            for (let i = 0; i < Math.min(op.len, 100); i++) {{ // Limit to prevent too many segments
                                const mx = op.x + i;
                                const my = op.y + (op.rev ? -i : i);
                                const mStartCoords = projectCoords(mx, my);
                                const mEndCoords = projectCoords(mx + 1, my + (op.rev ? -1 : 1));
                                
                                lineVertices.push(mStartCoords.x, mStartCoords.y);
                                lineVertices.push(mEndCoords.x, mEndCoords.y);
                                lineColors.push(mismatchColor.r, mismatchColor.g, mismatchColor.b, 1.0);
                                lineColors.push(mismatchColor.r, mismatchColor.g, mismatchColor.b, 1.0);
                            }}
                        }} else if (op.type === 'insertion') {{
                            const endCoords = projectCoords(op.x, op.y + (op.rev ? -op.len : op.len));
                            
                            lineVertices.push(startCoords.x, startCoords.y);
                            lineVertices.push(endCoords.x, endCoords.y);
                            lineColors.push(indelColor.r, indelColor.g, indelColor.b, 0.8);
                            lineColors.push(indelColor.r, indelColor.g, indelColor.b, 0.8);
                        }} else if (op.type === 'deletion') {{
                            const endCoords = projectCoords(op.x + op.len, op.y);
                            
                            lineVertices.push(startCoords.x, startCoords.y);
                            lineVertices.push(endCoords.x, endCoords.y);
                            lineColors.push(indelColor.r, indelColor.g, indelColor.b, 0.8);
                            lineColors.push(indelColor.r, indelColor.g, indelColor.b, 0.8);
                        }}
                    }});
                }});
            }} else {{
                // Mode 3: Standard alignment lines (1kbp - 50kbp viewport)
                alignments.forEach(alignment => {{
                    const startCoords = projectCoords(alignment.x, alignment.y);
                    const endX = alignment.x + (alignment.rev ? -alignment.len : alignment.len);
                    const endY = alignment.y + alignment.len;
                    const endCoords = projectCoords(endX, endY);
                    
                    // Simple frustum culling
                    const minX = Math.min(startCoords.x, endCoords.x);
                    const maxX = Math.max(startCoords.x, endCoords.x);
                    const minY = Math.min(startCoords.y, endCoords.y);
                    const maxY = Math.max(startCoords.y, endCoords.y);
                    
                    if (maxX >= viewLeft && minX <= viewRight && maxY >= viewTop && minY <= viewBottom) {{
                        lineVertices.push(startCoords.x, startCoords.y);
                        lineVertices.push(endCoords.x, endCoords.y);
                        lineColors.push(matchColor.r, matchColor.g, matchColor.b, 0.8);
                        lineColors.push(matchColor.r, matchColor.g, matchColor.b, 0.8);
                    }}
                }});
            }}
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
                    const endX = alignment.x + (alignment.rev ? -alignment.len : alignment.len);
                    const endY = alignment.y + alignment.len;
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
            
            // Calculate viewport info for display
            const viewLeft = -panX / zoom;
            const viewRight = (canvasWidth - panX) / zoom;
            const viewTop = -panY / zoom;
            const viewBottom = (canvasHeight - panY) / zoom;
            const viewWidth = viewRight - viewLeft;
            const viewHeight = viewBottom - viewTop;
            
            const viewportWidthBp = Math.max(
                usingZoom ? 
                    ((viewWidth / plotWidth) * (targetRange[1] - targetRange[0])) :
                    ((viewWidth / plotWidth) * targetLength),
                usingZoom ? 
                    ((viewHeight / plotHeight) * (queryRange[1] - queryRange[0])) :
                    ((viewHeight / plotHeight) * queryLength)
            );
            
            // Update rendering mode indicator
            let mode = 'Standard';
            if (viewportWidthBp > 50000) {{
                mode = 'Segments';
            }} else if (showDetails && viewportWidthBp < 1000) {{
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
        output_filename,
        if dark { "#1a1a1a" } else { "#ffffff" },
        if dark { "#ffffff" } else { "#000000" },
        if dark { "#444444" } else { "#cccccc" },
        if dark { "#2a2a2a" } else { "#f5f5f5" },
        if dark { "#444444" } else { "#cccccc" },
        if dark { "#2a2a2a" } else { "#f5f5f5" },  // config panel background
        if dark { "#444444" } else { "#cccccc" },  // config panel border
        if dark { "#2a2a2a" } else { "#f0f0f0" },  // minimap background
        if dark { "#666666" } else { "#999999" },  // minimap border
        if dark { "#00ff00" } else { "#0088ff" },  // viewport border
        if dark { "rgba(0,255,0,0.2)" } else { "rgba(0,136,255,0.2)" },  // viewport fill
        if dark { "#666666" } else { "#999999" },
        if dark { "rgba(26,26,26,0.8)" } else { "rgba(255,255,255,0.8)" },  // label background
        if dark { "#1a1a1a" } else { "#ffffff" },  // page background again
        if dark { "#00ff00" } else { "#0088ff" },  // render mode color
        alignments_json,
        summary_json,
        detailed_json,
        targets_json,
        queries_json,
        if using_zoom { ranges.0.1 - ranges.0.0 } else { paf.target_length as usize },
        if using_zoom { ranges.1.1 - ranges.1.0 } else { paf.query_length as usize },
        using_zoom,
        format!("[{}, {}]", ranges.0.0, ranges.0.1),
        format!("[{}, {}]", ranges.1.0, ranges.1.1),
        dark,
    );

    let html_filename = if output_filename.ends_with(".png") {
        output_filename.replace(".png", ".html")
    } else {
        format!("{}.html", output_filename)
    };
    std::fs::write(&html_filename, html_content)
        .expect("Failed to write HTML file");
    
    println!("Generated interactive HTML viewer: {}", html_filename);
}