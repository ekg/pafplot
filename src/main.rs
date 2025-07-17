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
        generate_html_viewer(&paf, &raw, axes, output_png, dark, using_zoom, (target_range, query_range));
    }
}

fn collect_alignment_data(paf: &PafFile) -> String {
    let mut alignments = Vec::new();
    
    paf.for_each_match_in_file(|_c, x, rev, y, len| {
        alignments.push(format!(
            r#"{{"x":{},"y":{},"len":{},"rev":{}}}"#,
            x, y, len, rev
        ));
    });
    
    format!("[{}]", alignments.join(","))
}

fn generate_html_viewer(
    paf: &PafFile,
    _raw_image: &[u8],
    axes: (usize, usize),
    png_filename: &str,
    dark: bool,
    using_zoom: bool,
    ranges: ((usize, usize), (usize, usize)),
) {
    // Collect alignment data for canvas rendering
    let alignments_json = collect_alignment_data(paf);

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
            padding: 20px;
        }}
        .container {{
            max-width: 100vw;
            overflow: auto;
        }}
        .plot-area {{
            position: relative;
            display: inline-block;
            margin: 50px;
        }}
        #plotCanvas {{
            display: block;
            border: 1px solid {};
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
            position: absolute;
            font-size: 9px;
            font-family: 'Monaco', 'Menlo', 'Ubuntu Mono', monospace;
            pointer-events: none;
        }}
        .target-labels {{
            bottom: -30px;
            left: 0;
            width: 100%;
            height: 25px;
        }}
        .query-labels {{
            right: -120px;
            top: 0;
            width: 115px;
            height: 100%;
        }}
        .label {{
            position: absolute;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
        }}
        .target-label {{
            transform-origin: left bottom;
            transform: rotate(-45deg);
            bottom: 0;
        }}
        .query-label {{
            text-align: right;
            right: 5px;
        }}
        .info-panel {{
            position: fixed;
            top: 20px;
            right: 20px;
            background-color: {};
            border: 1px solid {};
            padding: 15px;
            border-radius: 5px;
            min-width: 300px;
            font-size: 12px;
            z-index: 20;
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
    <div class="container">
        <h1>PAF Plot Viewer</h1>
        <p>File: {} | Alignments: {} | Target sequences: {} | Query sequences: {}</p>
        
        <div class="controls">
            <button onclick="resetZoom()">Reset Zoom</button>
            <span class="zoom-info">Scroll wheel to zoom, drag to pan</span>
        </div>
        
        <div class="plot-area" id="plotArea">
            <canvas id="plotCanvas" width="{}" height="{}"></canvas>
            <div class="tooltip" id="tooltip"></div>
            
            <div class="labels target-labels" id="targetLabels"></div>
            <div class="labels query-labels" id="queryLabels"></div>
        </div>
        
        <div class="info-panel">
            <h3>Current Position</h3>
            <div id="position">Move mouse over plot</div>
            <h3>Sequence Info</h3>
            <div id="sequenceInfo">Hover over plot for details</div>
            <div class="zoom-info">
                Zoom: <span id="zoomLevel">1.0x</span><br>
                Pan: <span id="panInfo">0, 0</span>
            </div>
        </div>
    </div>

    <script>
        const canvas = document.getElementById('plotCanvas');
        const ctx = canvas.getContext('2d');
        const tooltip = document.getElementById('tooltip');
        const position = document.getElementById('position');
        const sequenceInfo = document.getElementById('sequenceInfo');
        const zoomLevel = document.getElementById('zoomLevel');
        const panInfo = document.getElementById('panInfo');
        const targetLabels = document.getElementById('targetLabels');
        const queryLabels = document.getElementById('queryLabels');

        // Plot data
        const alignments = {};
        const targets = {};
        const queries = {};
        const canvasWidth = {};
        const canvasHeight = {};
        const targetLength = {};
        const queryLength = {};
        const usingZoom = {};
        const targetRange = {};
        const queryRange = {};
        const darkMode = {};

        // Zoom and pan state
        let zoom = 1.0;
        let panX = 0;
        let panY = 0;
        let isDragging = false;
        let dragStartX = 0;
        let dragStartY = 0;
        let dragStartPanX = 0;
        let dragStartPanY = 0;

        // Colors
        const backgroundColor = darkMode ? '#1a1a1a' : '#ffffff';
        const lineColor = darkMode ? '#ffffff' : '#000000';
        const borderColor = darkMode ? '#444444' : '#cccccc';

        function drawPlot() {{
            ctx.clearRect(0, 0, canvasWidth, canvasHeight);
            
            // Set background
            ctx.fillStyle = backgroundColor;
            ctx.fillRect(0, 0, canvasWidth, canvasHeight);
            
            // Draw sequence boundaries
            ctx.strokeStyle = borderColor;
            ctx.lineWidth = 0.5;
            
            // Target boundaries (vertical lines)
            targets.forEach(target => {{
                const x = ((target.offset / targetLength) * canvasWidth) * zoom + panX;
                if (x >= 0 && x <= canvasWidth) {{
                    ctx.beginPath();
                    ctx.moveTo(x, 0);
                    ctx.lineTo(x, canvasHeight);
                    ctx.stroke();
                }}
            }});
            
            // Query boundaries (horizontal lines)
            queries.forEach(query => {{
                const y = (((queryLength - query.offset - query.length) / queryLength) * canvasHeight) * zoom + panY;
                if (y >= 0 && y <= canvasHeight) {{
                    ctx.beginPath();
                    ctx.moveTo(0, y);
                    ctx.lineTo(canvasWidth, y);
                    ctx.stroke();
                }}
            }});
            
            // Draw alignments
            ctx.strokeStyle = lineColor;
            ctx.lineWidth = 1.0;
            ctx.globalAlpha = 0.7;
            
            alignments.forEach(alignment => {{
                const startX = ((alignment.x / targetLength) * canvasWidth) * zoom + panX;
                const startY = (((queryLength - alignment.y) / queryLength) * canvasHeight) * zoom + panY;
                const endX = (((alignment.x + (alignment.rev ? -alignment.len : alignment.len)) / targetLength) * canvasWidth) * zoom + panX;
                const endY = (((queryLength - alignment.y - alignment.len) / queryLength) * canvasHeight) * zoom + panY;
                
                // Only draw if visible
                if ((startX >= -10 && startX <= canvasWidth + 10) || (endX >= -10 && endX <= canvasWidth + 10)) {{
                    ctx.beginPath();
                    ctx.moveTo(startX, startY);
                    ctx.lineTo(endX, endY);
                    ctx.stroke();
                }}
            }});
            
            ctx.globalAlpha = 1.0;
        }}

        function pixelToGenomicCoords(pixelX, pixelY) {{
            const genomicX = ((pixelX - panX) / zoom) * (targetLength / canvasWidth);
            const genomicY = queryLength - ((pixelY - panY) / zoom) * (queryLength / canvasHeight);
            
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

        function updateLabels() {{
            // Clear existing labels
            targetLabels.innerHTML = '';
            queryLabels.innerHTML = '';
            
            // Target labels
            targets.forEach(target => {{
                const x = ((target.offset / targetLength) * canvasWidth) * zoom + panX;
                if (x >= -50 && x <= canvasWidth + 50) {{
                    const label = document.createElement('div');
                    label.className = 'label target-label';
                    label.textContent = target.name;
                    label.title = `${{target.name}} (${{target.length.toLocaleString()}} bp)`;
                    label.style.left = (x + 5) + 'px';
                    targetLabels.appendChild(label);
                }}
            }});
            
            // Query labels
            queries.forEach(query => {{
                const y = (((queryLength - query.offset - query.length) / queryLength) * canvasHeight) * zoom + panY;
                if (y >= -20 && y <= canvasHeight + 20) {{
                    const label = document.createElement('div');
                    label.className = 'label query-label';
                    label.textContent = query.name;
                    label.title = `${{query.name}} (${{query.length.toLocaleString()}} bp)`;
                    label.style.top = y + 'px';
                    queryLabels.appendChild(label);
                }}
            }});
        }}

        function updateDisplay() {{
            drawPlot();
            updateLabels();
            zoomLevel.textContent = zoom.toFixed(2) + 'x';
            panInfo.textContent = `${{Math.round(panX)}}, ${{Math.round(panY)}}`;
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
            isDragging = true;
            dragStartX = e.clientX - canvas.getBoundingClientRect().left;
            dragStartY = e.clientY - canvas.getBoundingClientRect().top;
            dragStartPanX = panX;
            dragStartPanY = panY;
            canvas.style.cursor = 'grabbing';
        }});

        canvas.addEventListener('mouseup', function() {{
            isDragging = false;
            canvas.style.cursor = 'crosshair';
        }});

        canvas.addEventListener('wheel', function(e) {{
            e.preventDefault();
            const rect = canvas.getBoundingClientRect();
            const mouseX = e.clientX - rect.left;
            const mouseY = e.clientY - rect.top;
            
            const zoomFactor = e.deltaY > 0 ? 0.9 : 1.1;
            const newZoom = Math.max(0.1, Math.min(10, zoom * zoomFactor));
            
            if (newZoom !== zoom) {{
                // Zoom towards mouse position
                const genomicBefore = pixelToGenomicCoords(mouseX, mouseY);
                zoom = newZoom;
                const genomicAfter = pixelToGenomicCoords(mouseX, mouseY);
                
                panX += (genomicAfter.target - genomicBefore.target) * (canvasWidth / targetLength) * zoom;
                panY += (genomicAfter.query - genomicBefore.query) * (canvasHeight / queryLength) * zoom;
                
                updateDisplay();
            }}
        }});

        function resetZoom() {{
            zoom = 1.0;
            panX = 0;
            panY = 0;
            updateDisplay();
        }}

        // Initialize
        updateDisplay();
    </script>
</body>
</html>"#,
        png_filename,
        if dark { "#1a1a1a" } else { "#ffffff" },
        if dark { "#ffffff" } else { "#000000" },
        if dark { "#444444" } else { "#cccccc" },
        if dark { "#2a2a2a" } else { "#f5f5f5" },
        if dark { "#444444" } else { "#cccccc" },
        if dark { "#2a2a2a" } else { "#f5f5f5" },
        if dark { "#444444" } else { "#cccccc" },
        if dark { "#666666" } else { "#999999" },
        png_filename,
        alignments_json.chars().filter(|c| *c == ',').count() + 1,
        paf.targets.len(),
        paf.queries.len(),
        axes.0, axes.1,
        alignments_json,
        targets_json,
        queries_json,
        axes.0,
        axes.1,
        if using_zoom { ranges.0.1 - ranges.0.0 } else { paf.target_length as usize },
        if using_zoom { ranges.1.1 - ranges.1.0 } else { paf.query_length as usize },
        using_zoom,
        format!("[{}, {}]", ranges.0.0, ranges.0.1),
        format!("[{}, {}]", ranges.1.0, ranges.1.1),
        dark,
    );

    let html_filename = png_filename.replace(".png", ".html");
    std::fs::write(&html_filename, html_content)
        .expect("Failed to write HTML file");
    
    println!("Generated interactive HTML viewer: {}", html_filename);
}