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
    // each query name in order of first appearance
    queries: Vec<AlignedSeq>,
    // each target name in order of first appearance
    targets: Vec<AlignedSeq>,
    // maps from sequence name to internal id
    query_mphf: Mphf<String>,
    // maps from sequence name to internal id
    target_mphf: Mphf<String>,
    // query axis length
    query_length: f64,
    // target axis length
    target_length: f64,
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

/*
fn paf_target_end(line: &str) -> usize {
    line.split('\t').nth(8).unwrap().parse::<usize>().unwrap()
}
*/

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
        let mut query_idx: usize = 0;
        let mut target_idx: usize = 0;
        let mut query_offset: usize = 0;
        let mut target_offset: usize = 0;
        for_each_line_in_file(filename, |l: &str| {
            let query_name: String = paf_query(l);
            let query_id = query_mphf.hash(&query_name) as usize;
            if !seen_queries[query_id] {
                seen_queries[query_id] = true;
                let query_length = paf_query_length(l);
                let mut query = &mut queries[query_id];
                query.name = query_name;
                query.rank = query_idx;
                query_idx += 1;
                query.length = query_length;
                query.offset = query_offset;
                query_offset += query_length;
            }
            let target_name: String = paf_target(l);
            let target_id = target_mphf.hash(&target_name) as usize;
            if !seen_targets[target_id] {
                seen_targets[target_id] = true;
                let target_length = paf_target_length(l);
                let mut target = &mut targets[target_id];
                target.name = target_name;
                target.rank = target_idx;
                target_idx += 1;
                target.length = target_length;
                target.offset = target_offset;
                target_offset += target_length;
            }
        });
        PafFile {
            filename: filename.to_string(),
            queries,
            targets,
            query_mphf,
            target_mphf,
            query_length: query_offset as f64,
            target_length: target_offset as f64,
        }
    }
    fn query_range(self: &PafFile, name: &str, start: usize, end: usize) -> (usize, usize) {
        println!("query_range {} {}", start, end);
        let query_id = self.query_mphf.hash(&name.into()) as usize;
        let length = self.query_length(query_id);
        let gstart = self.global_query_start(query_id);
        println!("global query start {}", gstart);
        println!("query length {}", length);
        let final_end = gstart + (length - start);
        let final_start = gstart + (length - end);
        println!("hmm {} {}", final_start, final_end);
        (final_start, final_end)
    }
    fn target_range(self: &PafFile, name: &str, start: usize, end: usize) -> (usize, usize) {
        println!("target_range {} {}", start, end);
        let target_id = self.target_mphf.hash(&name.into()) as usize;
        //let length = self.target_length(target_id);
        let gstart = self.global_target_start(target_id);
        println!("global target start {}", gstart);
        let final_start = gstart + start;
        let final_end = gstart + end;
        (final_start, final_end)
    }
    fn global_query_start(self: &PafFile, idx: usize) -> usize {
        self.queries[idx].offset
    }
    fn query_length(self: &PafFile, idx: usize) -> usize {
        self.queries[idx].length
    }
    fn global_target_start(self: &PafFile, idx: usize) -> usize {
        self.targets[idx].offset
    }
    fn target_length(self: &PafFile, idx: usize) -> usize {
        self.targets[idx].length
    }
    fn global_start(self: &PafFile, line: &str, query_rev: bool) -> (usize, usize) {
        let query_id = self.query_mphf.hash(&paf_query(line)) as usize;
        let target_id = self.target_mphf.hash(&paf_target(line)) as usize;
        (
            if query_rev {
                self.global_query_start(query_id)
                    + (self.query_length(query_id) - paf_query_end(line))
            } else {
                self.global_query_start(query_id)
                    + (self.query_length(query_id) - paf_query_begin(line))
            },
            //self.global_query_start(query_id) + paf_query_begin(line),
            self.global_target_start(target_id) + paf_target_begin(line),
        )
    }
    fn for_each_match<F>(self: &PafFile, line: &str, mut func: F)
    where
        F: FnMut(char, usize, bool, usize, usize),
    {
        let query_rev = paf_query_is_rev(line);
        let (x, y) = self.global_start(line, query_rev);
        let mut query_pos = x;
        let mut target_pos = y;
        // find and walk the cigar string
        //println!("{}", line);
        for cigar in line
            .split('\t')
            .skip_while(|s| !s.starts_with("cg:Z:"))
            .map(|s| s.strip_prefix("cg:Z:").unwrap())
            .collect::<Vec<&str>>()
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
                        query_pos += if query_rev { n } else { 0-n };
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
                        query_pos += if query_rev { n } else { 0-n };
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
        if self.query_length > self.target_length {
            let ratio = self.target_length as f64 / self.query_length as f64;
            (major_axis, (major_axis as f64 * ratio) as usize)
        } else {
            let ratio = self.query_length as f64 / self.target_length as f64;
            ((major_axis as f64 * ratio) as usize, major_axis)
        }
    }
    fn get_axes_zoom(self: &PafFile, major_axis: usize, zoom: ((usize, usize), (usize, usize))) -> (usize, usize) {
        let q_length = zoom.0.1 - zoom.0.0;
        let t_length = zoom.1.1 - zoom.1.0;
        if q_length > t_length {
            let ratio = t_length as f64 / q_length as f64;
            (major_axis, (major_axis as f64 * ratio) as usize)
        } else {
            let ratio = q_length as f64 / t_length as f64;
            ((major_axis as f64 * ratio) as usize, major_axis)
        }
    }
    fn project_xy(self: &PafFile, x: usize, y: usize, axes: (usize, usize)) -> (f64, f64) {
        //println!("axes {} {}", axes.0, axes.1);
        (
            (axes.0-1) as f64 * (x as f64 / self.query_length as f64),
            (axes.1-1) as f64 * (y as f64 / self.target_length as f64),
        )
    }
    fn project_xy_zoom(self: &PafFile, x: usize, y: usize, axes: (usize, usize), zoom: ((usize, usize), (usize, usize))) -> (f64, f64) {
        //println!("axes {} {}", axes.0, axes.1);
        (
            (axes.0-1) as f64 * (((x as f64)-(zoom.0.0 as f64)) / ((zoom.0.1 - zoom.0.0) as f64)),
            (axes.1-1) as f64 * (((y as f64)-(zoom.1.0 as f64)) / ((zoom.1.1 - zoom.1.0) as f64)),
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
                .help("Plot the given 2D range in query and target rather than the full matrix seqA:10-200,seqB:300-400"),
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
    let (query_range, target_range): ((usize, usize), (usize, usize)) =
        if using_zoom {
            let splitv = matches.value_of("range").unwrap().split(',').collect::<Vec<&str>>();
            if splitv.len() != 2 {
                panic!("[pafplot::main] invalid range specification {}", matches.value_of("range").unwrap());
            }
            let query_name = splitv[0].split(':').next().unwrap();
            let query_start = splitv[0].split(':').nth(1).unwrap().split('-').next().unwrap().parse::<usize>().unwrap();
            let query_end = splitv[0].split(':').nth(1).unwrap().split('-').nth(1).unwrap().parse::<usize>().unwrap();
            let query_range = paf.query_range(query_name, query_start, query_end);
            let target_name = splitv[1].split(':').next().unwrap();
            let target_start = splitv[1].split(':').nth(1).unwrap().split('-').next().unwrap().parse::<usize>().unwrap();
            let target_end = splitv[1].split(':').nth(1).unwrap().split('-').nth(1).unwrap().parse::<usize>().unwrap();
            let target_range = paf.target_range(target_name, target_start, target_end);

            //self.global_query_start(query_id)
            //    + (self.query_length(query_id) - paf_query_begin(line))

            /*
            let b = paf.project_xy(query_range.0, target_range.0, a);
            let c = paf.project_xy(query_range.1, target_range.1, a);
            query_range = (b.0.round() as usize, c.0.round() as usize);
            target_range = (b.1.round() as usize, c.1.round() as usize);
             */

            (query_range, target_range)
        } else {
            ((0, paf.query_length as usize), (0, paf.target_length as usize))
        };

    // colors we use
    let white = RGB8 {
        r: 255_u8,
        g: 255,
        b: 255,
    };
    let black = white.map(|ch| 255 - ch);

    println!("getting axes query=({} {}) target=({} {})", query_range.0, query_range.1, target_range.0, target_range.1);
    let axes = if using_zoom {
        paf.get_axes_zoom(major_axis, (query_range, target_range))
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

    let get_coords =
        |x: usize, y: usize| if using_zoom {
            paf.project_xy_zoom(x, y, axes, (query_range, target_range))
        } else {
            paf.project_xy(x, y, axes)
        };

    // for each match, we draw a line on our raster using Xiaolin Wu's antialiased line algorithm
    let draw_match = |_c, x: usize, rev: bool, y: usize, len: usize| {
        println!("draw_match {} {} {} {} {}", _c, x, rev, y, len);
        let start = get_coords(x, y);
        let end = get_coords(x + if rev { len } else { 0-len }, y + len);
        println!("start and end ({} {}) ({} {})", start.0, start.1, end.0, end.1);
        for ((x, y), val) in XiaolinWu::<f64, i64>::new(start, end) {
            if x >= 0 && x < (axes.0 as i64) && y >= 0 && y < (axes.1 as i64) {
                let i: usize = (x as usize) + (y as usize * axes.0);
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
}
