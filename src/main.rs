use std::fs::File;
//use std::collections::HashMap;
use std::io::{self, prelude::*, BufReader};

use boomphf::*;

use itertools::Itertools;

extern crate clap;
use clap::{Arg, App}; //, SubCommand};

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
            name: String::new(), length: 0, rank: 0, offset: 0
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
        let mut queries: Vec<AlignedSeq> = vec![AlignedSeq::new() ; query_names.len()];
        let mut targets: Vec<AlignedSeq> = vec![AlignedSeq::new() ; target_names.len()];
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
            target_mphf
        }
    }
    fn global_query_start(self: &PafFile, idx: usize) -> usize {
        self.queries[idx].offset
    }
    fn global_target_start(self: &PafFile, idx: usize) -> usize {
        self.targets[idx].offset
    }
    fn global_start(self: &PafFile, line: &str) -> (usize, usize) {
        let query_id = self.query_mphf.hash(&paf_query(line)) as usize;
        let target_id = self.target_mphf.hash(&paf_target(line)) as usize;
        (self.global_query_start(query_id), self.global_target_start(target_id))
    }
    fn process(self: &PafFile) {
        for_each_line_in_file(&self.filename, |line: &str| {
            let (x, y) = self.global_start(line);
            println!("{} {} {} {} {} {}",
                     paf_query(line),
                     paf_query_begin(line),
                     paf_target(line),
                     paf_target_begin(line),
                     x, y);
        })
    }
}

fn main() {
    let matches = App::new("pafplot")
        .version("0.1.0")
        .author("Erik Garrison <erik.garrison@gmail.com>")
        .about("Generate a dotplot from pairwise DNA alignments in PAF format")
        .arg(Arg::with_name("INPUT")
             .required(true)
             .takes_value(true)
             .index(1)
             .help("input PAF file"))
        .arg(Arg::with_name("png")
             .short("p")
             .long("png")
             .help("Save the dotplot to this file."))
        .get_matches();
    let filename = matches.value_of("INPUT").unwrap();
    let paf = PafFile::new(filename);
    paf.process()
}
