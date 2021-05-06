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

struct PafFile {
    // our input file
    paf_filename: String,
    // each query name in order of first appearance
    queries: Vec<AlignedSeq>,
    // each target name in order of first appearance
    targets: Vec<AlignedSeq>,
    // maps from sequence name to internal id
    query_mphf: Mphf<String>,
    // maps from sequence name to internal id
    target_mphf: Mphf<String>,
}

impl PafFile {
    fn new(paf_filename: &str) -> Self {
        let mut query_names: Vec<String> = Vec::new();
        let mut target_names: Vec<String> = Vec::new();
        for_each_line_in_file(paf_filename, |l: &str| {
            query_names.push(l.split('\t').next().unwrap().into());
            target_names.push(l.split('\t').nth(5).unwrap().into());
            
        });
        query_names = query_names.into_iter().sorted().dedup().collect();
        target_names = target_names.into_iter().sorted().dedup().collect();
        let query_mphf = Mphf::new(1.7, &query_names);
        let target_mphf = Mphf::new(1.7, &target_names);
        let mut seen_queries = vec![false; query_names.len()];
        let mut seen_targets = vec![false; target_names.len()];
        let mut queries: Vec<AlignedSeq> = Vec::<AlignedSeq>::new();
        let mut targets: Vec<AlignedSeq> = Vec::<AlignedSeq>::new();
        let mut query_idx: usize = 0;
        let mut target_idx: usize = 0;
        let mut query_offset: usize = 0;
        let mut target_offset: usize = 0;
        for_each_line_in_file(paf_filename, |l: &str| {
            let query_name: String = l.split('\t').next().unwrap().into();
            let query_id = query_mphf.hash(&query_name) as usize;
            if !seen_queries[query_id] {
                seen_queries[query_id] = true;
                let query_length = l.split('\t').nth(1).unwrap().parse::<usize>().unwrap();
                let mut query = &mut queries[query_id];
                query.name = query_name;
                query.rank = query_idx;
                query_idx += 1;
                query.length = query_length;
                query.offset = query_offset;
                query_offset += query_length;
            }
            let target_name: String = l.split('\t').nth(5).unwrap().into();
            let target_id = target_mphf.hash(&target_name) as usize;
            if !seen_targets[target_id] {
                seen_targets[target_id] = true;
                let target_length = l.split('\t').nth(1).unwrap().parse::<usize>().unwrap();
                let mut target = &mut targets[target_id];
                target.name = target_name;
                target.rank = target_idx;
                target_idx += 1;
                target.length = target_length;
                target.offset = target_offset;
                target_offset += target_length;
            }
        });
        //let mut paf = PafFile::new();
        PafFile {
            paf_filename: paf_filename.to_string(),
            queries,
            targets,
            query_mphf,
            target_mphf
        }
    }
    /*
    fn get_id(self: &PafFile, name: &str) -> u64 {
        self.seq_name_mphf.hash(&name.to_string())
    }
*/
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
}
