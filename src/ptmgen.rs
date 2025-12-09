use std::collections::HashSet;
use std::error::Error;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};

/*
Gaurav Sablok
codeprog@icloud.com
*/

pub fn ptmgenerate(pathfile: &str) -> Result<String, Box<dyn Error>> {
    let pathfileopen = File::open(pathfile).expect("file not present");
    let pathfileread = BufReader::new(pathfileopen);
    let mut head: Vec<String> = Vec::new();
    let mut seq: Vec<String> = Vec::new();
    for i in pathfileread.lines() {
        let line = i.expect("line not present");
        if line.starts_with("#") {
            head.push(line.clone())
        }
        if !line.starts_with("#") {
            seq.push(line)
        }
    }
    let mut sequencevec: HashSet<String> = HashSet::new();
    for i in seq.iter() {
        let kmervec = i
            .as_bytes()
            .windows(3)
            .map(|x| std::str::from_utf8(x).unwrap())
            .collect::<Vec<_>>();
        for kmer in kmervec.iter() {
            sequencevec.insert(kmer.to_string());
        }
    }
    let mut filewrite = File::create("kmerptms.txt").expect("file not present");
    for i in sequencevec.iter() {
        writeln!(filewrite, "{}\n", i).expect("file not present");
    }
    let returnstrin = String::from("file has been written");
    Ok(returnstrin)
}
