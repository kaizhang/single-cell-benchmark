use bed_utils::{bed, bed::{BED, GenomicRange, BEDLike}};
use flate2::{
    Compression,
    read::MultiGzDecoder,
    write::GzEncoder,
};
use std::fs::File;
use std::io::{BufRead, BufWriter, BufReader, Write};
use itertools::Itertools;
use std::env;

fn read_bed<'a, W: Write>(
    fl: &'a str,
    n: usize,
    writer: &'a mut bed::io::Writer<W>,
) {
    bed::io::Reader::new(
        MultiGzDecoder::new(File::open(fl).unwrap()),
        Some("#".to_string()),
    )
        .into_records()
        .map(Result::unwrap)
        .group_by(|x: &BED<5>| x.name().unwrap().to_string())
        .into_iter()
        .take(n)
        .flat_map(|x| x.1)
        .for_each(|r| writer.write_record(&r).unwrap());
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let input_file = args[1].as_str();
    let output_file = args[2].as_str();
    let n = args[3].parse().unwrap();

    let output = File::create(output_file).unwrap();
    let mut writer = bed::io::Writer::new(
        GzEncoder::new(BufWriter::new(output), Compression::default())
	);
    read_bed(input_file, n, &mut writer);
}
