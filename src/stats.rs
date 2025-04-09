use anyhow::{Context, Result};
use bio::io::fasta::Sequence;
use bio::utils;
use bstr::ByteSlice;
use flate2;
use itertools::{self, Itertools};
use log::{info, warn};
use noodles::fastq;
use noodles::sam::alignment::io::Write;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::Cigar;
use noodles::sam::header::record::value::map::header;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use serde::{de, Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::io::{BufRead};
use std::path::{Path, PathBuf};


#[derive(Debug, Serialize, Clone, Deserialize)]
struct LigationStats{
    n_cis: u64,
    n_trans: u64,
    n_total: u64,
}



pub fn get_ligation_stats(
    bam: &str,
    output: &str,
) -> Result<()>{


    let mut bam = noodles::bam::io::reader::Builder::default().build_from_path(bam)?;

    for record in bam.records() {
        let record = record?;

        // Extract the oligo coordinate tag from the read






        
    }













    Ok(())

}


