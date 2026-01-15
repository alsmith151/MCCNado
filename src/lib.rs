use log::info;
use pyo3::prelude::*;

mod bam_deduplicate;
mod fastq_deduplicate;
mod ligation_stats;
mod mcc_data_handler;
mod utils;
mod viewpoint_read_splitter;

/// Deduplicates a FASTQ or FASTQ pair.
#[pyfunction]
#[pyo3(signature = (fastq1, output1, fastq2=None, output2=None))]
#[pyo3(text_signature = "(fastq1, output1, fastq2=None, output2=None)")]
fn deduplicate_fastq(
    fastq1: &str,
    output1: &str,
    fastq2: Option<&str>,
    output2: Option<&str>,
) -> PyResult<fastq_deduplicate::FastqDeduplicationStats> {
    let deduplicator = fastq_deduplicate::DuplicateRemover::from_fastq_paths(fastq1, fastq2);
    let mut deduplicator = match deduplicator {
        Err(e) => {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                e.to_string(),
            ))
        }
        Ok(_) => deduplicator.unwrap(),
    };

    let res = deduplicator.deduplicate(output1, output2);

    match res {
        Err(e) => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            e.to_string(),
        )),
        Ok(_) => Ok(res.unwrap()),
    }
}

/// Deduplicates a BAM file based on segment coordinates.
#[pyfunction]
#[pyo3(signature = (bam, output))]
fn deduplicate_bam(bam: &str, output: &str) -> PyResult<bam_deduplicate::BamDeduplicationStats> {
    let res = bam_deduplicate::deduplicate_bam(bam, output);

    match res {
        Err(e) => {
            log::error!("{}", e);
            Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                e.to_string(),
            ))
        }
        Ok(stats) => Ok(stats),
    }
}

#[pyfunction]
#[pyo3(signature = (bam, output))]
fn split_viewpoint_reads(bam: &str, output: &str) -> PyResult<()> {
    let splitter_options = viewpoint_read_splitter::ReadSplitterOptions::default();
    let splitter = viewpoint_read_splitter::ReadSplitter::new(bam, splitter_options);
    let res = splitter.split_reads(output);

    match res {
        Err(e) => {
            log::error!("{}", e);
            Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                e.to_string(),
            ))
        }
        Ok(_) => Ok(()),
    }
}

#[pyfunction]
#[pyo3(signature = (bam, output_directory))]
fn identify_ligation_junctions(bam: &str, output_directory: &str) -> PyResult<()> {
    let res = mcc_data_handler::identify_ligation_junctions(bam, output_directory);

    match res {
        Err(e) => {
            log::error!("{}", e);
            Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                e.to_string(),
            ))
        }
        Ok(_) => Ok(()),
    }
}

#[pyfunction]
#[pyo3(signature = (bam, output))]
fn annotate_bam(bam: &str, output: &str) -> PyResult<()> {
    let res = mcc_data_handler::annotate_bam(bam, output);

    match res {
        Err(e) => {
            log::error!("{}", e);
            Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                e.to_string(),
            ))
        }
        Ok(_) => Ok(()),
    }
}

#[pyfunction]
#[pyo3(signature = (bam, stats))]
fn extract_ligation_stats(bam: &str, stats: &str) -> PyResult<()> {
    let res = ligation_stats::get_ligation_stats(bam, stats);

    match res {
        Err(e) => {
            log::error!("{}", e);
            Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                e.to_string(),
            ))
        }
        Ok(_) => Ok(()),
    }
}

/// Rust implementation of MCC code.
#[pymodule]
fn mccnado(m: &Bound<'_, PyModule>) -> PyResult<()> {
    pyo3_log::init();
    let ctrlc_handler = ctrlc::try_set_handler(move || {
        info!("Received SIGINT, exiting...");
        std::process::exit(0);
    });

    match ctrlc_handler {
        Ok(_) => (),
        Err(e) => {
            // if it errors then just ignore it
            info!("Failed to set up SIGINT handler: {}", e);
        }
    }

    m.add_function(wrap_pyfunction!(deduplicate_fastq, m)?)?;
    m.add_function(wrap_pyfunction!(deduplicate_bam, m)?)?;
    m.add_function(wrap_pyfunction!(split_viewpoint_reads, m)?)?;
    m.add_function(wrap_pyfunction!(annotate_bam, m)?)?;
    m.add_function(wrap_pyfunction!(identify_ligation_junctions, m)?)?;
    m.add_function(wrap_pyfunction!(extract_ligation_stats, m)?)?;
    m.add_class::<bam_deduplicate::BamDeduplicationStats>()?;
    Ok(())
}
