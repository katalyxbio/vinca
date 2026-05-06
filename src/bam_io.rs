use anyhow::{bail, Context, Result};
use noodles::bam;
use noodles::bgzf;
use noodles::sam;
use std::path::Path;

pub fn check_is_name_sorted(header: &sam::Header) -> Result<()> {
    // A more robust way to check if there's any map with SO:queryname.
    // But since API changes so much, we can just print the header using Sam writer into a string buffer
    let mut buf = Vec::new();
    let mut writer = sam::io::Writer::new(&mut buf);
    writer.write_header(header).context("Failed to write header to buffer")?;
    
    let raw_header = String::from_utf8_lossy(&buf);
    if raw_header.contains("SO:queryname") {
        return Ok(());
    }

    bail!("BAM file is not name sorted (`SO:queryname`). Please sort it by name (e.g., using `samtools sort -n`) before running Vinca.")
}

pub fn open_bam_reader<P: AsRef<Path>>(
    path: P,
) -> Result<(bam::io::Reader<bgzf::io::Reader<std::fs::File>>, sam::Header)> {
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(path)
        .context("Failed to open BAM file for reading")?;
    
    let header = reader.read_header().context("Failed to read BAM header")?;
    Ok((reader, header))
}

pub fn open_bam_writer<P: AsRef<Path>>(
    path: P,
    header: &sam::Header,
) -> Result<bam::io::Writer<bgzf::io::Writer<std::fs::File>>> {
    let mut writer = bam::io::writer::Builder::default()
        .build_from_path(path)
        .context("Failed to open BAM file for writing")?;
        
    writer
        .write_header(header)
        .context("Failed to write BAM header")?;
        
    Ok(writer)
}
