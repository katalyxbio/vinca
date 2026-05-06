use anyhow::{bail, Context, Result};
use std::path::PathBuf;
use indicatif::{ProgressBar, ProgressStyle};
use crossbeam_channel::bounded;
use std::collections::BTreeMap;
use std::thread;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::io::Write as _;
use noodles::bam;

use crate::bam_io::{open_bam_reader, open_bam_writer, check_is_name_sorted};
use crate::tags::transfer_tags;

pub struct Processor {
    pub donor_path: PathBuf,
    pub repair_path: PathBuf,
    pub output_path: PathBuf,
    pub threads: usize,
}

struct WorkItem {
    index: usize,
    donor_rec: bam::Record,
    repair_rec: bam::Record,
}

struct ProcessedItem {
    index: usize,
    rec: RecordBuf,
}

impl Processor {
    pub fn new(donor_path: PathBuf, repair_path: PathBuf, output_path: PathBuf, threads: usize) -> Self {
        Self {
            donor_path,
            repair_path,
            output_path,
            threads: threads.max(1),
        }
    }
    
    pub fn run(&self) -> Result<()> {
        println!("Opening BAM files...");
        let (mut donor_reader, donor_header) = open_bam_reader(&self.donor_path)?;
        let (mut repair_reader, repair_header) = open_bam_reader(&self.repair_path)?;
        
        check_is_name_sorted(&donor_header).context("Validating donor BAM")?;
        check_is_name_sorted(&repair_header).context("Validating repair BAM")?;

        println!("Using {} worker threads...", self.threads);
        
        let (tx_in, rx_in) = bounded::<WorkItem>(1000);
        let (tx_out, rx_out) = bounded::<Result<ProcessedItem>>(1000);
        let output_path = self.output_path.clone();
        
        let repair_header_clone = repair_header.clone();
        let writer_handle = thread::spawn(move || -> Result<usize> {
            let mut output_writer = open_bam_writer(&output_path, &repair_header_clone)?;
            let mut count = 0;
            let mut expected_index = 0usize;
            let mut pending = BTreeMap::<usize, RecordBuf>::new();
            
            while let Ok(res) = rx_out.recv() {
                let ProcessedItem { index, rec } = res?;

                if index == expected_index {
                    output_writer.write_alignment_record(&repair_header_clone, &rec)?;
                    count += 1;
                    expected_index += 1;

                    while let Some(next_rec) = pending.remove(&expected_index) {
                        output_writer.write_alignment_record(&repair_header_clone, &next_rec)?;
                        count += 1;
                        expected_index += 1;
                    }
                } else {
                    pending.insert(index, rec);
                }
            }

            if !pending.is_empty() {
                bail!("Writer terminated with {} pending records", pending.len());
            }

            Ok(count)
        });

        let mut worker_handles = Vec::with_capacity(self.threads);
        for _ in 0..self.threads {
            let rx_in = rx_in.clone();
            let tx_out = tx_out.clone();
            let repair_header = repair_header.clone();

            let handle = thread::spawn(move || -> Result<()> {
                while let Ok(item) = rx_in.recv() {
                    match transfer_tags(&item.donor_rec, &repair_header, &item.repair_rec) {
                        Ok(repaired_buf) => {
                            if tx_out
                                .send(Ok(ProcessedItem {
                                    index: item.index,
                                    rec: repaired_buf,
                                }))
                                .is_err()
                            {
                                break;
                            }
                        }
                        Err(e) => {
                            let _ = tx_out.send(Err(e));
                            break;
                        }
                    }
                }

                Ok(())
            });
            worker_handles.push(handle);
        }

        drop(rx_in);
        drop(tx_out);

        let pb = ProgressBar::new_spinner();
        pb.set_style(ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg} {pos} reads processed")?);
        pb.set_message("Processing reads...");

        let mut processed = 0;
        let mut index = 0usize;
        let mut donor_records = donor_reader.records();
        let mut repair_records = repair_reader.records();

        loop {
            let donor_item = donor_records.next();
            let repair_item = repair_records.next();
            
            match (donor_item, repair_item) {
                (Some(Ok(donor_rec)), Some(Ok(repair_rec))) => {
                    if donor_rec.name() != repair_rec.name() {
                        bail!("Mismatched read names: {:?}", donor_rec.name());
                    }

                    tx_in.send(WorkItem {
                        index,
                        donor_rec,
                        repair_rec,
                    })?;
                    index += 1;
                    
                    processed += 1;
                    if processed % 1000 == 0 {
                        pb.set_position(processed);
                    }
                }
                (None, None) => break,
                _ => bail!("Different number of reads in donor and repair BAMs!"),
            }
        }

        pb.set_position(processed);
        drop(tx_in);

        for handle in worker_handles {
            handle.join().expect("Worker thread panicked")?;
        }
        
        let final_count = writer_handle.join().expect("Writer thread panicked")?;
        
        pb.finish_with_message("Processing complete!");
        println!("Successfully processed {} reads.", final_count);
        Ok(())
    }
}
