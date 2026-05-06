use clap::Parser;
use std::num::NonZeroUsize;
use std::path::PathBuf;
use vinca::processor::Processor;
use anyhow::Result;

#[derive(Parser, Debug)]
#[command(
    name = "vinca",
    author,
    version,
    about = "Transfer and recalculate MM/ML methylation tags between BAM files",
    long_about = None
)]
struct Args {
    #[arg(short, long, value_name = "DONOR_BAM", help = "Donor BAM containing source MM/ML tags")]
    donor: PathBuf,
    
    #[arg(short, long, value_name = "REPAIR_BAM", help = "Repair BAM to receive recalculated MM/ML tags")]
    repair: PathBuf,
    
    #[arg(short, long, value_name = "OUTPUT_BAM", help = "Output BAM path")]
    output: PathBuf,

    #[arg(
        short = 't',
        long,
        value_name = "N",
        help = "Number of worker threads for tag transfer",
        default_value_t = default_threads()
    )]
    threads: usize,
}

fn default_threads() -> usize {
    std::thread::available_parallelism()
        .map(NonZeroUsize::get)
        .unwrap_or(1)
}

fn main() -> Result<()> {
    env_logger::init();
    let args = Args::parse();
    
    let processor = Processor::new(
        args.donor,
        args.repair,
        args.output,
        args.threads,
    );
    processor.run()?;
    
    println!("Done!");
    Ok(())
}
