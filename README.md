# Vinca

`Vinca` is a Rust command-line tool for transferring and recalculating methylation tags (`MM`/`ML`) from a donor BAM to a repair BAM.

It is implemented using pure Rust with the `noodles` BAM/SAM/BGZF stack and processes alignments in a streaming multithreaded pipeline.

## Features

- Transfers methylation metadata from a donor BAM to a repair BAM
- Recalculates `MM`/`ML` tag content when clipping or sequence offsets differ
- Uses a pure Rust implementation with `noodles`
- Multithreaded tag transfer workers with ordered BAM output
- CLI with donor, repair, output, and worker thread controls

## Requirements

- Rust toolchain (recommended: latest stable)
- `cargo`

## Installation

From the repository root:

```bash
cargo install --path .
```

This installs `vinca` to `$HOME/.cargo/bin` (or your configured Cargo bin directory).

For local development without installing:

```bash
cargo build --release
```

The compiled binary will be at `target/release/vinca`.

## Usage

```bash
vinca --donor <donor.bam> --repair <repair.bam> --output <output.bam>
```

Show help:

```bash
vinca --help
```

Tested help output:

```text
Transfer and recalculate MM/ML methylation tags between BAM files

Usage: vinca [OPTIONS] --donor <DONOR_BAM> --repair <REPAIR_BAM> --output <OUTPUT_BAM>
```

### CLI arguments

- `-d`, `--donor`   : Path to the donor BAM file containing the original methylation tags.
- `-r`, `--repair`  : Path to the repair BAM file that should receive the recalculated tags.
- `-o`, `--output`  : Path to the output BAM file to write.
- `-t`, `--threads` : Number of worker threads used for tag transfer. Defaults to available CPU parallelism.

## Example

```bash
vinca \
  --donor donor.bam \
  --repair repair.bam \
  --output repaired.bam
```

Run with an explicit worker count:

```bash
vinca \
  --donor donor.bam \
  --repair repair.bam \
  --output repaired.bam \
  --threads 8
```

If `vinca` is not found in your shell, add Cargo's bin directory to `PATH`:

```bash
export PATH="$HOME/.cargo/bin:$PATH"
```

## Notes

- The donor and repair BAM files should contain matching read names so the tool can transfer tags correctly.
- The output is written as a BAM file.
- Records are processed in parallel and written in the original input order.

## Project structure

- `src/bin/vinca.rs`  - CLI entrypoint
- `src/processor.rs`  - Streaming read/repair processing pipeline
- `src/tags.rs`       - Tag transfer and recalculation logic
- `src/bam_io.rs`     - BAM header and I/O helpers

## License

This repository does not include a license file. Add one if you plan to distribute or publish the tool.
