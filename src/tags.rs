use anyhow::{bail, Context, Result};
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::data::field::Value as RecordValue;
use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
use noodles::sam::alignment::record_buf::data::field::value::Array;
use noodles::sam::alignment::RecordBuf;
use noodles::sam;
use noodles::bam;

pub const MM_TAG: Tag = Tag::new(b'M', b'M');
pub const ML_TAG: Tag = Tag::new(b'M', b'L');

#[derive(Debug, Clone, PartialEq, Eq)]
struct RepairedMmGroup {
    header: String,
    skip_counts: Vec<usize>,
    ml_values: Vec<u8>,
}

fn find_repair_start(donor_bases: &[u8], repair_bases: &[u8]) -> Result<usize> {
    if repair_bases.is_empty() {
        return Ok(0);
    }

    if donor_bases.len() < repair_bases.len() {
        bail!("Donor sequence is shorter than repair sequence");
    }

    if donor_bases.len() == repair_bases.len() {
        if donor_bases == repair_bases {
            return Ok(0);
        }

        bail!("Donor and repair sequences differ but have the same length");
    }

    donor_bases
        .windows(repair_bases.len())
        .position(|window| window == repair_bases)
        .ok_or_else(|| anyhow::anyhow!("Repair sequence was not found as a contiguous substring of the donor sequence"))
}

fn canonical_positions(sequence: &[u8], base: u8) -> Vec<usize> {
    let base = base.to_ascii_uppercase();

    sequence
        .iter()
        .enumerate()
        .filter_map(|(index, observed)| {
            let observed = observed.to_ascii_uppercase();

            if base == b'N' || observed == base {
                Some(index)
            } else {
                None
            }
        })
        .collect()
}

fn parse_mm_group(group: &str) -> Result<(&str, Vec<usize>)> {
    let (header, counts) = group
        .split_once(',')
        .with_context(|| format!("MM group is missing skip counts: {group}"))?;

    let skip_counts = if counts.is_empty() {
        Vec::new()
    } else {
        counts
            .split(',')
            .map(|value| {
                value
                    .parse::<usize>()
                    .with_context(|| format!("Invalid MM skip count in group {group}: {value}"))
            })
            .collect::<Result<Vec<_>>>()?
    };

    Ok((header, skip_counts))
}

fn group_positions_from_skip_counts(
    canonical_positions: &[usize],
    skip_counts: &[usize],
    group: &str,
) -> Result<Vec<usize>> {
    let mut positions = Vec::with_capacity(skip_counts.len());
    let mut canonical_index = 0usize;

    for skip in skip_counts {
        canonical_index = canonical_index
            .checked_add(*skip)
            .ok_or_else(|| anyhow::anyhow!("MM skip count overflow in group {group}"))?;

        // Some real-world files can retain MM offsets that extend beyond the current
        // visible canonical-base list (e.g., after clipping/hard clipping). We keep
        // the valid prefix and drop the trailing calls we cannot represent.
        let Some(position) = canonical_positions.get(canonical_index).copied() else {
            break;
        };

        positions.push(position);
        canonical_index += 1;
    }

    Ok(positions)
}

fn repair_group(
    header: &str,
    donor_positions: &[usize],
    donor_ml_values: &[u8],
    repair_canonical_positions: &[Option<usize>],
    repair_start: usize,
    repair_len: usize,
) -> Result<Option<RepairedMmGroup>> {
    let mut repaired_positions = Vec::new();
    let mut repaired_ml_values = Vec::new();

    for (local_index, donor_position) in donor_positions.iter().copied().enumerate() {
        if donor_position < repair_start || donor_position >= repair_start + repair_len {
            continue;
        }

        let repair_position = donor_position - repair_start;
        // If the repaired read no longer has the expected canonical base at this
        // query position, drop the call for this site.
        let Some(canonical_index) = repair_canonical_positions
            .get(repair_position)
            .and_then(|index| *index)
        else {
            continue;
        };

        repaired_positions.push(canonical_index);

        if let Some(ml_value) = donor_ml_values.get(local_index) {
            repaired_ml_values.push(*ml_value);
        }
    }

    if repaired_positions.is_empty() {
        return Ok(None);
    }

    let mut skip_counts = Vec::with_capacity(repaired_positions.len());
    let mut previous = 0usize;

    for (index, canonical_index) in repaired_positions.iter().copied().enumerate() {
        if index == 0 {
            skip_counts.push(canonical_index);
        } else {
            skip_counts.push(canonical_index - previous - 1);
        }

        previous = canonical_index;
    }

    Ok(Some(RepairedMmGroup {
        header: header.to_string(),
        skip_counts,
        ml_values: repaired_ml_values,
    }))
}

fn repair_mm_ml_tags(
    donor_mm: &str,
    donor_ml_values: Option<&[u8]>,
    donor_bases: &[u8],
    repair_bases: &[u8],
    repair_start: usize,
) -> Result<Option<(String, Option<Vec<u8>>)>> {
    let donor_ml_values = donor_ml_values.unwrap_or(&[]);
    let mut ml_cursor = 0usize;
    let mut repaired_groups = Vec::new();
    let mut repaired_ml = Vec::new();

    for raw_group in donor_mm.split(';').filter(|group| !group.is_empty()) {
        let (header, skip_counts) = parse_mm_group(raw_group)?;
        let header_base = header
            .as_bytes()
            .first()
            .copied()
            .with_context(|| format!("MM group is missing a canonical base: {raw_group}"))?;

        let donor_positions = group_positions_from_skip_counts(
            &canonical_positions(donor_bases, header_base),
            &skip_counts,
            raw_group,
        )?;

        let mut repair_canonical_map = vec![None; repair_bases.len()];
        for (canonical_index, repair_position) in canonical_positions(repair_bases, header_base)
            .into_iter()
            .enumerate()
        {
            repair_canonical_map[repair_position] = Some(canonical_index);
        }

        let ml_end = ml_cursor.saturating_add(skip_counts.len());
        let group_ml_values = donor_ml_values
            .get(ml_cursor..ml_end)
            .unwrap_or(&[]);
        ml_cursor = ml_end;

        if let Some(repaired_group) = repair_group(
            header,
            &donor_positions,
            group_ml_values,
            &repair_canonical_map,
            repair_start,
            repair_bases.len(),
        )? {
            repaired_ml.extend(repaired_group.ml_values.iter().copied());
            repaired_groups.push(repaired_group);
        }
    }

    if repaired_groups.is_empty() {
        return Ok(None);
    }

    let repaired_mm = repaired_groups
        .into_iter()
        .map(|group| {
            let counts = group
                .skip_counts
                .into_iter()
                .map(|count| count.to_string())
                .collect::<Vec<_>>()
                .join(",");

            format!("{},{}", group.header, counts)
        })
        .collect::<Vec<_>>()
        .join(";")
        + ";";

    let repaired_ml = if donor_ml_values.is_empty() {
        None
    } else {
        Some(repaired_ml)
    };

    Ok(Some((repaired_mm, repaired_ml)))
}

pub fn transfer_tags(
    donor_record: &bam::Record,
    repair_header: &sam::Header,
    repair_record: &bam::Record,
) -> Result<RecordBuf> {
    let donor_bases: Vec<u8> = donor_record.sequence().iter().collect();
    let repair_bases: Vec<u8> = repair_record.sequence().iter().collect();

    let start_idx = find_repair_start(&donor_bases, &repair_bases)?;

    // Convert repair record to editable buffer
    let mut repaired_buf = RecordBuf::try_from_alignment_record(repair_header, repair_record)?;
    
    let donor_mm = match donor_record.data().get(&MM_TAG) {
        Some(Ok(RecordValue::String(s))) => std::str::from_utf8(s.as_ref())
            .context("Donor MM tag contains invalid UTF-8")?,
        Some(Ok(_)) => bail!("Donor MM tag is not a string"),
        Some(Err(e)) => return Err(e.into()),
        None => return Ok(repaired_buf),
    };

    let donor_ml_values = match donor_record.data().get(&ML_TAG) {
        Some(Ok(RecordValue::Array(array))) => match array {
            noodles::sam::alignment::record::data::field::value::Array::UInt8(values) => {
                let mut out = Vec::with_capacity(values.len());
                for value in values.iter() {
                    out.push(value?);
                }
                Some(out)
            }
            _ => bail!("Donor ML tag must be a uint8 array"),
        },
        Some(Ok(_)) => bail!("Donor ML tag is not an array"),
        Some(Err(e)) => return Err(e.into()),
        None => None,
    };

    if let Some((mm, ml)) = repair_mm_ml_tags(
        donor_mm,
        donor_ml_values.as_deref(),
        &donor_bases,
        &repair_bases,
        start_idx,
    )? {
        repaired_buf
            .data_mut()
            .insert(MM_TAG, BufValue::String(mm.into()));

        if let Some(ml) = ml {
            repaired_buf
                .data_mut()
                .insert(ML_TAG, BufValue::Array(Array::UInt8(ml)));
        }
    }
    
    Ok(repaired_buf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn repairs_mm_and_ml_after_prefix_trim() -> Result<()> {
        let donor_bases = b"AACCGGTT";
        let repair_bases = b"CCGG";
        let start_idx = find_repair_start(donor_bases, repair_bases)?;

        let repaired = repair_mm_ml_tags(
            "A+m,0;C+m,0",
            Some(&[11, 22]),
            donor_bases,
            repair_bases,
            start_idx,
        )?
        .expect("expected repaired tags");

        assert_eq!(repaired.0, "C+m,0;");
        assert_eq!(repaired.1, Some(vec![22]));
        Ok(())
    }

    #[test]
    fn drops_fully_clipped_groups() -> Result<()> {
        let donor_bases = b"AACCGGTT";
        let repair_bases = b"CCGG";
        let start_idx = find_repair_start(donor_bases, repair_bases)?;

        let repaired = repair_mm_ml_tags(
            "A+m,0;C+m,0",
            Some(&[1, 2]),
            donor_bases,
            repair_bases,
            start_idx,
        )?
        .expect("expected repaired tags");

        assert_eq!(repaired.0, "C+m,0;");
        assert_eq!(repaired.1, Some(vec![2]));
        Ok(())
    }

    #[test]
    fn returns_none_when_all_groups_are_clipped() -> Result<()> {
        let donor_bases = b"AACCGGTT";
        let repair_bases = b"GGTT";
        let start_idx = find_repair_start(donor_bases, repair_bases)?;

        let repaired = repair_mm_ml_tags(
            "A+m,0;C+m,0",
            Some(&[1, 2]),
            donor_bases,
            repair_bases,
            start_idx,
        )?;

        assert!(repaired.is_none());
        Ok(())
    }

    #[test]
    fn tolerates_out_of_range_skip_counts() -> Result<()> {
        let donor_bases = b"CCCCC";
        let repair_bases = b"CCCCC";
        let start_idx = find_repair_start(donor_bases, repair_bases)?;

        let repaired = repair_mm_ml_tags(
            "C+h?,0,1,10",
            Some(&[10, 20, 30]),
            donor_bases,
            repair_bases,
            start_idx,
        )?
        .expect("expected repaired tags");

        // The trailing skip count points past canonical C positions and is dropped.
        assert_eq!(repaired.0, "C+h?,0,1;");
        assert_eq!(repaired.1, Some(vec![10, 20]));
        Ok(())
    }
}
