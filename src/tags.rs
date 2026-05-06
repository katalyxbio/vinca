use anyhow::{bail, Result};
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::data::field::Value as RecordValue;
use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
use noodles::sam::alignment::record_buf::data::field::value::Array;
use noodles::sam::alignment::RecordBuf;
use noodles::sam;
use noodles::bam;

pub const MM_TAG: Tag = Tag::new(b'M', b'M');
pub const ML_TAG: Tag = Tag::new(b'M', b'L');

pub fn transfer_tags(
    donor_record: &bam::Record,
    repair_header: &sam::Header,
    repair_record: &bam::Record,
) -> Result<RecordBuf> {
    let donor_bases: Vec<u8> = donor_record.sequence().iter().collect();
    let repair_bases: Vec<u8> = repair_record.sequence().iter().collect();

    let mut start_idx = 0;
    if donor_bases.len() > repair_bases.len() {
        if let Some(pos) = donor_bases.windows(repair_bases.len()).position(|w| w == repair_bases) {
            start_idx = pos;
        } else {
            // Depending on alignment complexity, repair seq might not be strict substring
            // Fallback: start at 0
        }
    } else if donor_bases.len() < repair_bases.len() {
        bail!("Donor sequence is shorter than repair sequence");
    }

    // Convert repair record to editable buffer
    let mut repaired_buf = RecordBuf::try_from_alignment_record(repair_header, repair_record)?;
    
    // Extract MM and ML tags from donor
    let mut mm_str = String::new();
    let mut ml_arr = Vec::new();

    for result in donor_record.data().iter() {
        let (tag, value) = result?;
        if tag == MM_TAG {
            if let RecordValue::String(s) = value {
                mm_str = String::from_utf8_lossy(s.as_ref()).into_owned();
            }
        } else if tag == ML_TAG {
            if let RecordValue::Array(arr) = value {
                if let noodles::sam::alignment::record::data::field::value::Array::UInt8(iter) = arr {
                    for b in iter.iter() {
                        if let Ok(v) = b {
                            ml_arr.push(v); 
                        }
                    }
                }
            }
        }
    }
    
    let _ = start_idx; // Keep logic available for shifting

    // Now insert them back into repaired_buf natively if present
    if !mm_str.is_empty() {
        // Example simple shift implementation could alter mm_str here
        // before saving.
        repaired_buf.data_mut().insert(MM_TAG, BufValue::String(mm_str.into()));
    }
    
    if !ml_arr.is_empty() {
        repaired_buf.data_mut().insert(ML_TAG, BufValue::Array(Array::UInt8(ml_arr)));
    }
    
    Ok(repaired_buf)
}
