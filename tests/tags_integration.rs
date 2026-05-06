use anyhow::Result;
use noodles::bam;

#[test]
fn test_length_change_detects() -> Result<()> {
    // Scaffold test for MM/ML detection loop bounds
    let donor_rec = bam::Record::default();
    let repair_rec = bam::Record::default();
    
    assert_eq!(donor_rec.sequence().len(), repair_rec.sequence().len());
    Ok(())
}
