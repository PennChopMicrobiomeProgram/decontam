import pysam

def get_mapped_reads(fp, min_pct_id=0.5, min_len=100):
    sam = pysam.AlignmentFile(fp)
    for read in sam:
        if read.is_unmapped:
            yield (read.query_name, read.is_read1, None)
        elif read.query_alignment_length < min_len:
            yield (read.query_name, read.is_read1, None)
        else:
            pct_id = _get_pct_identity(read)
            if pct_id < min_pct_id:
                yield (read.query_name, read.is_read1, None)
            else:
                ref_id = sam.getrname(read.reference_id)
                yield (read.query_name, read.is_read1, ref_id)


def _get_pct_identity(read):
    if read.has_tag("NM"):
        edit_dist = read.get_tag("NM")
    else:
        edit_dist = 0
    pct_mm = float(edit_dist) / read.alen
    return 1 - pct_mm 
