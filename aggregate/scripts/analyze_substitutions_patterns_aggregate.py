#!/usr/bin/env python3
"""
Analyze substitutions in BAM/SAM files.

This script:
1. Checks that there is only one alignment per read
2. Applies strand correction to convert reverse reads to original orientation
3. Calculates observed substitutions by comparing reads to reference
4. Reports base counts/frequencies for observed and reference bases
5. Displays substitution tables ordered by frequency
6. Splits analysis by tax ID (from read.reference_name / RNAME field)
"""

import argparse
import sys
import os
from collections import defaultdict, Counter
import pysam


BASES = ['A', 'C', 'G', 'T']


def complement(base):
    """Return complement of a base."""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return comp.get(base, 'N')


def _make_taxid_stats():
    """Return a fresh per-taxid stats dict."""
    return {
        'num_reads': 0,
        'num_reverse': 0,
        'substitutions': defaultdict(int),   # (ref_base, obs_base) -> count (includes matches)
        'observed_bases': Counter(),
        'reference_bases': Counter(),
    }


def _csv_header_parts(include_taxid=False):
    """Return ordered list of CSV column names for the summary CSV."""
    parts = ['filename']
    if include_taxid:
        parts.append('tax_id')
    parts += ['total_bases', 'A_obs', 'C_obs', 'G_obs', 'T_obs', 'total_mismatches']
    parts += ['A_ref', 'C_ref', 'G_ref', 'T_ref']
    for ref in BASES:
        for obs in BASES:
            if ref != obs:
                parts.append(f"{ref}>{obs}")
    return parts


def _reads_csv_header_parts():
    """Return ordered list of CSV column names for the per-read CSV."""
    parts = ['read_name', 'tax_id', 'ref_A', 'ref_C', 'ref_G', 'ref_T']
    for ref in BASES:
        for obs in BASES:
            if ref != obs:
                parts.append(f"{ref}>{obs}")
    return parts


def analyze_bam(bam_file, error_log_file=None, strand_correction=True, reads_out=None):
    """
    Analyze BAM/SAM file for substitutions and base frequencies, split by tax ID.

    Args:
        bam_file: Path to BAM or SAM file
        error_log_file: Optional path to write failed reads
        strand_correction: If True, convert reverse-strand reads to original orientation
        reads_out: Optional file handle for streaming per-read CSV output

    Returns:
        (global_stats, stats_by_taxid)
        global_stats keys: filename, num_alignments, reads_seen, duplicate_reads,
                           failed_reads, strand_correction
        stats_by_taxid: dict keyed by tax_id string -> per-taxid stats dict
    """
    global_stats = {
        'filename': bam_file,
        'num_alignments': 0,
        'reads_seen': set(),
        'duplicate_reads': [],
        'failed_reads': 0,
        'strand_correction': strand_correction,
    }
    stats_by_taxid = {}

    # Determine file type and open mode
    file_ext = os.path.splitext(bam_file)[1].lower()
    if file_ext == '.sam':
        mode = 'r'
    elif file_ext == '.bam':
        mode = 'rb'
    else:
        print(f"Warning: Unrecognized file extension '{file_ext}'. Attempting to auto-detect format.", file=sys.stderr)
        mode = 'rb'

    try:
        samfile = pysam.AlignmentFile(bam_file, mode)
    except Exception as e:
        print(f"Error opening alignment file: {e}", file=sys.stderr)
        sys.exit(1)

    # Open error log file if specified
    error_log = None
    if error_log_file:
        try:
            error_log = open(error_log_file, 'w')
            error_log.write("# Reads where reference sequence could not be retrieved\n")
            error_log.write("# Format: read_name\tchromosome\tposition\tCIGAR\tMD_tag\terror_reason\n")
        except Exception as e:
            print(f"Warning: Could not open error log file: {e}", file=sys.stderr)
            error_log = None

    # Write per-read CSV header
    if reads_out is not None:
        reads_out.write(','.join(_reads_csv_header_parts()) + '\n')

    for read in samfile:
        global_stats['num_alignments'] += 1

        # Skip unmapped, secondary, and supplementary reads
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        read_name = read.query_name
        tax_id = read.reference_name

        # Ensure taxid entry exists
        if tax_id not in stats_by_taxid:
            stats_by_taxid[tax_id] = _make_taxid_stats()

        # Check for duplicate read names
        if read_name in global_stats['reads_seen']:
            global_stats['duplicate_reads'].append(read_name)
        else:
            global_stats['reads_seen'].add(read_name)
            stats_by_taxid[tax_id]['num_reads'] += 1

        # Track reverse strand reads
        is_reverse = read.is_reverse
        if is_reverse:
            stats_by_taxid[tax_id]['num_reverse'] += 1

        # Get query sequence
        query_seq = read.query_sequence
        if query_seq is None:
            if error_log:
                md_tag = read.get_tag('MD') if read.has_tag('MD') else 'N/A'
                cigar = read.cigarstring if read.cigarstring else 'N/A'
                error_log.write(f"{read_name}\t{tax_id}\t{read.reference_start}\t{cigar}\t{md_tag}\tno_query_sequence\n")
            global_stats['failed_reads'] += 1
            continue

        # Get reference sequence for the aligned portion
        ref_seq = None
        error_reason = None
        try:
            ref_seq = read.get_reference_sequence()
        except Exception as e:
            error_reason = f"exception:{str(e)}"

        if ref_seq is None:
            if error_reason is None:
                error_reason = "returned_None"
            if error_log:
                md_tag = read.get_tag('MD') if read.has_tag('MD') else 'N/A'
                cigar = read.cigarstring if read.cigarstring else 'N/A'
                error_log.write(f"{read_name}\t{tax_id}\t{read.reference_start}\t{cigar}\t{md_tag}\t{error_reason}\n")
            global_stats['failed_reads'] += 1
            continue

        # Get aligned query sequence (excludes soft-clipped bases)
        query_aligned = read.query_alignment_sequence

        if query_aligned is None or len(query_aligned) != len(ref_seq):
            if error_log:
                md_tag = read.get_tag('MD') if read.has_tag('MD') else 'N/A'
                cigar = read.cigarstring if read.cigarstring else 'N/A'
                q_len = len(query_aligned) if query_aligned else 0
                r_len = len(ref_seq)
                error_log.write(f"{read_name}\t{tax_id}\t{read.reference_start}\t{cigar}\t{md_tag}\tlength_mismatch:query_aligned={q_len},ref={r_len}\n")
            global_stats['failed_reads'] += 1
            continue

        taxid_stats = stats_by_taxid[tax_id]

        # Per-read accumulators (written to reads_out before aggregating)
        read_subs = defaultdict(int)                        # (ref, obs) -> count
        read_ref_bases = defaultdict(int)                   # ref_base -> count (all aligned positions)

        # Compare sequences to find substitutions
        for i in range(len(ref_seq)):
            ref_base = ref_seq[i].upper()
            query_base = query_aligned[i].upper()

            # Apply strand correction if needed
            if strand_correction and is_reverse:
                ref_base = complement(ref_base)
                query_base = complement(query_base)

            # Count bases directly into taxid stats (all positions)
            taxid_stats['reference_bases'][ref_base] += 1
            taxid_stats['observed_bases'][query_base] += 1
            read_ref_bases[ref_base] += 1

            # Record all base pairs (including matches) per read
            read_subs[(ref_base, query_base)] += 1

        # Write per-read CSV row (substitutions only, no matches)
        if reads_out is not None:
            row = [read_name, tax_id,
                   str(read_ref_bases['A']), str(read_ref_bases['C']),
                   str(read_ref_bases['G']), str(read_ref_bases['T'])]
            for ref in BASES:
                for obs in BASES:
                    if ref != obs:
                        row.append(str(read_subs[(ref, obs)]))
            reads_out.write(','.join(row) + '\n')

        # Aggregate per-read counters into taxid stats
        for key, cnt in read_subs.items():
            taxid_stats['substitutions'][key] += cnt

    samfile.close()

    if error_log:
        error_log.close()
        print(f"Error log written to: {error_log_file}", file=sys.stderr)

    return global_stats, stats_by_taxid


def write_stats(global_stats, stats_by_taxid, out=sys.stdout):
    """Write formatted stats report with one section per tax ID."""

    # Global duplicate-read warning
    if global_stats['duplicate_reads']:
        print("WARNING: Found reads with multiple alignments!", file=out)
        print(f"Number of duplicate read names: {len(global_stats['duplicate_reads'])}", file=out)
        print("First few duplicates:", global_stats['duplicate_reads'][:5], file=out)
        print(file=out)
    else:
        print("All reads have only one alignment", file=out)
        print(file=out)

    for tax_id in sorted(stats_by_taxid.keys()):
        ts = stats_by_taxid[tax_id]

        print("=" * 60, file=out)
        print(f"TAX ID: {tax_id}", file=out)
        print("=" * 60, file=out)
        print(file=out)

        # Basic statistics
        print("=" * 60, file=out)
        print("BASIC STATISTICS", file=out)
        print("=" * 60, file=out)
        print(f"Total number of reads: {ts['num_reads']}", file=out)
        print(f"Total number of alignments: {global_stats['num_alignments']}", file=out)
        print(f"Reverse strand reads: {ts['num_reverse']}", file=out)
        print(f"Failed reads (no reference sequence): {global_stats['failed_reads']}", file=out)
        if global_stats['failed_reads'] > 0:
            print(f"  (See error log for details)", file=out)
        print(file=out)

        # Strand correction info
        if global_stats['strand_correction']:
            print("Strand correction: ENABLED (reverse reads converted to original orientation)", file=out)
        else:
            print("Strand correction: DISABLED (reads counted as-is)", file=out)
        print(file=out)

        # Observed base counts and frequencies
        print("=" * 60, file=out)
        print("OBSERVED BASE COUNTS (from reads)", file=out)
        print("=" * 60, file=out)
        total_obs = sum(ts['observed_bases'].values())
        print(f"{'Base':<10} {'Count':<15} {'Frequency':<15}", file=out)
        print("-" * 60, file=out)
        for base in ['A', 'C', 'G', 'T', 'N']:
            count = ts['observed_bases'][base]
            freq = count / total_obs if total_obs > 0 else 0
            print(f"{base:<10} {count:<15} {freq:<15.4f}", file=out)
        print(f"{'Total':<10} {total_obs:<15}", file=out)
        print(file=out)

        # Reference base counts and frequencies
        print("=" * 60, file=out)
        print("REFERENCE BASE COUNTS (from reference sequence)", file=out)
        print("=" * 60, file=out)
        total_ref = sum(ts['reference_bases'].values())
        print(f"{'Base':<10} {'Count':<15} {'Frequency':<15}", file=out)
        print("-" * 60, file=out)
        for base in ['A', 'C', 'G', 'T', 'N']:
            count = ts['reference_bases'][base]
            freq = count / total_ref if total_ref > 0 else 0
            print(f"{base:<10} {count:<15} {freq:<15.4f}", file=out)
        print(f"{'Total':<10} {total_ref:<15}", file=out)
        print(file=out)

        # Substitution table
        print("=" * 60, file=out)
        print("SUBSTITUTION TABLE (ordered by count)", file=out)
        print("=" * 60, file=out)

        substitution_list = []
        total_subs = 0
        for (ref_base, obs_base), count in ts['substitutions'].items():
            if ref_base != obs_base:
                ref_count = ts['reference_bases'][ref_base]
                freq = (count / ref_count * 100) if ref_count > 0 else 0
                substitution_list.append((f"{ref_base}>{obs_base}", count, freq))
                total_subs += count

        substitution_list.sort(key=lambda x: x[1], reverse=True)

        print(f"{'Substitution':<15} {'Count':>15} {'Frequency (%)':>20}", file=out)
        print("-" * 60, file=out)
        for sub_type, count, freq in substitution_list:
            print(f"{sub_type:<15} {count:>15} {freq:>19.4f}%", file=out)
        print("-" * 60, file=out)
        print(f"{'Total':<15} {total_subs:>15}", file=out)
        print(file=out)


def write_csv(global_stats, stats_by_taxid, out):
    """Write CSV summary with one row per tax ID."""
    basename = os.path.basename(global_stats['filename'])
    filename_no_ext = os.path.splitext(basename)[0]

    out.write(','.join(_csv_header_parts(include_taxid=True)) + '\n')

    for tax_id in sorted(stats_by_taxid.keys()):
        ts = stats_by_taxid[tax_id]
        total_obs = sum(ts['observed_bases'].values())
        total_subs = sum(cnt for (r, o), cnt in ts['substitutions'].items() if r != o)

        row = [filename_no_ext, tax_id, str(total_obs)]

        for base in BASES:
            row.append(str(ts['observed_bases'][base]))

        row.append(str(total_subs))

        for base in BASES:
            row.append(str(ts['reference_bases'][base]))

        for ref in BASES:
            for obs in BASES:
                if ref != obs:
                    row.append(str(ts['substitutions'][(ref, obs)]))

        out.write(','.join(row) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Analyze substitutions in BAM/SAM files by comparing reads to reference sequence',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s input.bam
  %(prog)s input.sam -o failed_reads.tsv
  %(prog)s aligned_reads.bam -o failed_reads.tsv --no-strand-correction
  %(prog)s input.bam -s stats_report.txt -c summary.csv -r per_read.csv

Note: Reference bases are reconstructed from the MD tag (no reference FASTA needed); add
      missing MD tags with: samtools calmd -b in.bam ref.fa > out.bam
      Automatically detects BAM or SAM format based on file extension.
      By default, applies strand correction to convert reverse reads to original orientation.
      The RNAME field (read.reference_name) is used as the tax ID for multi-taxa BAM files.
        """
    )
    parser.add_argument('alignment_file', help='Input BAM or SAM file')
    parser.add_argument('-o', '--error-log', dest='error_log',
                        help='Output file for reads where reference sequence could not be retrieved')
    parser.add_argument('--no-strand-correction', action='store_true',
                        help='Do NOT convert reverse reads to original orientation (default: apply correction)')
    parser.add_argument('-s', '--stats-output', dest='stats_output',
                        help='Path for human-readable stats report (default: stdout)')
    parser.add_argument('-c', '--csv-output', dest='csv_output',
                        help='Path for CSV summary file (default: stdout)')
    parser.add_argument('-r', '--reads-output', dest='reads_output',
                        help='Path for per-read CSV file (optional; skipped if not given)')

    args = parser.parse_args()

    print(f"Analyzing alignment file: {args.alignment_file}", file=sys.stderr)
    if args.no_strand_correction:
        print("Strand correction: DISABLED", file=sys.stderr)
    else:
        print("Strand correction: ENABLED (converting reverse reads to original orientation)", file=sys.stderr)
    print(file=sys.stderr)

    # Open output file handles
    stats_fh = open(args.stats_output, 'w') if args.stats_output else sys.stdout
    csv_fh = open(args.csv_output, 'w') if args.csv_output else sys.stdout
    reads_fh = open(args.reads_output, 'w') if args.reads_output else None

    try:
        global_stats, stats_by_taxid = analyze_bam(
            args.alignment_file,
            args.error_log,
            strand_correction=not args.no_strand_correction,
            reads_out=reads_fh,
        )
        write_stats(global_stats, stats_by_taxid, out=stats_fh)
        write_csv(global_stats, stats_by_taxid, out=csv_fh)
    finally:
        if args.stats_output:
            stats_fh.close()
        if args.csv_output:
            csv_fh.close()
        if reads_fh is not None:
            reads_fh.close()


if __name__ == '__main__':
    main()
