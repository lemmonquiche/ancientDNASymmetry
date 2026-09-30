#!/usr/bin/env python3
"""
Calculate substitution rates from BAM/SAM using MD tags (no reference FASTA needed).

Fixes two bugs in calc_noOrien.py:
  1. GTR parsing: reads the 6th rate parameter (GT) from its own line instead of
     silently using 1.0.
  2. Expected rates: uses data-derived base frequencies instead of model frequencies,
     so obs/exp is symmetric for complementary pairs when data composition is balanced.
"""

import sys
import argparse
import csv
from collections import defaultdict
import pysam


BASES = 'ACGT'

COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# Maps each substitution type to its index in the rate_matrix list [AC,AG,AT,CG,CT,GT]
SUB_TO_IDX = {
    'A→C': 0, 'C→A': 0,
    'A→G': 1, 'G→A': 1,
    'A→T': 2, 'T→A': 2,
    'C→G': 3, 'G→C': 3,
    'C→T': 4, 'T→C': 4,
    'G→T': 5, 'T→G': 5,
}


def complement(base):
    return COMP.get(base, 'N')


def parse_gtr_parameters(gtr_file_path):
    """
    Parse GTR parameter file.

    File format (4 non-empty lines):
        line 0: gamma shape parameter
        line 1: A C G T base frequencies
        line 2: AC AG AT CG CT  (5 exchangeability rates)
        line 3: GT              (6th exchangeability rate, own line)

    Returns dict with keys 'base_frequencies' and 'rate_matrix' (list of 6 floats,
    order: AC AG AT CG CT GT), or raises on parse error.
    """
    with open(gtr_file_path, 'r') as f:
        lines = f.readlines()

    base_freq_line = lines[1].strip().split()
    base_frequencies = {
        'A': float(base_freq_line[0]),
        'C': float(base_freq_line[1]),
        'G': float(base_freq_line[2]),
        'T': float(base_freq_line[3]),
    }

    rate_line = lines[2].strip().split()
    rate_matrix = [float(x) for x in rate_line]  # 5 values

    sixth_param = float(lines[3].strip())
    rate_matrix.append(sixth_param)  # now 6 values: AC AG AT CG CT GT

    return {
        'base_frequencies': base_frequencies,
        'rate_matrix': rate_matrix,
    }


def get_expected_gtr_rate_from_data(sub_type, gtr_params, ref_base_counts):
    """
    Expected GTR rate for sub_type using data-derived target-base frequency.

    rate(i→j) = rate_matrix[ij] × data_freq[j]

    This is symmetric for complementary pairs (e.g. A→C and C→A share the same
    exchangeability parameter; when data is balanced the obs/exp ratios match).
    """
    total_bases = sum(ref_base_counts.values())
    if total_bases == 0:
        return None

    data_freq = {b: ref_base_counts.get(b, 0) / total_bases for b in BASES}

    matrix_idx = SUB_TO_IDX.get(sub_type)
    if matrix_idx is None:
        return None

    target_base = sub_type.split('→')[1]
    return gtr_params['rate_matrix'][matrix_idx] * data_freq[target_base]


def compute_norm_factor(obs_counts, ref_counts, gtr_params):
    """
    Compute k such that sum(exp_count * k) == sum(obs_count) over all 12 substitution types.
    Returns k, or 1.0 if gtr_params is None or total raw expected is zero.
    """
    if gtr_params is None:
        return 1.0

    total_obs = 0
    total_exp_raw = 0.0

    for ref in BASES:
        for obs in BASES:
            if ref == obs:
                continue
            sub_type = f"{ref}→{obs}"
            total_obs += obs_counts[ref][obs]
            ref_count = ref_counts[ref]
            exp_rate = get_expected_gtr_rate_from_data(sub_type, gtr_params, ref_counts)
            if exp_rate is not None and ref_count > 0:
                total_exp_raw += exp_rate * ref_count

    if total_exp_raw == 0:
        return 1.0
    return total_obs / total_exp_raw


def analyze_alignment(input_file, strand_specific=True, min_baseq=0, min_mapq=0, max_pos=0):
    """
    Count reference bases and substitutions from a SAM/BAM file using MD tags.

    Args:
        max_pos: When > 0, also accumulate per-position counts for first/last max_pos bases.

    Returns:
        ref_counts           – dict[base] -> int
        obs_counts           – dict[ref_base][obs_base] -> int
        position_ref_counts  – dict[pos][base] -> int, or None when max_pos == 0
        position_sub_counts  – dict[pos][(ref,obs)] -> int, or None when max_pos == 0
    """
    ref_counts = defaultdict(int)
    obs_counts = defaultdict(lambda: defaultdict(int))

    if max_pos > 0:
        position_ref_counts = defaultdict(lambda: defaultdict(int))
        position_sub_counts = defaultdict(lambda: defaultdict(int))
    else:
        position_ref_counts = None
        position_sub_counts = None

    # "r" auto-detects SAM vs BAM; check_sq=False allows files without @SQ headers
    aln = pysam.AlignmentFile(input_file, "r", check_sq=False)

    reads_processed = 0
    reads_reverse = 0
    bases_processed = 0
    md_error_shown = False

    print("Processing alignment file...", file=sys.stderr)
    if strand_specific:
        print("  Orientation: converting reverse reads to original strand", file=sys.stderr)
    else:
        print("  Orientation: using alignment strand as-is", file=sys.stderr)
    print(f"  Filters: min_baseq={min_baseq}, min_mapq={min_mapq}", file=sys.stderr)

    for read in aln:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue

        reads_processed += 1
        if reads_processed % 100000 == 0:
            print(f"  {reads_processed} reads, {bases_processed} bases...", file=sys.stderr)

        is_reverse = read.is_reverse
        if is_reverse:
            reads_reverse += 1

        qual_array = read.query_qualities
        read_seq = read.query_sequence
        if read_seq is None:
            continue
        read_seq = read_seq.upper()

        read_length = read.query_length

        # get_aligned_pairs with_seq=True uses the MD tag to provide ref bases
        for query_pos, ref_pos, ref_base in read.get_aligned_pairs(
                matches_only=True, with_seq=True):

            if ref_base is None:
                if not md_error_shown:
                    print(
                        "\nERROR: MD tag is missing or incomplete — ref_base is None.\n"
                        "Please add MD tags with:\n"
                        "  samtools calmd -b <in.bam> <ref.fa> > <out.bam>\n"
                        "then re-run this script.",
                        file=sys.stderr,
                    )
                    aln.close()
                    sys.exit(1)
                continue  # unreachable after exit, kept for clarity

            # pysam convention: lowercase = match, uppercase = mismatch
            ref_base = ref_base.upper()

            read_base = read_seq[query_pos]

            if ref_base not in BASES or read_base not in BASES:
                continue

            if qual_array is not None and qual_array[query_pos] < min_baseq:
                continue

            if strand_specific and is_reverse:
                ref_base = complement(ref_base)
                read_base = complement(read_base)

            ref_counts[ref_base] += 1
            obs_counts[ref_base][read_base] += 1
            bases_processed += 1

            if max_pos > 0:
                if strand_specific and is_reverse:
                    pos_from_start = read_length - query_pos
                    pos_from_end = -query_pos - 1
                else:
                    pos_from_start = query_pos + 1
                    pos_from_end = query_pos - read_length

                if pos_from_start <= max_pos:
                    position_ref_counts[pos_from_start][ref_base] += 1
                    position_sub_counts[pos_from_start][(ref_base, read_base)] += 1

                if pos_from_end >= -max_pos:
                    position_ref_counts[pos_from_end][ref_base] += 1
                    position_sub_counts[pos_from_end][(ref_base, read_base)] += 1

    aln.close()
    print(
        f"\nProcessed {reads_processed} reads ({reads_reverse} reverse), "
        f"{bases_processed} bases",
        file=sys.stderr,
    )
    return ref_counts, obs_counts, position_ref_counts, position_sub_counts


def print_results(ref_counts, obs_counts, gtr_params=None, output_file=None, norm_factor=1.0):
    """Print analysis results to stdout (and optionally a TSV file)."""

    total_ref = sum(ref_counts[b] for b in BASES)

    # ------------------------------------------------------------------
    # 1. Reference base coverage
    # ------------------------------------------------------------------
    print("=" * 70)
    print("REFERENCE BASE COVERAGE")
    print("=" * 70)
    print(f"{'Base':<10} {'Count':>15} {'Proportion':>12}")
    print("-" * 40)
    for base in BASES:
        count = ref_counts[base]
        prop = count / total_ref if total_ref > 0 else 0.0
        print(f"{base:<10} {count:>15,} {prop:>12.4f}")
    print(f"{'Total':<10} {total_ref:>15,}")
    print()

    # ------------------------------------------------------------------
    # 2. Substitution matrix
    # ------------------------------------------------------------------
    print("=" * 70)
    print("SUBSTITUTION MATRIX (Reference → Observed)")
    print("=" * 70)
    header = f"{'Ref\\Obs':<8}" + "".join(f"{b:>12}" for b in BASES)
    header += f"{'Errors':>12}{'ErrorRate':>12}"
    print(header)
    print("-" * (8 + 4 * 12 + 12 + 12))

    total_matches = 0
    total_errors = 0

    for ref in BASES:
        row = f"{ref:<8}"
        for obs in BASES:
            row += f"{obs_counts[ref][obs]:>12,}"
        matches = obs_counts[ref][ref]
        errors = sum(obs_counts[ref][obs] for obs in BASES if obs != ref)
        total_matches += matches
        total_errors += errors
        err_rate = errors / ref_counts[ref] if ref_counts[ref] > 0 else 0.0
        row += f"{errors:>12,}{err_rate:>12.6f}"
        print(row)

    overall_err = total_errors / (total_matches + total_errors) if (total_matches + total_errors) > 0 else 0.0
    print("-" * (8 + 4 * 12 + 12 + 12))
    print(f"{'TOTAL':<8}{'':>{4*12}}{total_errors:>12,}{overall_err:>12.6f}")
    print()

    # ------------------------------------------------------------------
    # 3. Per-substitution table
    # ------------------------------------------------------------------
    print("=" * 70)
    print("PER-SUBSTITUTION RATES")
    print("=" * 70)

    has_gtr = gtr_params is not None

    if has_gtr:
        print(f"{'Sub':<8} {'ObsCount':>10} {'RefCount':>10} {'ObsRate':>12}"
              f" {'ExpRate':>12} {'ExpCount':>12} {'Obs/Exp':>10} {'piQ':>12}")
        print("-" * 91)
    else:
        print(f"{'Sub':<8} {'ObsCount':>10} {'RefCount':>10} {'ObsRate':>12} {'piQ':>12}")
        print("-" * 57)

    # First pass: compute all per-type values before printing, so we can find the
    # minimum obs/exp ratio and rescale so that the smallest ratio equals 1.
    raw_rows = []
    for ref in BASES:
        for obs in BASES:
            if ref == obs:
                continue
            sub_type = f"{ref}→{obs}"
            obs_count = obs_counts[ref][obs]
            ref_count = ref_counts[ref]
            obs_rate = obs_count / ref_count if ref_count > 0 else 0.0
            pi_q = obs_count / total_ref if total_ref > 0 else 0.0

            if has_gtr:
                exp_rate = get_expected_gtr_rate_from_data(sub_type, gtr_params, ref_counts)
                if exp_rate is not None:
                    exp_rate = exp_rate * norm_factor
                exp_count = exp_rate * ref_count if (exp_rate is not None and ref_count > 0) else None
                obs_exp = obs_count / exp_count if (exp_count is not None and exp_count > 0) else None
            else:
                exp_rate = None
                exp_count = None
                obs_exp = None

            raw_rows.append((sub_type, obs_count, ref_count, obs_rate,
                              exp_rate, exp_count, obs_exp, pi_q))

    # Scale ratios so the minimum is 1.
    if has_gtr:
        valid_ratios = [r[6] for r in raw_rows if r[6] is not None and r[6] > 0]
        ratio_min = min(valid_ratios) if valid_ratios else 1.0
    else:
        ratio_min = 1.0

    tsv_rows = []
    for sub_type, obs_count, ref_count, obs_rate, exp_rate, exp_count, obs_exp, pi_q in raw_rows:
        if has_gtr:
            scaled_obs_exp = obs_exp / ratio_min if obs_exp is not None else None

            exp_rate_s  = f"{exp_rate:.8f}"       if exp_rate       is not None else "NA"
            exp_count_s = f"{exp_count:.2f}"       if exp_count      is not None else "NA"
            obs_exp_s   = f"{scaled_obs_exp:.4f}"  if scaled_obs_exp is not None else "NA"

            print(f"{sub_type:<8} {obs_count:>10,} {ref_count:>10,} {obs_rate:>12.8f}"
                  f" {exp_rate_s:>12} {exp_count_s:>12} {obs_exp_s:>10} {pi_q:>12.8f}")

            tsv_rows.append((sub_type, obs_count, ref_count, obs_rate,
                             exp_rate_s, exp_count_s, obs_exp_s, pi_q))
        else:
            print(f"{sub_type:<8} {obs_count:>10,} {ref_count:>10,} {obs_rate:>12.8f}"
                  f" {pi_q:>12.8f}")
            tsv_rows.append((sub_type, obs_count, ref_count, obs_rate,
                             "NA", "NA", "NA", pi_q))

    print()

    # ------------------------------------------------------------------
    # 4. Summary
    # ------------------------------------------------------------------
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Overall error rate : {overall_err:.6f} ({overall_err * 100:.4f}%)")
    print(f"Total matches      : {total_matches:,}")
    print(f"Total mismatches   : {total_errors:,}")
    print("=" * 70)

    # ------------------------------------------------------------------
    # Optional TSV output
    # ------------------------------------------------------------------
    if output_file:
        with open(output_file, 'w') as fh:
            fh.write("sub_type\tobs_count\tref_count\tobs_rate\t"
                     "expected_rate\texpected_count\tobs_exp_ratio\tpi_q\n")
            for row in tsv_rows:
                sub_type, obs_count, ref_count, obs_rate, er, ec, oe, pi_q = row
                fh.write(f"{sub_type}\t{obs_count}\t{ref_count}\t{obs_rate:.8f}\t"
                         f"{er}\t{ec}\t{oe}\t{pi_q:.8f}\n")
        print(f"\nTSV written to: {output_file}", file=sys.stderr)

    return ratio_min


def write_positional_csv(position_ref_counts, position_sub_counts, gtr_params,
                         pos_output, bases_output=None, max_pos=15, norm_factor=1.0,
                         ratio_scale=1.0):
    """
    Write positional substitution CSV with same column format as calc_pos_mine.py.

    Uses calc_subs.py's GTR dict format (rate_matrix list + SUB_TO_IDX) and
    MD-tag-derived bases rather than FASTA.
    When gtr_params is None, expected_rate/expected_count are 0, obs_exp_ratio is inf/0.
    """
    all_positions = list(range(1, max_pos + 1)) + list(range(-max_pos, 0))

    bases_fh = open(bases_output, 'w', newline='') if bases_output else None
    base_writer = csv.writer(bases_fh) if bases_fh else None

    try:
        with open(pos_output, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                'position', 'position_type', 'substitution',
                'ref_base', 'obs_base',
                'observed_count', 'ref_count', 'observed_rate',
                'expected_rate', 'expected_count', 'obs_exp_ratio', 'pi_q',
            ])

            if base_writer:
                base_writer.writerow([
                    'position', 'ref_A', 'ref_C', 'ref_G', 'ref_T', 'ref_total',
                    'obs_A', 'obs_C', 'obs_G', 'obs_T', 'obs_total',
                ])

            for pos in all_positions:
                if pos not in position_ref_counts:
                    continue

                pos_type = '5prime' if pos > 0 else '3prime'
                ref_pos_total = sum(position_ref_counts[pos].values())

                # Per-position local base frequencies and expected rates (GTR)
                position_expected_rate = {}
                local_freq = {b: 0.0 for b in BASES}
                if ref_pos_total > 0:
                    local_freq = {b: position_ref_counts[pos].get(b, 0) / ref_pos_total
                                  for b in BASES}

                if gtr_params:
                    for ref in BASES:
                        if position_ref_counts[pos].get(ref, 0) == 0:
                            continue
                        for obs in BASES:
                            if obs != ref:
                                idx = SUB_TO_IDX[f"{ref}→{obs}"]
                                position_expected_rate[(ref, obs)] = (
                                    gtr_params['rate_matrix'][idx] * local_freq[obs] * norm_factor
                                )

                # Base frequency row
                if base_writer:
                    total_obs = sum(position_sub_counts[pos].values())
                    pos_obs_counts = defaultdict(int)
                    for (r, o), val in position_sub_counts[pos].items():
                        pos_obs_counts[o] += val
                    obs_freq = ({b: pos_obs_counts[b] / total_obs for b in BASES}
                                if total_obs > 0 else {b: 0.0 for b in BASES})
                    base_writer.writerow([
                        pos,
                        local_freq['A'], local_freq['C'], local_freq['G'], local_freq['T'],
                        ref_pos_total,
                        obs_freq['A'], obs_freq['C'], obs_freq['G'], obs_freq['T'],
                        total_obs,
                    ])

                # Substitution rows
                for ref in BASES:
                    total_ref = position_ref_counts[pos].get(ref, 0)
                    if total_ref == 0:
                        continue
                    for obs in BASES:
                        if obs == ref:
                            continue

                        obs_count = position_sub_counts[pos].get((ref, obs), 0)
                        obs_rate = obs_count / total_ref if total_ref > 0 else 0.0

                        if gtr_params:
                            exp_rate = position_expected_rate.get((ref, obs), 0.0)
                            exp_count = exp_rate * total_ref
                        else:
                            exp_rate = 0.0
                            exp_count = 0.0

                        pi_q = (total_ref / ref_pos_total) * obs_rate

                        if exp_rate > 0:
                            obs_exp = (obs_rate / exp_rate) / ratio_scale
                        else:
                            obs_exp = float('inf') if obs_rate > 0 else 0.0

                        writer.writerow([
                            pos, pos_type, f"{ref}>{obs}",
                            ref, obs,
                            obs_count, total_ref, f"{obs_rate:.8f}",
                            f"{exp_rate:.8f}", f"{exp_count:.6f}",
                            f"{obs_exp:.6f}", f"{pi_q:.6f}",
                        ])
    finally:
        if bases_fh:
            bases_fh.close()


def main():
    parser = argparse.ArgumentParser(
        description=(
            'Substitution analysis from SAM/BAM using MD tags (no reference FASTA). '
            'Fixes GTR parsing and expected-rate formula vs calc_noOrien.py.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('input_file', help='Input SAM or BAM file')
    parser.add_argument('--gtr-params', metavar='FILE',
                        help='GTR parameter file for expected substitution rates')
    parser.add_argument('--no-strand-correction', action='store_true',
                        help='Skip complementing reverse-strand bases')
    parser.add_argument('--min-baseq', type=int, default=0,
                        metavar='INT', help='Minimum base quality (default: 0)')
    parser.add_argument('--min-mapq', type=int, default=0,
                        metavar='INT', help='Minimum mapping quality (default: 0)')
    parser.add_argument('--output', metavar='FILE',
                        help='Write results to TSV file in addition to stdout')
    parser.add_argument('--pos-output', metavar='FILE',
                        help='Positional substitution CSV (same columns as calc_pos_mine.py -o)')
    parser.add_argument('--pos-bases', metavar='FILE',
                        help='Positional base frequency CSV (same columns as calc_pos_mine.py -b)')
    parser.add_argument('--max-pos', type=int, default=15, metavar='INT',
                        help='Positions from each end to analyze (default: 15)')

    args = parser.parse_args()

    # Load GTR parameters
    gtr_params = None
    if args.gtr_params:
        try:
            gtr_params = parse_gtr_parameters(args.gtr_params)
            bf = gtr_params['base_frequencies']
            rm = gtr_params['rate_matrix']
            print(f"Loaded GTR parameters from {args.gtr_params}", file=sys.stderr)
            print(f"  Base frequencies: A={bf['A']:.6f} C={bf['C']:.6f} "
                  f"G={bf['G']:.6f} T={bf['T']:.6f}", file=sys.stderr)
            print(f"  Rate matrix (AC AG AT CG CT GT): "
                  + " ".join(f"{r:.6f}" for r in rm), file=sys.stderr)
        except Exception as exc:
            print(f"Error loading GTR parameters: {exc}", file=sys.stderr)
            sys.exit(1)

    max_pos = args.max_pos if args.pos_output else 0

    ref_counts, obs_counts, position_ref_counts, position_sub_counts = analyze_alignment(
        args.input_file,
        strand_specific=not args.no_strand_correction,
        min_baseq=args.min_baseq,
        min_mapq=args.min_mapq,
        max_pos=max_pos,
    )

    norm_factor = compute_norm_factor(obs_counts, ref_counts, gtr_params)
    if gtr_params:
        print(f"  Normalization factor (k): {norm_factor:.6f}", file=sys.stderr)

    ratio_min = print_results(ref_counts, obs_counts, gtr_params=gtr_params,
                              output_file=args.output, norm_factor=norm_factor)

    if args.pos_output:
        write_positional_csv(
            position_ref_counts, position_sub_counts, gtr_params,
            args.pos_output,
            bases_output=args.pos_bases,
            max_pos=args.max_pos,
            norm_factor=norm_factor,
            ratio_scale=ratio_min,
        )
        print(f"\nWrote positional analysis to: {args.pos_output}", file=sys.stderr)
        if args.pos_bases:
            print(f"Wrote positional base frequencies to: {args.pos_bases}", file=sys.stderr)


if __name__ == '__main__':
    main()
