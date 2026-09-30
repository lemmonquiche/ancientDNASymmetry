#!/usr/bin/env python3
"""
Calculate substitution rates from BAM/SAM using MD tags (no reference FASTA needed),
scoring expected substitutions under a 12-parameter UNREST (non-reversible) model.

This is the UNREST counterpart of calc_subs.py. Where calc_subs.py uses a GTR model
(6 symmetric exchangeabilities × target-base frequency), this script uses the 12
independent directional rates from an UNREST fit (e.g. IQ-TREE UNREST+FO).

Key difference from GTR:
  - GTR off-diagonal rate:    Q[i→j] = exchangeability_ij × pi_j   (target freq factor)
  - UNREST off-diagonal rate: Q[i→j] = c × R[i→j]                  (directed param, NO
    target-freq factor; verified against IQ-TREE's reported Q matrix)

Expected counts therefore follow:
    exp_count(i→j) = pi_i · Q_ij · N = ref_count[i] × R[i→j] × (norm)
where pi_i (the SOURCE base frequency) is carried by ref_count[i], R[i→j] is the directed
UNREST parameter, and the global constant c is absorbed by the normalization factor.
"""

import sys
import argparse
import csv
from collections import defaultdict
import pysam


BASES = 'ACGT'

COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# The 12 directed substitution types under UNREST (no shared parameters).
SUB_TYPES = [f"{i}→{j}" for i in BASES for j in BASES if i != j]


def complement(base):
    return COMP.get(base, 'N')


def parse_unrest_parameters(unrest_file_path):
    """
    Parse an UNREST parameter file.

    File format (12 non-empty lines, one per directed rate):
        A-C: 0.997
        A-G: 2.274
        ...
        T-G: 1.000
    where each line is "SOURCE-TARGET: rate". No base-frequency lines are expected;
    UNREST equilibrium frequencies are a function of Q, not supplied here.

    Returns dict with key 'rates' -> dict[(src, tgt)] = float, covering all 12 directed
    base pairs. Raises on parse error or if any of the 12 pairs is missing.
    """
    rates = {}
    with open(unrest_file_path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            # Expect "X-Y: value"
            try:
                pair, value = line.split(':')
                src, tgt = pair.strip().split('-')
            except ValueError:
                raise ValueError(f"Cannot parse UNREST line: {raw!r} "
                                 "(expected 'SOURCE-TARGET: rate')")
            src = src.strip().upper()
            tgt = tgt.strip().upper()
            if src not in BASES or tgt not in BASES or src == tgt:
                raise ValueError(f"Invalid base pair in UNREST line: {raw!r}")
            rates[(src, tgt)] = float(value.strip())

    expected_pairs = {(i, j) for i in BASES for j in BASES if i != j}
    missing = expected_pairs - set(rates)
    if missing:
        missing_str = ", ".join(f"{s}-{t}" for s, t in sorted(missing))
        raise ValueError(f"UNREST file missing directed rate(s): {missing_str}")

    return {'rates': rates}


def get_expected_unrest_rate(sub_type, unrest_params):
    """
    Expected per-source-base UNREST rate for sub_type i→j.

    Under UNREST the off-diagonal rate matrix entry is the directed parameter itself
    (Q_ij = c × R_ij), so the per-source-base rate is just R[i→j]. There is NO
    multiplication by the target-base frequency (that is GTR-only). The source-base
    frequency enters later as ref_count[i] in the expected-count product, and the global
    constant c is absorbed by the normalization factor.
    """
    src, tgt = sub_type.split('→')
    return unrest_params['rates'].get((src, tgt))


def compute_norm_factor(obs_counts, ref_counts, unrest_params):
    """
    Compute k such that sum(exp_count * k) == sum(obs_count) over all 12 substitution types.
    Returns k, or 1.0 if unrest_params is None or total raw expected is zero.
    """
    if unrest_params is None:
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
            exp_rate = get_expected_unrest_rate(sub_type, unrest_params)
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


def print_results(ref_counts, obs_counts, unrest_params=None, output_file=None, norm_factor=1.0):
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
    ref_obs_label = 'Ref\\Obs'
    header = f"{ref_obs_label:<8}" + "".join(f"{b:>12}" for b in BASES)
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

    has_unrest = unrest_params is not None

    if has_unrest:
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

            if has_unrest:
                exp_rate = get_expected_unrest_rate(sub_type, unrest_params)
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
    if has_unrest:
        valid_ratios = [r[6] for r in raw_rows if r[6] is not None and r[6] > 0]
        ratio_min = min(valid_ratios) if valid_ratios else 1.0
    else:
        ratio_min = 1.0

    tsv_rows = []
    for sub_type, obs_count, ref_count, obs_rate, exp_rate, exp_count, obs_exp, pi_q in raw_rows:
        if has_unrest:
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


def write_positional_csv(position_ref_counts, position_sub_counts, unrest_params,
                         pos_output, bases_output=None, max_pos=15, norm_factor=1.0,
                         ratio_scale=1.0):
    """
    Write positional substitution CSV with same column format as calc_pos_mine.py.

    Uses the UNREST directed rates: expected_rate(i→j) = R[i→j] × norm_factor (no
    target-frequency factor). MD-tag-derived bases are used rather than a FASTA.
    When unrest_params is None, expected_rate/expected_count are 0, obs_exp_ratio is inf/0.
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

                # Per-position local base frequencies (used for the base-frequency row;
                # the UNREST expected rate itself does NOT use target frequency).
                local_freq = {b: 0.0 for b in BASES}
                if ref_pos_total > 0:
                    local_freq = {b: position_ref_counts[pos].get(b, 0) / ref_pos_total
                                  for b in BASES}

                # Per-position expected rates (UNREST): directed param × norm_factor
                position_expected_rate = {}
                if unrest_params:
                    for ref in BASES:
                        if position_ref_counts[pos].get(ref, 0) == 0:
                            continue
                        for obs in BASES:
                            if obs != ref:
                                position_expected_rate[(ref, obs)] = (
                                    unrest_params['rates'][(ref, obs)] * norm_factor
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

                        if unrest_params:
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
            'Substitution analysis from SAM/BAM using MD tags (no reference FASTA), '
            'scoring expected substitutions under a 12-parameter UNREST model. '
            'UNREST counterpart of calc_subs.py.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('input_file', help='Input SAM or BAM file')
    parser.add_argument('--unrest-params', metavar='FILE',
                        help='UNREST parameter file (12 lines "SOURCE-TARGET: rate") '
                             'for expected substitution rates')
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

    # Load UNREST parameters
    unrest_params = None
    if args.unrest_params:
        try:
            unrest_params = parse_unrest_parameters(args.unrest_params)
            rates = unrest_params['rates']
            print(f"Loaded UNREST parameters from {args.unrest_params}", file=sys.stderr)
            print("  Directed rates R[i→j]:", file=sys.stderr)
            for src in BASES:
                row = "    " + "  ".join(
                    f"{src}→{tgt}={rates[(src, tgt)]:.4f}"
                    for tgt in BASES if tgt != src)
                print(row, file=sys.stderr)
        except Exception as exc:
            print(f"Error loading UNREST parameters: {exc}", file=sys.stderr)
            sys.exit(1)

    max_pos = args.max_pos if args.pos_output else 0

    ref_counts, obs_counts, position_ref_counts, position_sub_counts = analyze_alignment(
        args.input_file,
        strand_specific=not args.no_strand_correction,
        min_baseq=args.min_baseq,
        min_mapq=args.min_mapq,
        max_pos=max_pos,
    )

    norm_factor = compute_norm_factor(obs_counts, ref_counts, unrest_params)
    if unrest_params:
        print(f"  Normalization factor (k): {norm_factor:.6f}", file=sys.stderr)

    ratio_min = print_results(ref_counts, obs_counts, unrest_params=unrest_params,
                              output_file=args.output, norm_factor=norm_factor)

    if args.pos_output:
        write_positional_csv(
            position_ref_counts, position_sub_counts, unrest_params,
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
