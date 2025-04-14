import primer3 # type: ignore
import requests
import time
from bs4 import BeautifulSoup

def bisulfite_convert(sequence: str):
    """
    Perform in silico bisulfite conversion.
    
    Args:
        sequence (str): Original DNA sequence (assumed uppercase).
        
    Returns:
        Tuple[str, str]: (methylated_sequence, unmethylated_sequence)
    """
    sequence = sequence.upper()
    methylated = []
    i = 0
    while i < len(sequence):
        if sequence[i] == 'C':
            # If this C is part of a CpG site
            if i + 1 < len(sequence) and sequence[i + 1] == 'G':
                methylated.append('C')  # leave C unchanged
            else:
                methylated.append('T')  # convert non-CpG C to T
        else:
            methylated.append(sequence[i])
        i += 1

    # Unmethylated: convert all C to T
    unmethylated = sequence.replace('C', 'T')
    
    return ''.join(methylated), unmethylated

def count_cpgs(seq):
    return sum(1 for i in range(len(seq) - 1) if seq[i:i+2] == 'CG')

def reverse_complement(seq):
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]

def find_cpg_excluded_regions(seq):
    ranges = []
    for i in range(len(seq) - 1):
        if seq[i] == 'C' and seq[i + 1] == 'G':
            ranges.append((i, i + 1))

    # Merge overlapping or adjacent regions
    merged = []
    for start, end in sorted(ranges):
        if not merged:
            merged.append([start, end])
        else:
            last_start, last_end = merged[-1]
            if start <= last_end + 1:
                merged[-1][1] = max(last_end, end)
            else:
                merged.append([start, end])

    # Convert to [start, length]
    result = []
    for start, end in merged:
        length = end - start + 1
        if start >= 0 and (start + length) <= len(seq):
            result.append([start, length])
    return result

def highlight_cpg(seq):
    """Highlight CpG dinucleotides in red bold font."""
    highlighted = ''
    i = 0
    while i < len(seq) - 1:
        if seq[i:i+2] == 'CG':
            highlighted += '<span style="color:red;"><strong>CG</strong></span>'
            i += 2
        else:
            highlighted += seq[i]
            i += 1
    if i < len(seq):
        highlighted += seq[i]
    return highlighted

def highlight_original_cpgs_in_unmeth_probe(original_seq, unmeth_probe_seq, probe_start):
    highlighted = ''
    i = 0
    while i < len(unmeth_probe_seq) - 1:
        orig_dinuc = original_seq[probe_start + i : probe_start + i + 2]
        probe_dinuc = unmeth_probe_seq[i:i+2]

        if orig_dinuc == 'CG':
            # Highlight both T and G in the TG
            highlighted += '<span style="color:red;"><strong>' + probe_dinuc + '</strong></span>'
            i += 2
        else:
            highlighted += unmeth_probe_seq[i]
            i += 1

    # If one base remains at the end, add it
    if i == len(unmeth_probe_seq) - 1:
        highlighted += unmeth_probe_seq[-1]

    return highlighted

def degenerate_primer(seq, strand='forward'):
    result = ''
    i = 0
    while i < len(seq) - 1:
        c, g = seq[i], seq[i + 1]
        if strand == 'forward' and c == 'C' and g == 'G':
            result += '<span style="color:blue;"><strong>Y</strong></span>'
            i += 1  # move to 'G'
        elif strand == 'reverse' and c == 'C' and g == 'G':
            result += 'C' + '<span style="color:blue;"><strong>R</strong></span>'
            i += 2  # skip both
        else:
            result += c
            i += 1
    if i == len(seq) - 1:
        result += seq[-1]
    return result

def is_probe_valid(seq):
    # Filter out probe if it starts with G or contains 5+ consecutive Gs
    if seq.startswith('G'):
        return False
    if 'GGGGG' in seq:
        return False
    return True

def design_probe_from_fixed_primers(seq, left_primer, right_primer_rc, config=None):
    
    config = config or {}
    
    result = primer3.bindings.design_primers(
        {
            'SEQUENCE_ID': 'probe_check',
            'SEQUENCE_TEMPLATE': seq,
            'SEQUENCE_PRIMER': left_primer,
            'SEQUENCE_PRIMER_REVCOMP': right_primer_rc
        },
        {
            'PRIMER_PICK_LEFT_PRIMER': 0,
            'PRIMER_PICK_RIGHT_PRIMER': 0,
            'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_INTERNAL_OPT_SIZE': 25,
            'PRIMER_INTERNAL_MIN_SIZE': 18,
            'PRIMER_INTERNAL_MAX_SIZE': 30,
            'PRIMER_INTERNAL_OPT_TM': config.get('probe_opt_tm', 60.0),
            'PRIMER_INTERNAL_MIN_TM': config.get('probe_min_tm', 57.0),
            'PRIMER_INTERNAL_MAX_TM': config.get('probe_max_tm', 63.0),
            'PRIMER_INTERNAL_OPT_GC_PERCENT': 50.0,
            'PRIMER_INTERNAL_MIN_GC': 20.0,
            'PRIMER_INTERNAL_MAX_GC': 80.0
        }
    )

    if result.get('PRIMER_INTERNAL_NUM_RETURNED', 0) > 0:
        probe_seq = result['PRIMER_INTERNAL_0_SEQUENCE']
        probe_start, _ = result['PRIMER_INTERNAL_0']
        return {
            'sequence': probe_seq,
            'tm': result['PRIMER_INTERNAL_0_TM'],
            'gc': result['PRIMER_INTERNAL_0_GC_PERCENT'],
            'probe_self_end': result.get('PRIMER_INTERNAL_0_SELF_END_TH'),
            'probe_hairpin': result.get('PRIMER_INTERNAL_0_HAIRPIN_TH'),
            'cpgs': count_cpgs(probe_seq),
            'start': probe_start 
        }

    return None

def run_primer3(seq, excluded=None, unmethylated_seq=None, original_seq=None, config=None):
    
    excluded = excluded or []
    config = config or {}

    result = primer3.bindings.design_primers(
        {
            'SEQUENCE_ID': 'test_sequence',
            'SEQUENCE_TEMPLATE': seq,
            'SEQUENCE_EXCLUDED_REGION': excluded
        },
        {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 30,
            'PRIMER_OPT_TM': config.get('primer_opt_tm', 55.0),
            'PRIMER_MIN_TM': config.get('primer_min_tm', 52.0),
            'PRIMER_MAX_TM': config.get('primer_max_tm', 58.0),
            'PRIMER_MAX_DIFF_TM': 2.0,
            'PRIMER_OPT_GC_PERCENT': 50.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_PRODUCT_SIZE_RANGE': config.get('product_size_range', [[70, 150]]),
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_SALT_MONOVALENT': config.get('salt_mono', 50.0),
            'PRIMER_SALT_DIVALENT': config.get('salt_div', 3.0),
            'PRIMER_DNTP_CONC': config.get('dntp_conc', 0.8),
            'PRIMER_NUM_RETURN': 10,
            'PRIMER_MIN_THREE_PRIME_DISTANCE': 3,
            'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_INTERNAL_OPT_SIZE': 25,
            'PRIMER_INTERNAL_MIN_SIZE': 18,
            'PRIMER_INTERNAL_MAX_SIZE': 30,
            'PRIMER_INTERNAL_OPT_TM': config.get('probe_opt_tm', 60.0),
            'PRIMER_INTERNAL_MIN_TM': config.get('probe_min_tm', 57.0),
            'PRIMER_INTERNAL_MAX_TM': config.get('probe_max_tm', 63.0),
            'PRIMER_INTERNAL_OPT_GC_PERCENT': 50.0,
            'PRIMER_INTERNAL_MIN_GC': 20.0,
            'PRIMER_INTERNAL_MAX_GC': 80.0
        }

    )

    primers = []
    for i in range(result['PRIMER_PAIR_NUM_RETURNED']):
        
        left_seq = result[f'PRIMER_LEFT_{i}_SEQUENCE']
        right_seq = result[f'PRIMER_RIGHT_{i}_SEQUENCE']
        left_pos, left_len = result[f'PRIMER_LEFT_{i}']
        right_pos, right_len = result[f'PRIMER_RIGHT_{i}']

        internal_seq = result.get(f'PRIMER_INTERNAL_{i}_SEQUENCE', '')
        internal_cpgs = count_cpgs(internal_seq)
        probe_m_tm = result.get(f'PRIMER_INTERNAL_{i}_TM')
        probe_m_gc = result.get(f'PRIMER_INTERNAL_{i}_GC_PERCENT')
        
        # Unmethylated probe logic (match positions)
        probe_u_info = None
        if unmethylated_seq:
            unmeth_left = unmethylated_seq[left_pos:left_pos + left_len]
            unmeth_right = unmethylated_seq[right_pos:right_pos + right_len]
            unmeth_right_rc = reverse_complement(unmeth_right)

            probe_u_info = design_probe_from_fixed_primers(
                unmethylated_seq, unmeth_left, unmeth_right_rc, config
            )

        highlighted_u = ''
        if probe_u_info and original_seq:
            highlighted_u = highlight_original_cpgs_in_unmeth_probe(
                original_seq,
                probe_u_info['sequence'],
                probe_u_info['start'] 
            )

        # Validate methylated probe
        if not is_probe_valid(internal_seq):
            continue  # Skip this primer pair

        # Validate unmethylated probe if found
        if probe_u_info and not is_probe_valid(probe_u_info['sequence']):
            probe_u_info = None  # Ignore the probe but keep the primer

        primers.append({
            'left': left_seq,
            'right': right_seq,
            'left_tm': result[f'PRIMER_LEFT_{i}_TM'],
            'right_tm': result[f'PRIMER_RIGHT_{i}_TM'],
            'left_gc': result[f'PRIMER_LEFT_{i}_GC_PERCENT'],
            'right_gc': result[f'PRIMER_RIGHT_{i}_GC_PERCENT'],
            'product_size': result[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
            'cpg_left': count_cpgs(left_seq),
            'cpg_right': count_cpgs(right_seq),
            'cpg_count': count_cpgs(left_seq) + count_cpgs(right_seq),
            'has_cpg_left': 'CG' in left_seq,
            'has_cpg_right': 'CG' in right_seq,
            'left_start': left_pos,
            'degenerate_left': degenerate_primer(left_seq, strand='forward'),
            'degenerate_right': degenerate_primer(right_seq, strand='reverse'),

            'probe_m': internal_seq,
            'probe_m_tm': probe_m_tm,
            'probe_m_gc': probe_m_gc,
            'probe_m_cpgs': internal_cpgs,
            'highlighted_probe_m': highlight_cpg(internal_seq),

            'probe_u': probe_u_info['sequence'] if probe_u_info else '',
            'probe_u_tm': probe_u_info['tm'] if probe_u_info else None,
            'probe_u_gc': probe_u_info['gc'] if probe_u_info else None,
            'probe_u_cpgs': probe_u_info['cpgs'] if probe_u_info else 0,
            'highlighted_probe_u': highlighted_u,

            'left_self_end': result.get(f'PRIMER_LEFT_{i}_SELF_END_TH'),
            'right_self_end': result.get(f'PRIMER_RIGHT_{i}_SELF_END_TH'),
            'probe_m_self_end': result.get(f'PRIMER_INTERNAL_{i}_SELF_END_TH'),
            'probe_u_self_end': probe_u_info['probe_self_end'] if probe_u_info else None,

            'left_hairpin': result.get(f'PRIMER_LEFT_{i}_HAIRPIN_TH'),
            'right_hairpin': result.get(f'PRIMER_RIGHT_{i}_HAIRPIN_TH'),
            'probe_m_hairpin': round(result.get(f'PRIMER_INTERNAL_{i}_HAIRPIN_TH'),2),
            'probe_u_hairpin': probe_u_info['probe_hairpin'] if probe_u_info else None,

            'pair_compl_end': result.get(f'PRIMER_PAIR_{i}_COMPL_END')

        })

    return primers

def design_primers(seq, unmethylated_seq=None, original_seq=None, config=None):

    excluded = find_cpg_excluded_regions(seq)

    # Phase 1: CpG-excluded search
    primers = run_primer3(seq, excluded=excluded, unmethylated_seq=unmethylated_seq, original_seq=original_seq, config=config)
    if primers:
        print("[INFO] Using CpG-excluded primers.")
    else:
        print("[INFO] No CpG-free primers found. Using relaxed fallback.")
        primers = run_primer3(seq, excluded=[], unmethylated_seq=unmethylated_seq, original_seq=original_seq)

    # Sort by primer CpG count, then probe CpGs
    primers.sort(key=lambda p: (p['cpg_count'], -p['probe_m_cpgs']))

    return primers

