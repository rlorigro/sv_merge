"""
A minimal version of $resolve.py$ by Adam English, adapted by Fabio Cunial
"""

import sys
import pysam

REF_CLEAN = True # Set to false if you're working with the right reference


RC = str.maketrans("ATCG", "TAGC")
def do_rc(s):
    """
    Reverse complement a sequence
    """
    return s.translate(RC)[::-1]


def resolve(entry, ref):
    if entry.start > ref.get_reference_length(entry.chrom):
        return None
    if entry.alts[0] == '<INS>':
        # Removing symbolic INS, since they carry no information.
        return None
    
    # Replacing degenerate DNA characters with N, if any.
    seq_prime = ref.fetch(entry.chrom, entry.start, entry.stop)
    seq = ""
    for i in range(len(seq_prime)):
        c = seq_prime[i].upper()
        if c!='A' and c!='C' and c!='G' and c!='T':
            seq = seq + 'N'
        else:
            seq = seq + seq_prime[i]
    
    # Fixing symbolic ALTs
    if entry.alts[0] == '<DEL>':
        entry.ref = seq
        entry.alts = [seq[0]]
    elif entry.alts[0] == '<INV>':
        entry.ref = seq
        entry.alts = [do_rc(seq)]
    elif entry.alts[0] == '<DUP>' or entry.alts[0] == '<CNV>':  
        entry.info['SVTYPE'] = 'INS'
        entry.ref = seq[0]
        entry.alts = [seq]
        entry.stop = entry.start + 1
    
    return entry


if __name__ == '__main__':
    vcf = pysam.VariantFile(sys.argv[1])
    ref = pysam.FastaFile(sys.argv[2])
    n_header = vcf.header.copy()
    if REF_CLEAN:
        for ctg in vcf.header.contigs.keys():
            if ctg not in ref.references:
                n_header.contigs[ctg].remove_header()
    out = pysam.VariantFile("/dev/stdout", 'w', header=n_header)
    for entry in vcf:
        sys.stderr.write(str(entry.pos))
        sys.stderr.write('\n')

        if REF_CLEAN and entry.chrom not in ref.references:
            continue
        
        # Fixing symbolic ALTs
        if entry.alts[0].startswith("<"):
            entry = resolve(entry, ref)
        if entry is None or set(entry.alts[0]) == {'N'}:
            continue
        
        # Uppercasing REF and ALT
        entry.ref = entry.ref.upper()
        entry.alts = [entry.alts[0].upper()]
        
        # Fixing SVLEN
        if 'SVTYPE' in entry.info:
            if entry.info['SVTYPE'] == 'INS' or entry.info['SVTYPE'] == 'DEL':
                entry.info['SVLEN'] = abs(len(entry.ref) - len(entry.alts[0]))
            elif entry.info['SVTYPE'] != 'BND':
                entry.info['SVLEN'] = len(entry.ref)
        
        # Fixing GTs
        n_gt = tuple([_ if _ is not None else 0 for _ in entry.samples[0]['GT']])
        is_phased = entry.samples[0].phased
        entry.samples[0]['GT'] = n_gt
        entry.samples[0].phased = is_phased
        
        # Outputting
        entry.translate(n_header)
        try:
            out.write(entry)
        except Exception as e:
            sys.stderr.write(f"Error: {str(e)}\n")
            sys.stderr.write(f"{entry}\n{type(entry)}\n")
