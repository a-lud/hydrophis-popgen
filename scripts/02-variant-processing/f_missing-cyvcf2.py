#!/usr/bin/env python3

from sys import stdin
from cyvcf2 import VCF, Writer
import numpy as np
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="""Filter VCF file for variants with at most 1 missing
        genotype per species group. This script was written as BCFtools +fill-tags
        is not working correctly at the time of writing.""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        'vcf',
        nargs='?',
        default='-',
        help='VCF file path (use "-" for stdin, default: stdin)'
    )
    
    parser.add_argument(
        '--species1', '-s1',
        required=True,
        help='First species accounting for samples 1-5 (HCU/HEL are valid)'
    )
    
    parser.add_argument(
        '--species2', '-s2', 
        required=True,
        help='Second species accounting for samples 6-10 (HEL when S1 is HCU, or HMA)'
    )

    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Filtered VCF file name'
    )
    
    args = parser.parse_args()
    
    # Handle stdin case
    if args.vcf == '-':
        args.vcf = stdin
    
    return args

if __name__ == "__main__":
    args = parse_arguments()
    
    samples = {
        "HCU": ["HCU-1VT", "HCU-2VT", "HCU-3VT", "HCU-4VT", "HCU-7VT"],
        "HEL": ["HEL-KLS1024", "HEL-KLS1025", "HEL-KLS1027", "HEL-KLS1034", "HEL-KLS1121"],
        "HMA": ["HMA-KLS0629", "HMA-KLS0630", "HMA-KLS0631", "HMA-KLS0647", "HMA-KLS1020"]
    }

    S1 = samples[args.species1]
    S2 = samples[args.species2]

    vcf = VCF(args.vcf)
    
    # Adding INFO header
    vcf.add_info_to_header({
        "ID": f"F_MISSING_{args.species1}", 
        "Description": f"Added by cyvcf2 expression F_MISSING=F_MISSING in {args.species1}",
        "Type": "Float",
        "Number": "1"
    })
    vcf.add_info_to_header({
        "ID": f"F_MISSING_{args.species2}", 
        "Description": f"Added by cyvcf2 expression F_MISSING=F_MISSING in {args.species2}",
        "Type": "Float",
        "Number": "1"
    })
    vcf.add_info_to_header({
        "ID": "F_MISSING", 
        "Description": "Added by cyvcf2 expression F_MISSING=F_MISSING",
        "Type": "Float",
        "Number": "1"
    })

    w = Writer(args.output, vcf)

    for variant in vcf:
        s1_idx = [vcf.samples.index(s) for s in S1]
        s2_idx = [vcf.samples.index(s) for s in S2]

        gt_s1 = variant.gt_types[s1_idx]
        gt_s2 = variant.gt_types[s2_idx]
        
        nmiss_s1 = np.sum(gt_s1 == 2)
        fmiss_s1 = nmiss_s1 / len(S1)

        nmiss_s2 = np.sum(gt_s2 == 2)
        fmiss_s2 = nmiss_s2 / len(S2)

        # genotype 2 is 'UNKNOWN' (i.e missing)
        if nmiss_s1 <= 1 and nmiss_s2 <= 1:
            variant.INFO[f"F_MISSING_{args.species1}"] = fmiss_s1
            variant.INFO[f"F_MISSING_{args.species2}"] = fmiss_s2
            variant.INFO["F_MISSING"] = (nmiss_s1 + nmiss_s2) / (len(S1) + len(S2))
            w.write_record(variant)

    w.close()