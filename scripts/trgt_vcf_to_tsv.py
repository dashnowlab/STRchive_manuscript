"""Convert a TRGT output VCF to two TSVs. One for allele lengths and one for motif counts.
"""

import argparse
import gzip
import re
from tqdm import tqdm

def main():
    p = argparse.ArgumentParser(
        description="Convert a TRGT output VCF to two TSVs. One for allele lengths and one for motif counts."
    )
    p.add_argument("--verbose", action="store_true", help="Print verbose output")
    p.add_argument("--sample-id",
                   help="If not specified, the sample id will be parsed from the last column of the vcf header.")
    p.add_argument("--tsv-prefix", help="Output tsv path prefix (.tsv will be added)", default="all_samples")
    p.add_argument("vcf_paths", help="TRGT vcf path", nargs="+")
    
    args = p.parse_args()

    all_locus_results = []

    for vcf_path in args.vcf_paths:

        print(f"Processing {vcf_path}")
        locus_results = process_trgt_vcf(
            vcf_path,
            sample_id=args.sample_id,
            verbose=args.verbose,
        )
        all_locus_results.append(locus_results)

        AL_tsv_path = f"{args.tsv_prefix}.allele_lengths.tsv"
        MC_tsv_path = f"{args.tsv_prefix}.motif_counts.tsv"

        with open(AL_tsv_path, "w") as AL_tsv, open(MC_tsv_path, "w") as MC_tsv:
            ALheader = [
                "Sample", "Sex", 
                "Chrom", "Pos", "Locus",  
                "Structure", "AlleleIndex",  
                "AlleleLengthBP", "AlleleLengthRange"
                      ]
            MCheader = [
                "Sample", "Sex", 
                "Chrom", "Pos", "Locus", 
                "LocusMotifId", "AlleleIndex",  
                "Motif", "MotifCount"
                      ]
            AL_tsv.write("\t".join(ALheader) + "\n")
            MC_tsv.write("\t".join(MCheader) + "\n")
            for sampledict in all_locus_results:
                sample = sampledict['Sample']
                for locus in sampledict['LocusResults']:
                    locusdict = sampledict['LocusResults'][locus]
                    for variant in locusdict['Variants']:
                        variantdict = locusdict['Variants'][variant]
                        for i in range(len(variantdict['AL'])):
                            ALoutlist = [
                                sample, sampledict['Sex'], 
                                locusdict['Chrom'], locusdict['Pos'], locusdict['LocusId'], 
                                variantdict['STRUC'], i,
                                variantdict['AL'][i], variantdict['ALLR'][i]
                                       ]
                            AL_tsv.write("\t".join(map(str, ALoutlist)) + "\n")
                            for motif in locusdict['MotifIds']:
                                motifdict = locusdict['MotifIds'][motif]
                                MCoutlist = [
                                    sample, sampledict['Sex'], 
                                    locusdict['Chrom'], locusdict['Pos'], locusdict['LocusId'],
                                    motifdict['MotifId'], i,
                                    motifdict['Motif'], motifdict['MotifCount'][i]
                                        ]
                                MC_tsv.write("\t".join(map(str, MCoutlist)) + "\n")

# Convert these to tests:
# example_format, example_genotype = "GT:AL:ALLR:SD:MC:MS:AP:AM       1/2:29,38:29-31,38-41:13,16:10,13:0(0-29),0(0-38):0.966667,0.974359:0.01,.".split()
# example_format, example_genotype = "GT:AL:ALLR:SD:MC:MS:AP:AM       1/1:49,49:47-50,48-50:27,31:0_10,0_10:1(0-49),1(0-49):0.980000,0.980000:.,.".split()
# example_format, example_genotype = "GT:AL:ALLR:SD:MC:MS:AP:AM   0:48:47-51:32:16:0(0-48):0.666667:0.22".split()
# print(parse_TRGT_genotype(example_format, example_genotype))
def parse_TRGT_genotype(format, genotype):
    """Parse the TRGT genotype fields into a dictionary.
    """
    genotype_fields = format.split(":")
    genotype_values = genotype.split(":")
    genotype_dict = dict(zip(genotype_fields, genotype_values))
    
    for key, value in genotype_dict.items():
        if "," in value:
            genotype_dict[key] = value.split(",")
        else:
            genotype_dict[key] = [value]
    if "MC" in genotype_dict:
        for mc_i, mc in enumerate(genotype_dict["MC"]):
            genotype_dict["MC"][mc_i] = [c for c in mc.split("_")]
    return genotype_dict

def process_trgt_vcf(vcf_path, sample_id=None, verbose=False):
    """Convert a TRGT output VCF to a list of dictionaries, one for each allele.

    Output format:
    {
        "LocusResults": {
            "locus_id": {
                "AlleleCount": int,
                "LocusId": str,
                "Variants": {
                    "locus_id": {
                        "STRUC": str,
                        "GT": str,
                        "AL": list,
                        "ALLR": list,
                        "SD": list,
                        "MC": list,
                        "MS": list,
                        "AP": list,
                        "AM": list,
    """
    locus_results = {
        "LocusResults": {}, # Store all the locus results for this sample in a dictionary
        "Sample": sample_id,
        "Sex": None,
    }

    fopen = gzip.open if vcf_path.endswith("gz") else open

    with fopen(vcf_path, "rt") as vcf:
        for line in vcf:
            if line.startswith("##trgtCommand"):
                if 'karyotype' in line:
                    # find the value after the karyotype flag
                    locus_results['Sex'] = re.findall(r'karyotype (XY|XX)', line)[0]
            if line.startswith("#CHROM"):
                header_fields = line.strip().split("\t")
                if sample_id is None and len(header_fields) == 10:
                    print(f"Got sample id '{header_fields[9]}' from the VCF header")
                    locus_results["Sample"] = header_fields[9]
            if not line.startswith("#"):
                break

    with fopen(vcf_path, "rt") as vcf:
        line_counter = 0
        if verbose:
            vcf = tqdm(vcf, unit=" vcf records", unit_scale=True, unit_divisor=1000)

        for line in vcf:
            if line.startswith("#"):
                continue

            line_counter += 1
            fields = line.strip().split("\t")
            chrom = fields[0]
            start_1based = int(fields[1])
            info = fields[7]
            if not fields[9] or fields[9] == ".":  # no genotype
                continue

            info_dict = dict([key_value.split("=") for key_value in info.split(";")])
            genotype_dict = parse_TRGT_genotype(fields[8], fields[9])

            # GT:AL:ALLR:SD:MC:MS:AP:AM
            if genotype_dict["AL"] == ".":   # no genotype
                continue

            try:
                locus_id = info_dict["TRID"]
                end_1based = int(info_dict["END"])
                    
                locus_results["LocusResults"][locus_id] = {
                    "AlleleCount": len(genotype_dict["AL"]),
                    "LocusId": locus_id,
                    "Variants": {
                        locus_id: info_dict | genotype_dict
                    },
                    "MotifIds": {},
                    "Chrom": chrom,
                    "Pos": start_1based,
                }

                motifs = info_dict["MOTIFS"].split(",")
                for motif_i, motif in enumerate(motifs):
                    motif_id = f"{locus_id}-m{motif_i}-{motif}"
                    locus_results["LocusResults"][locus_id]["MotifIds"][motif_id] = {
                            "MotifId": motif_id,
                            "MotifCount": [x[motif_i] for x in genotype_dict["MC"]],
                            #"MotifSpans": [x[motif_i] for x in genotype_dict["MS"]],
                            "Motif": motif,
                    }

            except Exception as e:
                print(f"Error on vcf record #{line_counter}: {e}")
                print(line)
                print(genotype_dict)
                import traceback
                traceback.print_exc()

    return locus_results


if __name__ == "__main__":
    main()
