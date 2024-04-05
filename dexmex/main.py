#!/usr/bin/env python

import os
import argparse
import subprocess

from dexmex import feature_to_mag
from dexmex import convert_featureCounts


def main():
    parser = argparse.ArgumentParser(description='Calculating local (per-bin-normalized) differential expression using DESeq2.')
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # local differential expression command
    local_diff_exp_parser = subparsers.add_parser('localdiffexp', help='Run the local_diffExp.R script')
    local_diff_exp_parser.add_argument('--outdir', required=True, help='Path to the output directory.')
    local_diff_exp_parser.add_argument('--coldata_path', required=True, help='Path to the coldata file.')
    local_diff_exp_parser.add_argument('--counts_path', required=True, help='Path to the counts file.')
    local_diff_exp_parser.add_argument('--feature_to_mag_path', required=True, help='Feature to MAG data path.')
    local_diff_exp_parser.add_argument('--reference_level', required=True, help='Reference level for condition comparison.')
        
    # convertfc command
    convert_fc_parser = subparsers.add_parser('convertfc', help='Access functionality of convert_featureCounts.py')
    convert_fc_parser.add_argument('-f', '--featureCountsFile', nargs='+', help='Paths to featureCounts files. You can provide multiple files, separated by spaces.', required=True)
    convert_fc_parser.add_argument('-s', '--samplenames', nargs='*', help='You can overwrite the samplenames in the featureCounts tables with these names. Provide the names in the same order as the files, separated by spaces. When providing multiple names for one file (because it has multiple count columns), separate them with a comma (no space).' )
    convert_fc_parser.add_argument('-c', '--coldata', help='Path to the a coldata file as used by deseq. If provided, output columns will be in the same order as the samples in the coldata file.', required=False)
    convert_fc_parser.add_argument('-o', '--output', help='Path to the output file.', required=True)
    
    # featuretomag command
    feature_to_mag_parser = subparsers.add_parser('featuretomag', help='Access functionality of feature_to_mag.py')
    feature_to_mag_parser.add_argument('-o', '--output', help='Path to the output file.', required=True)
    feature_to_mag_parser.add_argument('-b', '--bins_directory', help='Path to the directory containing the MAG bins.', required=True)
    feature_to_mag_parser.add_argument('-s', '--skip', help='Filenames of MAG bins to skip. Separate multiple filenames with a space.', nargs='*', default=[])
    group = feature_to_mag_parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-f', '--featureCounts', help='Output of the subread featureCounts tool. If you have multiple featureCounts files, just provide one of them.')
    group.add_argument('-g', '--gff', help='A gff file containing gene annotations for all MAG bins.')
    feature_to_mag_parser.add_argument('-i', '--gene_id', default='gene_id', help='The attribute name of the target features in the gff file.')
    
    args = parser.parse_args()
    
    if args.command == 'convertfc':
        print("Tryin to convert featureCounts files to a single count table...")
        convert_featureCounts.convert(args.featureCountsFile,
                                      args.output,
                                      args.samplenames,
                                      args.coldata)
        print("\n\t...done!")
    elif args.command == 'featuretomag':
        print("Starting feature to MAG mapping..")
        feature_to_mag.build(args.output,
                             args.bins_directory,
                             args.skip,
                             args.featureCounts,
                             args.gff,
                             args.gene_id)
        print("\n\t...done!")
    elif args.command == 'localdiffexp':
        print("Starting local differential expression analysis...")
        mainpy_dir = os.path.dirname(os.path.abspath(__file__))
        rscript_path = os.path.join(mainpy_dir, 'local_diffExp.R')
        command = [
            'Rscript', rscript_path,
            '--outdir', args.outdir,
            '--coldata_path', args.coldata_path,
            '--counts_path', args.counts_path,
            '--feature_to_mag_path', args.feature_to_mag_path,
            '--reference_level', args.reference_level
        ]
        subprocess.run(command)
        print("\n\t...done!")
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
