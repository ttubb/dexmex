#!/usr/bin/env python

import argparse

from dexmex import feature_to_mag
from dexmex import convert_featureCounts


def main():
    parser = argparse.ArgumentParser(description='dexmex: A tool to access various functionalities for processing and analyzing metagenomic data.')
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # convertfc command
    convert_fc_parser = subparsers.add_parser('convertfc', help='Access functionality of convert_featureCounts.py')
    convert_fc_parser.add_argument('-f', '--featureCountsFile', nargs='+', help='Paths to featureCounts files. You can provide multiple files, separated by spaces.', required=True)
    convert_fc_parser.add_argument('-s', '--samplenames', nargs='*', help='You can overwrite the samplenames in the featureCounts tables with these names. Provide the names in the same order as the files, separated by spaces. When providing multiple names for one file (because it has multiple count columns), separate them with a comma (no space).' )
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
        convert_featureCounts.convert(args.featureCountsFile,
                                      args.output,
                                      args.samplenames)
    elif args.command == 'featuretomag':
        feature_to_mag.build(args.output,
                             args.bins_directory,
                             args.skip,
                             args.featureCounts,
                             args.gff,
                             args.gene_id)
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
