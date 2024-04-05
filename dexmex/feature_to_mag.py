import os
import csv
import argparse


def check_input(gff, bins_directory, skip, output):
    """
    Checks the input arguments for errors.

    Args:
    gff: Path to the gff file.
    bins_directory: Path to the directory containing the MAG bins.
    skip: Filenames of MAG bins to skip.
    output: Path to the output file.
    """
    if gff is not None and not os.path.isfile(gff):
        raise ValueError('The provided gff file does not exist.')

    if not os.path.isdir(bins_directory):
        raise ValueError('The provided bins_directory does not exist.')

    for bin_file in skip:
        if not os.path.isfile(os.path.join(bins_directory, bin_file)):
            raise ValueError(f'The filename {bin_file} provided using --skip does not exist in the bins_directory {bins_directory}.')
        
    if os.path.isfile(output):
        raise ValueError(f"The output file {output} already exists. Please provide a new output file path.")
        

def get_contig_to_mag(bins_directory: str,
                      bins_to_skip: list) -> dict:
    """
    Returns a dictionary matching contig names (as keys) to MAG names (which
    are derived from filenames in the bins_directory). The bins_to_skip list
    contains filenames of MAG bins to skip.

    Args:
    bins_directory: Path to the directory with FASTAs of MAG bisn
    bins_to_skip: Filenames of MAG bins to skip.

    Returns:
    Dictionary matching contig names (as keys) to MAG names.
    """
    contig_to_mag = {}
    for bin_file in os.listdir(bins_directory):
        if bin_file in bins_to_skip:
            continue

        with open(os.path.join(bins_directory, bin_file), 'r') as f:
            for line in f:
                if line.startswith('>'):
                    contig = line[1:].strip()
                    contig = contig.split(' ')[0]
                    contig_to_mag[contig] = bin_file

    return contig_to_mag


def get_feature_to_mag_from_gff(gff_file: str,
                                feature_fieldname: str,
                                contig_to_mag: dict,
                                gffAttributeDelimiter: str=';',
                                gffAttributeIndex: int=8,
                                contigIndex: int=0) -> dict:
    """
    Returns a dictionary matching feature names (as keys) to MAG names. The
    feature names are derived from the 'Name' attribute in the gff file.

    Args:
    gff_file: Path to the gff file.
    feature_fieldname: The attribute type of the feature names in the gff file.
    contig_to_mag: Dictionary matching contig names (as keys) to MAG names.
    gffAttributeDelimiter: The delimiter used in the attribute field of the gff
        file.
    gffAttributeIndex: The index of the attribute field in the gff file.
    contigIndex: The index of the contig field in the gff file.

    Returns:
    Dictionary matching feature names (as keys) to MAG names.
    """
    feature_to_mag = {}
    with open(gff_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0].startswith('#'):
                continue
            contig = row[contigIndex].strip()
            attribute_strings = row[gffAttributeIndex].split(gffAttributeDelimiter)
            attributes = {}
            for attribute_string in attribute_strings:
                if '=' not in attribute_string:
                    continue
                key, value = attribute_string.split('=')
                attributes[key] = value
            if feature_fieldname not in attributes:
                print(f"ERROR: The target attribute was set to {feature_fieldname}, but it was not found in the gff file.")
                print(f"Please check the gff file. Make sure to set the correct attribute name using the --gene_id argument.")
                exit(1)
            feature = attributes[feature_fieldname]
            feature_to_mag[feature] = contig_to_mag.get(contig, None)    
    return feature_to_mag


def get_feature_to_mag_from_featureCounts(featureCounts_file: str,
                                          contig_to_mag: dict,
                                          geneId_index: int = 0,
                                          contig_index: int = 1) -> dict:
    """
    Returns a dictionary matching feature names (as keys) to MAG names. The
    feature names are derived from the first column of the featureCounts file.

    Args:
    featureCounts_file: Path to the featureCounts file.
    contig_to_mag: Dictionary matching contig names (as keys) to MAG names.
    geneId_index: The index of the gene ID field in the featureCounts file.
    contig_index: The index of the contig field in the featureCounts file.

    Returns:
    Dictionary matching feature names (as keys) to MAG names.
    """
    with open(featureCounts_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # Skip 2 header rows
        next(reader)
        feature_to_mag = {}
        for row in reader:
            feature = row[geneId_index].strip()
            contig = row[contig_index].strip()
            feature_to_mag[feature] = contig_to_mag.get(contig, None)
    return feature_to_mag


def write_output(outfile: str,
                 feature_to_mag: dict,
                 header=['feature_id', 'mag_id']):
    """
    Writes the feature_to_mag dictionary to a file.

    Args:
    outfile: Path to the output file.
    feature_to_mag: Dictionary matching feature names (as keys) to MAG names.
    """
    outdir = os.path.dirname(outfile)
    if not os.path.isdir(outdir):
        if not outdir == '':
            os.makedirs(outdir)
    with open(outfile, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(header)
        for feature, mag in feature_to_mag.items():
            writer.writerow([feature, mag])

    
def build(output: str,
          bins_directory: str,
          bins_to_skip: list,
          featureCounts: str = None,
          gff: str = None,
          gene_id: str = 'gene_id'):
    """
    Creates a feature_to_mag.tsv file as required by local_diffExp.R.

    Args:
    output: Path to the output file.
    bins_directory: Path to the directory containing the MAG bins.
    bins_to_skip: Filenames of MAG bins to skip.
    featureCounts: Path to the featureCounts file.
    gff: Path to the gff file.
    gene_id: The attribute name of the target features in the gff file.
    """
    check_input(gff, bins_directory, bins_to_skip, output)

    if (featureCounts is None) == (gff is None):
        raise ValueError('Exactly one of featureCounts or gff must be provided.')
    
    contig_to_mag = get_contig_to_mag(bins_directory, bins_to_skip)
    if featureCounts is None:
        feature_to_mag = get_feature_to_mag_from_gff(gff, gene_id,
                                                     contig_to_mag)
    else:
        feature_to_mag = get_feature_to_mag_from_featureCounts(featureCounts,
                                                               contig_to_mag)
    write_output(output, feature_to_mag)


def main():
    parser = argparse.ArgumentParser(description='Creates a feature_to_mag.tsv as required by local_diffExp.R, either from featureCounts output or form a gff file.')

    parser.add_argument('-o', '--output', help='Path to the output file.', required=True)
    parser.add_argument('-b', '--bins_directory', help='Path to the directory containing the MAG bin.s', required=True)
    parser.add_argument('-s', '--skip', help='Filenames of MAG bins to skip. When providing multiple filenames, separate them with a space.', nargs='*', default=[])

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-f', '--featureCounts', help='Output of the subread featureCounts tool. If you have multiple featureCounts files, it is sufficient to provide just one.')
    group.add_argument('-g', '--gff', help='A gff file containing gene annotations for all MAG bins.')

    parser.add_argument('-i', '--gene_id', default='gene_id', help='The attribute name of the target features in the gff file.')

    args = parser.parse_args()

    build(args.output,
          args.bins_directory,
          args.skip,
          args.featureCounts,
          args.gff,
          args.gene_id)


if __name__ == '__main__':
    main()
