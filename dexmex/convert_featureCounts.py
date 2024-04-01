import os
import csv
import argparse


def check_input(featureCounts_files: list,
                sample_names: list,
                output_file: str):
    """
    Checks the input arguments for errors.

    Args:
    featureCounts_files: List of paths to the featureCounts files.
    output_file: Path to the output file.
    """
    if sample_names:
        if not len(featureCounts_files) == len(sample_names):
            print("files",featureCounts_files)
            print("snames",sample_names)
            raise ValueError("The number of featureCounts files and sample name lists must be equal.")

    for file in featureCounts_files:
        if not os.path.isfile(file):
            raise ValueError(f"The provided featureCounts file {file} does not exist.")

    if os.path.isfile(output_file):
        raise ValueError(f"The output file {output_file} already exists. Please provide a new output file path.")


def parse_feature_counts(feature_counts_file: str,
                         sample_names: list,
                         first_sample_index: int = 6,
                         feature_name_index: int = 0) -> dict:
    """
    Parses a featureCounts file and builds a nested dictionary. The outer dict
    has samples as keys, the inner dicts contains feature names and their
    respective counts.

    Args:
    feature_counts_file: Path to the featureCounts file.
    first_sample_index: The index of the first sample name column in the file.
    feature_name_index: The index of the feature name column in the file.

    Returns:
    Nested dictionary with samples as keys, containing feature names and
        their respective counts.
    """
    feature_counts = {}
    if sample_names is None:
        newnames = None
    else:
        newnames = sample_names.split(',')
    with open(feature_counts_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        if header[0].startswith('#'):
            header = next(reader)
        sample_names = header[first_sample_index:]
        if newnames:
            if not (len(newnames) == len(sample_names)):
                raise ValueError(f"The number of provided sample names does not match the number of samples in the file {feature_counts_file}\nIn file ({len(sample_names)}): {sample_names}\nNames you provided ({len(newnames)}): {newnames}")
            sample_names = newnames
        for sname in sample_names:
            feature_counts[sname] = {}
        for row in reader:
            feature_name = row[feature_name_index].strip()
            for i, sname in enumerate(sample_names):
                count = float(row[i + first_sample_index])
                feature_counts[sname][feature_name] = count

    return feature_counts


def merge_feature_counts(feature_counts_list: list) -> dict:
    """
    Validates and merges multiple feature counts dictionaries by ensuring:
    - No duplicate samples across files.
    - The set of feature names is identical across all samples and all files.

    Args:
    feature_counts_list: List of feature counts dictionaries.

    Returns:
    A single merged dictionary with samples as keys, containing feature names
    and their respective counts, if the above conditions are met.

    Raises:
    ValueError: If the conditions are not met.
    """
    all_features_set = set()    # stores all unique feature names
    merged_counts = {}          # stores all feature counts

    # First pass to collect all unique feature names and check for duplicate samples
    for feature_counts in feature_counts_list:
        for sample, features in feature_counts.items():
            if sample in merged_counts:
                raise ValueError(f"Duplicate condition found: {sample}. Each condition must be unique across files.")
            merged_counts[sample] = features
            all_features_set.update(features.keys())

    # Second pass to ensure each condition has the identical set of feature names
    for sample, features in merged_counts.items():
        if set(features.keys()) != all_features_set:
            raise ValueError("Feature set mismatch. All conditions must have the same set of features.")

    return merged_counts


def write_output(merged_counts: dict,
                 output_file: str):
    """
    Writes the merged feature counts to an output file as tsv.

    Args:
    merged_counts: The merged dictionary containing the feature counts.
    output_file: Path to the output file.
    """
    output_dir = os.path.dirname(output_file)
    if not os.path.exists(output_dir):
        if not output_dir == '':
            os.makedirs(output_dir)
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        samples = list(merged_counts.keys())
        features = list(merged_counts[samples[0]].keys())
        header = ['feature_id'] + samples
        writer.writerow(header)
        for feature in features:
            row = [feature] + [merged_counts[sample].get(feature, 0) for sample in samples]
            writer.writerow(row)


def convert(fc_files: list,
            outfile: str,
            sample_names: list = None):
    """
    Converts multiple featureCounts files to a single count table.

    Args:
    fc_files: List of paths to the featureCounts files.
    outfile: Path to the output file.
    """
    check_input(fc_files, sample_names, outfile)

    feature_counts_list = []

    if sample_names:
        for file, newnames in zip(fc_files, sample_names):
            file_results = parse_feature_counts(file, newnames)
            feature_counts_list.append(file_results)
    else:
        for file in fc_files:
            file_results = parse_feature_counts(file, None)
            feature_counts_list.append(file_results)

    merged_counts = merge_feature_counts(feature_counts_list)
    write_output(merged_counts, outfile)


def main():
    parser = argparse.ArgumentParser(description='Creates a count table compatible to local_diffExp.R from multiple featureCounts tables.')
    parser.add_argument('-f', '--featureCountsFile', nargs='+', help='Paths to the featureCounts files.', required=True)
    parser.add_argument('-s', '--samplenames', nargs='*', help='You can overwrite the samplenames in the featureCounts tables with these names. Provide the names in the same order as the files, separated by spaces. When providing multiple names for one file (because it has multiple count columns), separate them with a comma (no space).' )
    parser.add_argument('-o', '--output', help='Path to the output file.', required=True)

    args = parser.parse_args()

    convert(args.featureCountsFile, args.output)


if __name__ == "__main__":
    main()
