import csv
import random

debug = True

counts_file = 'counts_sim.tsv'
feature_to_mag_file = 'feature_to_mag_sim.tsv'
coldata_file = 'coldata_sim.tsv'

mags = ['zebra', 'lion', 'sow', 'otter', 'horse']  # the unique identifiers of our MAGs
nr_of_features = 1000  # How many features to generate per MAG

samples    = ['reactor_a_replicate_1',
              'reactor_a_replicate_2',
              'reactor_b_replicate_1',
              'reactor_b_replicate_2',
              'reactor_b_replicate_3',
              'reactor_c_replicate_1',
              'reactor_c_replicate_2',
              ]

sample_to_condition =  {'reactor_a_replicate_1': 'normal',
                        'reactor_a_replicate_2': 'normal',
                        'reactor_b_replicate_1': 'speiseeis',
                        'reactor_b_replicate_2': 'speiseeis',
                        'reactor_b_replicate_3': 'speiseeis',
                        'reactor_c_replicate_1': 'fasting',
                        'reactor_c_replicate_2': 'fasting',
                        }

counts_header = ['feature_id']
for s in samples:
    counts_header.append(s)
counts_rows = [counts_header]

feature_to_mags_header = ['feature_id', 'mag_id']
feature_to_mag_rows = [feature_to_mags_header]


def generate_mag_abundance(mags=mags,
                           sample_to_condition=sample_to_condition):
    """
    Generates random relative abundances for each MAG in each condition.
    Conditions are the unique values in sample_to_condition.
    The sum of abundances of all MAGs in a single condition has to be 1.
    """
    conditions = set(sample_to_condition.values())
    mag_abundance = {mag: {condition: 0.0 for condition in conditions} for mag in mags}
    
    # Generate random abundances such that the sum of abundances in each condition is 1
    for condition in conditions:
        # Generate initial random values
        random_values = [random.uniform(0.0, 1.0) for _ in mags]
        total = sum(random_values)
        
        # Normalize the values so their sum is 1
        normalized_values = [value / total for value in random_values]
        
        # Assign normalized values to each mag for the current condition
        for mag, value in zip(mags, normalized_values):
            mag_abundance[mag][condition] = value

    return mag_abundance


mag_abundance = generate_mag_abundance()


def generate_counts(mag,
                    feature_name="unknown",
                    samples=samples,
                    sample_to_condition=sample_to_condition,
                    mag_abundance=mag_abundance):
    
    conditions = set(sample_to_condition.values())
    condition_status = {}
    for condition in conditions:
        condition_status[condition] = random.choices(["normal", "increased", "decreased"],
                                                     weights=(.9,.05,.05))
    
    baseline_value = random.uniform(0.0,1000.0)
  
    sample_values = {}
    for sample, condition in sample_to_condition.items():
        status = condition_status[condition][0]

        if status == "normal":
            sample_value = baseline_value * random.uniform(0.5,  1.5)
        elif status == "increased":
            sample_value = baseline_value * random.uniform(2.0, 10.0)
        elif status == "decreased":
            sample_value = baseline_value * random.uniform(0.1,  0.5)
        else:
            print("ILLEGAL STATUS", status)
            exit(1)

        mag_abudannce_factor = mag_abundance[mag][condition]
        sample_value *= mag_abudannce_factor
        
        sample_value = round(sample_value, 2)
        sample_values[sample] = sample_value
        #print(f"Status {status}\tval {sample_value}")

    if debug:
        print(f"\nFEATURE::     {feature_name}")
        print(f"BASELINE::    {baseline_value}", end="")
        # MAG data
        print("\nMAG::".ljust(15), end="")
        for k, v in mag_abundance[mag].items():
            print(f"{k}: {str(round(v,2))}".ljust(20), end="")
        # Expression levels
        print("\nEXPRESSION::".ljust(15), end="")
        for k, v in condition_status.items():
            print(f"{k}: {str(v[0])}".ljust(20), end="")
        # Results
        print("\nRESULTS::".ljust(15), end="")
        for smp, val in sample_values.items():
            c = sample_to_condition[smp]
            print(f"{c}: {str(int(val))} | ", end="")
        print("")

    return sample_values, baseline_value



for mag in mags:
    for fnumber in range(nr_of_features+1):
        fname = f"{mag}_{fnumber}"
        mapping_row = [fname, mag]
        feature_to_mag_rows.append(mapping_row)

        row = [fname]
        sample_to_count, baseline = generate_counts(mag, feature_name=fname) 
        for sample in samples:
            row.append(sample_to_count[sample])
        counts_rows.append(row)

with open(counts_file, 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    for row in counts_rows:
        writer.writerow(row)

with open(feature_to_mag_file, 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    for row in feature_to_mag_rows:
        writer.writerow(row)

with open(coldata_file, 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['sample', 'condition'])
    for sample, condition in sample_to_condition.items():
        writer.writerow([sample, condition])


print("")
print("\n################# MAG ABUNDANCES ############################")
for mag, abddata in mag_abundance.items():
    print(mag)
    for condition, value in abddata.items():
        print(f"{condition.rjust(15)}: {str(round(value,2)).rjust(5)}")

