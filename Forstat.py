import pandas as pd
import json
from collections import defaultdict

# Read from a JSON file
with open('kits_metadata.json', 'r') as f:
    kits_metadata = json.load(f)

# Convert lists back to sets
for kit, kit_info in kits_metadata.items():
    for locus_type, loci in kit_info.items():
        kits_metadata[kit][locus_type] = set(loci)


class FORSTAT:
    def __init__(self, kit, population_group, sample_number, autosomalSTR=None, sexLinkedXSTR=None, sexLinkedYSTR=None, sexLinkedXInDels=None, sexLinkedYInDels=None, mitochondrialSNPs=None, mitochondrialInDels=None, mitochondrialSequence=None):
        self.kit = kit
        self.population_group = population_group
        self.sample_number = sample_number
        self.autosomalSTR = autosomalSTR if autosomalSTR is not None else {}
        self.sexLinkedXSTR = sexLinkedXSTR if sexLinkedXSTR is not None else {}
        self.sexLinkedYSTR = sexLinkedYSTR if sexLinkedYSTR is not None else {}
        self.sexLinkedXInDels = sexLinkedXInDels if sexLinkedXInDels is not None else {}
        self.sexLinkedYInDels = sexLinkedYInDels if sexLinkedYInDels is not None else {}
        self.mitochondrialSNPs = mitochondrialSNPs if mitochondrialSNPs is not None else {}
        self.mitochondrialInDels = mitochondrialInDels if mitochondrialInDels is not None else {}
        self.mitochondrialSequence = mitochondrialSequence
        self.analyseObject = True

    def check_and_correct_alleles(self, consider_multialleles):
        for attr in ['autosomalSTR', 'sexLinkedXSTR', 'sexLinkedYSTR']:
            str_data = getattr(self, attr)
            for locus, alleles in str_data.items():
                if len(alleles) == 1:
                    str_data[locus] = alleles * 2
                elif len(alleles) > 2 and not consider_multialleles:
                    self.analyseObject = False

#Determine if the data is an STR or mitochondrial input (e.g. convert to float or sting)
def process_alleles(allele_data, is_mitochondrial=False):
        # If allele_data is a single number (integer or float), convert it to a list of floats
        if isinstance(allele_data, (int, float)):
            return [float(allele_data), float(allele_data)]
        # If allele_data is a string of comma-separated values, split it into a list of floats
        elif isinstance(allele_data, str):
            # If this is mitochondrial data, return the string as is
            if is_mitochondrial:
                return allele_data
            # Otherwise, convert to floats and ignore empty strings
            else:
                alleles = [float(allele) for allele in allele_data.split(',') if allele]
                # If there is only one allele, duplicate it
                if len(alleles) == 1:
                    alleles *= 2
                return alleles

def read_data_excel(file_path, kit_name, consider_multialleles):
    data = pd.read_excel(file_path)
    data.rename(columns={'S No.': 'SampleNumber', 'Population': 'PopulationGroup'}, inplace=True)

    result = []  # initialize result as an empty list

    for index, row in data.iterrows():
        population_group = row['PopulationGroup']
        sample_number = row['SampleNumber']

        autosomalSTR = {}
        sexLinkedXSTR = {}
        sexLinkedYSTR = {}
        sexLinkedXInDels = {}
        sexLinkedYInDels = {}
        mitochondrialSNPs = {}
        mitochondrialInDels = {}
        mitochondrialSequence = {}

        for column_name, allele_data in row.drop(['SampleNumber', 'PopulationGroup']).items():
            if not pd.isna(allele_data):
                is_mitochondrial = column_name in kits_metadata[kit_name]['mitochondrialSNPs'] or column_name in kits_metadata[kit_name]['mitochondrialInDels']
                alleles = process_alleles(allele_data, is_mitochondrial)

                if column_name in kits_metadata[kit_name]['autosomalSTR']:
                    autosomalSTR[column_name] = alleles
                elif column_name in kits_metadata[kit_name]['sexLinkedXSTR']:
                    sexLinkedXSTR[column_name] = alleles
                elif column_name in kits_metadata[kit_name]['sexLinkedYSTR']:
                    sexLinkedYSTR[column_name] = alleles
                elif column_name in kits_metadata[kit_name]['sexLinkedXInDels']:
                    sexLinkedXInDels[column_name] = alleles
                elif column_name in kits_metadata[kit_name]['sexLinkedYInDels']:
                    sexLinkedYInDels[column_name] = alleles
                elif column_name in kits_metadata[kit_name]['mitochondrialSNPs']:
                    mitochondrialSNPs[column_name] = alleles
                elif column_name in kits_metadata[kit_name]['mitochondrialInDels']:
                    mitochondrialInDels[column_name] = alleles
                elif column_name in kits_metadata[kit_name]['mitochondrialSequence']:
                    mitochondrialSequence[column_name] = alleles

        forensic_kit = FORSTAT(kit_name, population_group, sample_number, autosomalSTR, sexLinkedXSTR, sexLinkedYSTR, sexLinkedXInDels, sexLinkedYInDels, mitochondrialSNPs, mitochondrialInDels, mitochondrialSequence)
        forensic_kit.check_and_correct_alleles(consider_multialleles)

        result.append(forensic_kit)  # append the FORSTAT object to the list
    return result

def read_genepop_shorttandemrepeats(file_path, kit_name, digit_format, consider_multialleles):
    # GENEPOP data is in a text file and we use a reader to access the data
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Dynamically extract loci from the file
    loci = []
    # First line is the file title so we skip that line
    for line in lines[1:]:
        # We want each line to be a unique locus until we find pop. Using lowercase to evaluate these until it is found
        if 'pop' not in line.lower():
            loci.append(line.strip())
        else:
            break

    # The remaining lines contain the data
    data_lines = lines[len(loci) + 2:]
    forensic_kits = []

    # Initialize a counter for populations
    pop_counter = 1
    # Initialize a counter for samples
    sample_counter = 1
    # Initialise the first pop group
    population_group = f'pop_{pop_counter}'

    # Moving through each line we want to know if there is a new population or if we should record the line as a sample
    for line in data_lines:
        # ID a new population group
        if 'pop' in line.lower():
            pop_counter += 1
            sample_counter = 1  # Reset the sample counter for the new population group
            population_group = f'pop_{pop_counter}'
        # If we find a comma in the line it is likely a sample that separates sample name and alleles
        elif ',' in line:
            line_data = line.split(',')
            # Use the sample counter as the key instead of the original sample name. This is to avoid duplicate sample names
            sample_number = f'sample_{sample_counter}'
            # Increment the sample counter
            sample_counter += 1
            # Split the string at white spaces and remove other characters (whitespaces)
            allele_data = line_data[1].strip().split()

            # Initialize the dictionaries
            autosomalSTR = {}
            sexLinkedXSTR = {}
            sexLinkedYSTR = {}
            sexLinkedXInDels = {}
            sexLinkedYInDels = {}
            mitochondrialSNPs = {}
            mitochondrialInDels = {}

            # Determine what type of marker and place it in the correct attribute
            for i, locus in enumerate(loci):
                alleles = [int(allele_data[i][j:j + digit_format]) for j in range(0, len(allele_data[i]), digit_format)]
                if locus in kits_metadata[kit_name]['autosomalSTR']:
                    autosomalSTR[locus] = alleles
                elif locus in kits_metadata[kit_name]['sexLinkedXSTR']:
                    sexLinkedXSTR[locus] = alleles
                elif locus in kits_metadata[kit_name]['sexLinkedYSTR']:
                    sexLinkedYSTR[locus] = alleles
                elif locus in kits_metadata[kit_name]['sexLinkedXInDels']:
                    sexLinkedXInDels[locus] = alleles
                elif locus in kits_metadata[kit_name]['sexLinkedYInDels']:
                    sexLinkedYInDels[locus] = alleles
                elif locus in kits_metadata[kit_name]['mitochondrialSNPs']:
                    mitochondrialSNPs[locus] = alleles
                elif locus in kits_metadata[kit_name]['mitochondrialInDels']:
                    mitochondrialInDels[locus] = alleles

            # Make a FORSTAT object with the data above
            forensic_kit = FORSTAT(kit_name, population_group, sample_number, autosomalSTR, sexLinkedXSTR,
                                   sexLinkedYSTR, sexLinkedXInDels, sexLinkedYInDels, mitochondrialSNPs,
                                   mitochondrialInDels)

            # Check and correct alleles according to the 'consider_multialleles' flag
            forensic_kit.check_and_correct_alleles(consider_multialleles)

            # Add the object to the list
            forensic_kits.append(forensic_kit)

    # Return the list of FORSTAT objects after reading the whole file
    return forensic_kits

@staticmethod
def calculate_allele_frequencies(data_list, locus_type_input):
    allele_counts = {}
    sample_counts = {}

    # Define what types of loci we can input
    locus_types = ['autosomalSTR', 'sexLinkedXSTR', 'sexLinkedYSTR', 'sexLinkedXInDels', 'sexLinkedYInDels', 'mitochondrialSNPs', 'mitochondrialInDels']
    if locus_type_input != 'All':
        # If the input is not all, then set the input to the specific loci the user specified
        locus_types = [locus_type_input]

    # Loop through the list of FORSTAT objects
    for sample_object in data_list:
        if sample_object.analyseObject == True:
            population_group = sample_object.population_group  # Get the population group of the current FORSTAT object
            if population_group not in allele_counts:
                # If this population group is not yet in our dictionaries, add it with empty / 0 values
                allele_counts[population_group] = {}
                sample_counts[population_group] = 0

            sample_counts[population_group] += 1

            for locus_type in locus_types:
                loci = getattr(sample_object, locus_type)  # Access the attribute of the current locus type
                for locus, alleles in loci.items():
                    if locus not in allele_counts[population_group]:
                        allele_counts[population_group][locus] = {}

                    for allele in alleles:
                        if allele is None:
                            continue

                        if allele not in allele_counts[population_group][locus]:
                            allele_counts[population_group][locus][allele] = 0

                        allele_counts[population_group][locus][allele] += 1


    # Calculation of frequencies
    allele_frequencies = {}
    for population_group, loci in allele_counts.items():
        allele_frequencies[population_group] = {}
        for locus, alleles in loci.items():
            allele_frequencies[population_group][locus] = {}
            total_allele_count = sum(alleles.values())
            for allele, count in alleles.items():
                frequency = count / (total_allele_count)
                allele_frequencies[population_group][locus][allele] = frequency
    return allele_frequencies

@staticmethod
def calculate_homozygote_heterozygote_frequencies(data_list, locus_type_input):
    homozygote_counts = {}
    allele_counts = {}
    sample_counts = {}

    locus_types = ['autosomalSTR', 'sexLinkedXSTR', 'sexLinkedYSTR', 'sexLinkedXInDels', 'sexLinkedYInDels',
                   'mitochondrialSNPs', 'mitochondrialInDels']
    if locus_type_input != 'All':
        locus_types = [locus_type_input]

    # Loop through the list of FORSTAT objects
    for sample_object in data_list:
        if sample_object.analyseObject:
            population_group = sample_object.population_group  # Get the population group of the current FORSTAT object
            if population_group not in homozygote_counts:
                homozygote_counts[population_group] = {}
                allele_counts[population_group] = {}
                sample_counts[population_group] = 0

            sample_counts[population_group] += 1

            for locus_type in locus_types:
                loci = getattr(sample_object, locus_type)
                for locus, alleles in loci.items():
                    if locus not in homozygote_counts[population_group]:
                        homozygote_counts[population_group][locus] = 0
                        allele_counts[population_group][locus] = {}

                    # If all alleles are the same, count as homozygote
                    if len(set(alleles)) == 1:
                        homozygote_counts[population_group][locus] += 1

                    # Count alleles
                    for allele in alleles:
                        if allele not in allele_counts[population_group][locus]:
                            allele_counts[population_group][locus][allele] = 0
                        allele_counts[population_group][locus][allele] += 1

    frequencies = {}
    for population_group, loci in homozygote_counts.items():
        frequencies[population_group] = {}
        for locus, homozygote_count in loci.items():
            homozygote_frequency = homozygote_count / sample_counts[population_group]
            heterozygote_frequency = 1 - homozygote_frequency

            # Calculate expected heterozygosity
            total_alleles = sum(allele_counts[population_group][locus].values())
            allele_frequencies = {allele: count / total_alleles for allele, count in
                                  allele_counts[population_group][locus].items()}
            expected_heterozygosity = 1.0 - sum(freq ** 2 for freq in allele_frequencies.values())

            # Apply correction factor if sample size < 30
            if sample_counts[population_group] < 30:
                expected_heterozygosity = expected_heterozygosity * 2 * sample_counts[population_group] / (
                            2 * sample_counts[population_group] - 1)

            frequencies[population_group][locus] = {
                'homozygote_frequency': homozygote_frequency,
                'heterozygote_frequency': heterozygote_frequency,
                'expected_heterozygosity': expected_heterozygosity
            }
    return frequencies

@staticmethod
def calculate_match_probability_and_power_of_discrimination(data_list, locus_type_input):
    locus_genotype_counts = {}
    sample_counts = {}

    locus_types = ['autosomalSTR', 'sexLinkedXSTR', 'sexLinkedYSTR', 'sexLinkedXInDels', 'sexLinkedYInDels', 'mitochondrialSNPs', 'mitochondrialInDels']
    if locus_type_input != 'All':
        locus_types = [locus_type_input]

    locus_results = {}

    # Loop through the list of FORSTAT objects
    for sample_object in data_list:
        if sample_object.analyseObject:
            population_group = sample_object.population_group  # Get the population group of the current FORSTAT object
            if population_group not in locus_genotype_counts:
                locus_genotype_counts[population_group] = {}
                sample_counts[population_group] = 0

            sample_counts[population_group] += 1

            for locus_type in locus_types:
                loci = getattr(sample_object, locus_type)
                for locus, alleles in loci.items():
                    if locus not in locus_genotype_counts[population_group]:
                        locus_genotype_counts[population_group][locus] = {}

                    locus_genotype = tuple(sorted(alleles))
                    if locus_genotype not in locus_genotype_counts[population_group][locus]:
                        locus_genotype_counts[population_group][locus][locus_genotype] = 0

                    locus_genotype_counts[population_group][locus][locus_genotype] += 1

    for population_group, loci in locus_genotype_counts.items():
        for locus, genotypes in loci.items():
            genotype_frequencies = {}
            for locus_genotype, count in genotypes.items():
                genotype_frequency = count / sample_counts[population_group]
                genotype_frequencies[locus_genotype] = genotype_frequency

            match_probability = sum(genotype_frequency ** 2 for genotype_frequency in genotype_frequencies.values())
            power_of_discrimination = 1 - match_probability

            if population_group not in locus_results:
                locus_results[population_group] = {}

            locus_results[population_group][locus] = {
                'Match Probability': match_probability,
                'Power of Discrimination': power_of_discrimination
            }

    return locus_results

@staticmethod
def calculate_polymorphic_information_criteria(match_probabilities, allele_frequencies):
    for population_group, loci in match_probabilities.items():
        for locus, results in loci.items():
            if population_group in allele_frequencies and locus in allele_frequencies[population_group]:
                allele_frequency = allele_frequencies[population_group][locus]
                num_alleles = len(allele_frequency)

                if num_alleles > 1:
                    #This will square every frequency value in this locus per population and sum all the values
                    sum_squared_frequencies = sum(frequency ** 2 for frequency in allele_frequency.values())
                    sum_fourth_power_frequencies = sum(frequency ** 4 for frequency in allele_frequency.values())
                    pic = 1 - sum_squared_frequencies + sum_fourth_power_frequencies - (sum_squared_frequencies ** 2)
                else:
                    pic = 0

                results['Polymorphic Information Content'] = pic

    return match_probabilities

@staticmethod
def calculate_power_of_exclusion(match_probabilities, homo_hetero_frequencies):
    for population_group, loci in match_probabilities.items():
        for locus, results in loci.items():
            homozygote_frequency = homo_hetero_frequencies[population_group][locus]['homozygote_frequency']
            heterozygote_frequency = homo_hetero_frequencies[population_group][locus]['heterozygote_frequency']
            pe = (heterozygote_frequency ** 2) * (1 - (2 * ( heterozygote_frequency * homozygote_frequency ** 2)))
            results['Power of Exclusion'] = pe

    return match_probabilities

@staticmethod
def calculate_paternity_index(match_probabilities, homozygosity):
    for population_group, loci in match_probabilities.items():
        for locus, results in loci.items():
            if population_group in homozygosity and locus in homozygosity[population_group]:
                homozygote_frequency = homozygosity[population_group][locus]['homozygote_frequency']
                pi = 1 / (2 * homozygote_frequency)
                results['Paternity Index'] = pi

    return match_probabilities

@staticmethod
def calculate_combined_statistics(match_probabilities):
    combined_statistics = {}

    for population_group, loci in match_probabilities.items():
        cmp = 1.0
        cpd = 1.0
        cpe = 1.0
        cpi = 1.0

        for locus, results in loci.items():
            match_probability = results['Match Probability']
            power_discrimination = results['Power of Discrimination']
            power_exclusion = results['Power of Exclusion']
            paternity_index = results['Paternity Index']

            cmp *= match_probability
            cpd *= (1 - power_discrimination)
            cpe *= (1 - power_exclusion)
            cpi *= paternity_index

        combined_statistics[population_group] = {
            'Combined Match Probability': cmp,
            'Combined Power of Discrimination': 1 - cpd,
            'Combined Power of Exclusion': 1 - cpe,
            'Combined Paternity Index': cpi
        }

    return combined_statistics

@staticmethod
def combine_outputs(homo_hetero_outputs, match_probabilities, unique_alleles):
    combined_outputs = {}

    # Iterate over population groups
    for population_group in homo_hetero_outputs:
        combined_outputs[population_group] = {}

        # Iterate over loci
        for locus in homo_hetero_outputs[population_group]:
            combined_outputs[population_group][locus] = {}

            # Combine homozygote and heterozygote frequencies
            combined_outputs[population_group][locus]['Unique Alleles'] = unique_alleles[population_group][locus]
            combined_outputs[population_group][locus]['Homozygote frequency'] = homo_hetero_outputs[population_group][locus]['homozygote_frequency']
            combined_outputs[population_group][locus]['Heterozygote frequency'] = homo_hetero_outputs[population_group][locus]['heterozygote_frequency']
            combined_outputs[population_group][locus]['Expected heterozygote frequency'] = homo_hetero_outputs[population_group][locus]['expected_heterozygosity']
            combined_outputs[population_group][locus]['Match Probability'] = match_probabilities[population_group][locus]['Match Probability']
            combined_outputs[population_group][locus]['Power of Discrimination'] = match_probabilities[population_group][locus]['Power of Discrimination']
            combined_outputs[population_group][locus]['Polymorphic Information Content'] = match_probabilities[population_group][locus]['Polymorphic Information Content']
            combined_outputs[population_group][locus]['Power of Exclusion'] = match_probabilities[population_group][locus]['Power of Exclusion']
            combined_outputs[population_group][locus]['Paternity Index'] = match_probabilities[population_group][locus]['Paternity Index']

    return combined_outputs

@staticmethod
def unique_alleles_per_locus(forensic_kits):
    unique_alleles = {}
    unique_allele_counts = {}

    for sample_object in forensic_kits:  # iterate over the list of FORSTAT objects
        if sample_object.analyseObject == True:
            population_group = sample_object.population_group

            if population_group not in unique_alleles:
                unique_alleles[population_group] = {}
                unique_allele_counts[population_group] = {}

            for locus, alleles in sample_object.autosomalSTR.items():
                if locus not in unique_alleles[population_group]:
                    unique_alleles[population_group][locus] = set()

                unique_alleles[population_group][locus].update(alleles)

            # Convert the sets of alleles to counts of unique alleles
            for locus in unique_alleles[population_group]:
                unique_allele_counts[population_group][locus] = len(unique_alleles[population_group][locus])
    return unique_allele_counts

@staticmethod
def find_common_genotypes(datainput, kits_metadata):
    # Generate a hashable form of genotype for each sample.
    # The genotype is considered as a tuple of sorted alleles from all loci.
    sample_genotypes = defaultdict(list)
    for sample in datainput:
        total_genotype = []
        for locus_type, loci in kits_metadata[sample.kit].items():
            for locus in loci:
                total_genotype.append(tuple(sorted(getattr(sample, locus_type).get(locus, []))))
        # Create a hashable form of total_genotype
        hashable_genotype = tuple(total_genotype)
        sample_genotypes[hashable_genotype].append((sample.sample_number, sample.population_group))

    # Check for genotypes that appear more than once
    duplicate_genotypes = {genotype: samples for genotype, samples in sample_genotypes.items() if len(samples) > 1}

    return duplicate_genotypes

@staticmethod
def export_metrics_to_excel(metrics, population_groups, shared_genotypes):
    for population_group in population_groups:
        with pd.ExcelWriter(f"{population_group}.xlsx", engine='xlsxwriter') as writer:
            workbook = writer.book

            # Export allele frequencies
            allele_frequency = metrics['allele_frequency'][population_group]
            ordered_allele_frequency = {
                locus: {allele: allele_frequency[locus][allele] for allele in sorted(allele_frequency[locus])} for locus
                in allele_frequency}
            allele_frequency_df = pd.DataFrame(ordered_allele_frequency)

            # Reset index to convert alleles to a separate column
            allele_frequency_df = allele_frequency_df.reset_index()

            # Set the first column name as 'Allele'
            allele_frequency_df.columns = allele_frequency_df.columns.astype(str)
            allele_frequency_df.columns.values[0] = 'Allele'

            # Sort the dataframe by the 'Allele' column
            sorted_allele_frequency_df = allele_frequency_df.sort_values(by='Allele', ascending=True)

            # Write the sorted allele frequency dataframe to Excel
            sorted_allele_frequency_df.to_excel(writer, sheet_name='Allele Frequency', index=False, startcol=0)

            # Export additional metrics
            additional_metrics = metrics['combined_metrics'][population_group]
            metrics_data = {}
            for metric_type, metric_data in additional_metrics.items():
                # Directly convert metric_data dictionary to DataFrame and set index
                metric_df = pd.DataFrame(metric_data, index=[metric_type])
                metrics_data[metric_type] = metric_df

            combined_metrics_df = pd.concat(metrics_data)

            # Reset the second level of index
            combined_metrics_df.reset_index(level=1, inplace=True, drop=True)

            # Rename index column to 'Locus type'
            combined_metrics_df.index.name = 'Locus type'

            # Create a new DataFrame with the combined statistics
            combined_statistics = pd.DataFrame({
                'Locus type': ['Combined'],
                'Unique Alleles': ['NA'],
                'Homozygote frequency': ['NA'],
                'Heterozygote frequency': ['NA'],
                'Expected heterozygote frequency': ['NA'],
                'Match Probability': metrics['combined_statistics'][population_group]['Combined Match Probability'],
                'Power of Discrimination': metrics['combined_statistics'][population_group]['Combined Power of Discrimination'],
                'Polymorphic Information Content': ['NA'],
                'Power of Exclusion': metrics['combined_statistics'][population_group]['Combined Power of Exclusion'],
                'Paternity Index': metrics['combined_statistics'][population_group]['Combined Paternity Index']
            }).set_index('Locus type')

            # Append the new DataFrame to the existing one
            combined_metrics_df = pd.concat([combined_metrics_df, combined_statistics])

            # Write the combined metrics dataframe to Excel
            combined_metrics_df.to_excel(writer, sheet_name='metrics', index=True, startcol=0)

            shared_genotypes_df = pd.DataFrame(
                [(genotype, sample_name, sample_population_group)
                 for genotype, samples in shared_genotypes.items()
                 for sample_name, sample_population_group in samples],
                columns=['Sample Name', 'Population Group', 'Genotype']
            )

            # Write the shared genotypes dataframe to Excel
            shared_genotypes_df.to_excel(writer, sheet_name='Shared Genotypes', index=False, startcol=0)

@staticmethod
def autosomalSTR_population_forensic_metrics(datainput):
    #Ask the user what type of analysis they will perform
    # locus_type_input = input(
    #     "Enter the locus type (autosomalSTR, sexLinkedXSTR, sexLinkedYSTR, sexLinkedXInDels, sexLinkedYInDels, mitochondrialSNPs, mitochondrialInDels, All): ")
    locus_type_input = "autosomalSTR"
    unique_alleles_locus = unique_alleles_per_locus(datainput)

    shared_genotypes = find_common_genotypes(datainput, kits_metadata)
    
    allele_frequency = calculate_allele_frequencies(data_list=datainput,
                                                    locus_type_input=locus_type_input
                                                    )
    # Ask the user for the locus type

    homo_hetero_outputs = calculate_homozygote_heterozygote_frequencies(data_list=datainput,
                                                                        locus_type_input=locus_type_input
                                                                        )
    match_probabilities = calculate_match_probability_and_power_of_discrimination(data_list=datainput,
                                                                                  locus_type_input=locus_type_input
                                                                                  )

    #This code takes the match probability and allele frequencies and calculates the PIC then combines it to the matchprobability dictionery
    calculate_polymorphic_information_criteria(match_probabilities=match_probabilities,
                                               allele_frequencies=allele_frequency
                                               )

    match_probabilities = calculate_power_of_exclusion(homo_hetero_frequencies=homo_hetero_outputs,
                                                       match_probabilities=match_probabilities
                                                       )

    match_probabilities = calculate_paternity_index(match_probabilities=match_probabilities,
                                                        homozygosity=homo_hetero_outputs
                                                    )

    combined_statistics = calculate_combined_statistics(match_probabilities=match_probabilities)
    combined_metrics = combine_outputs(homo_hetero_outputs,
                                       match_probabilities,
                                       unique_alleles_locus
                                       )

    metrics = {
        'allele_frequency': allele_frequency,
        'combined_metrics': combined_metrics,
        'combined_statistics': combined_statistics
    }

    population_groups = {sample_object.population_group for sample_object in datainput}  # Get the list of population groups from the readdatafile
    export_metrics_to_excel(metrics, population_groups, shared_genotypes)


