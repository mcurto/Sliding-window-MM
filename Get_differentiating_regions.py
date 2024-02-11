#Script calculating mismtaches statistics for each sliding window of a given size across an algnement
#Mismatches are calculated between a user defined target species (present in the alignement) and all remaining sequences
#Alignement should be in the fasta format with sequence titles with fileds sperated by spaces with species names in the second and third position:
#Example: >OP689439.1 Etheostoma vitreum isolate EvitH cytochrome b (cytb) gene, partial cds; mitochondrial


import sys
import numpy as np

#Read alignement file and save information into a dictionary where:
#keys = sequence titles
#items = sequences
def parse_fasta(fasta):
    result = {}
    with open(fasta) as f:
        for l in f:
            led = l.rstrip("\n")
            if l.startswith(">"):
                name = l
                result[name] = ""
            else:
                result[name] = result[name] + led.upper()
    return result

#Go though each position of the alugnement and count the number of sequences containg gap simbols (-);
#If all sequences contain gaps sympbols, remove that position
def remove_gap_only_positions(fasta_parsed, alignement_length):
    nr_seqs = len(fasta_parsed)
    result = {}
    for i in range(alignement_length):
        temp_positions = {}
        count_gap = 0
        for name, seq in fasta_parsed.items():
            nucl = seq[i]
            temp_positions[name] = nucl
            if nucl == "-":
                count_gap += 1
        if count_gap < nr_seqs:
            for name, seq in fasta_parsed.items():
                result[name] = result.get(name, "") + seq[i]
    return result

#Join all sequences from the same species into the same key in the dictionary:
#Key = sepecies name
#Items = list of sequences
#Example: {species1 : [seq1, seq4, seq24], species2: [seq5, seq,7], ...}
#Species name information are in the second and third fields of the title name seperated by spaces
def sort_by_name(sequences_all):
    result = {}
    for name, seq in sequences_all.items():
        species = " ".join(name.split()[1:3])
        result[species] = result.get(species, []) + [seq]
    return result

#Write fasta file sorted based on a species list
# It takes the sorted by species dictionry data format
def write_fasta_sorted(sequences_per_species, out_fasta):
    with open(out_fasta, "w") as out:
        species_list = list(sequences_per_species.keys())
        species_list.sort()
        for species in species_list:
            data = sequences_per_species[species]
            for seq in data:
                out.write(">{}\n{}\n".format(species, seq))

#Get number of positions in an alignement
def get_alignement_length(sequences_all):
    return len(list(sequences_all.values())[0])

#Partition the data for into sliding windows for a specified species saving only windows with no variation into a dictionary where
#Keys = Alignemt position
#Items = list of: list of motifs across all sequences, percentage of missing data, number of sequences
def get_equal_sliding_windows(sequences_per_species, target_species, alignement_length, window_size):
    sequence_partitions_equal = {}
    target_species_sequences = sequences_per_species[target_species]
    nr_seq = len(target_species)
    for i in range(0, alignement_length):
        missing = 0
        seq_partition = []
        for seq in target_species_sequences:
            motif = seq[i: i + window_size]
            nr_missing = motif.count("-") + motif.count("N")
            if nr_missing / window_size < 1:
                seq_partition.append(motif)
            else:
                missing += 1
        if len(set(seq_partition)) == 1:
            sequence_partitions_equal[i] = [seq_partition[0], missing / len(seq_partition), nr_seq]
    return sequence_partitions_equal

#Calculate numeber of mismactes between two sequences
def mm(seq1, seq2):
    mm = 0
    for i in range(min(len(seq1), len(seq2))):
        if seq1[i] != seq2[i]:
            mm += 1
    return mm

#Calculate number of mismatches between the sliding windows of the target secies and the other species
#Results are saved in a dictionary with the information of:
#Key = algnemnt position
#items = list containing: number of sequences, minimum number ofmismactes,maximum number of mismatches, sd of number of mismatches
def compare_motifs(sequence_partitions_equal, sequence_other_species, window_size):
    species_result = {}
    for i, data in sequence_partitions_equal.items():
        motif = data[0]
        mm_s = []
        for seq in sequence_other_species:
            target_motif = seq[i: i + window_size]
            nr_missing = target_motif.count("-") + target_motif.count("N")
            if nr_missing / window_size < 0.25:
                mm_s.append(mm(motif, seq[i: i + window_size]))
        n_seqs_all = len(data)
        nr_seqs_mm =  len(mm_s)
        if nr_seqs_mm > 0:
            max_mm = max(mm_s)
            min_mm = min(mm_s)
            mean_mm = np.mean(mm_s)
            sd_mm = np.std(mm_s)
            species_result[i] = [nr_seqs_mm, min_mm, max_mm, sd_mm]
        else:
            species_result[i] = ["NA"] * 4
    return species_result

#Calculate mismtatch statistics between the target species and the remaining.
#Its saves the results in a tab seperated text files where
#each line is the alignement position
#All four parameters per are saved in collumns per species
def calculate_all(sequences_per_species, target_species, out_file, alignement_length, window_size):
    sequence_partitions_equal = get_equal_sliding_windows(sequences_per_species, target_species, alignement_length, window_size)
    results_per_species = {}
    species_list = []
    for species, sequences in sequences_per_species.items():
        if species != target_species:
            species_list.append(species)
            results_per_species[species] = compare_motifs(sequence_partitions_equal, sequences, window_size)
    position_list = list(sequence_partitions_equal.keys())
    position_list.sort()
    species_list.sort()
    with open(out_file, "w") as out:
        out.write("Position\tMotif\t% missing\tnr. seq\t" + "\t\t\t\t".join(species_list) + "\n")
        out.write("\t"*3 + ("\t" + "\t".join(["nr. sequences", "min. mm", "max. mm", "sd. mm"])) * len(species_list) + "\n")
        for i, data in sequence_partitions_equal.items():
            to_write = [str(i)] + [str(j) for j in data]
            for species in species_list:
                species_data = results_per_species[species]
                to_write += [str(k) for k in species_data[i]]
            out.write("\t".join(to_write) + "\n")

#Main fuction requiring the folowing positional inputs:
#1) alignement file
#2) target species name
#3) window size
#4) output prefix
def main():
    print("Read alignement")
    fasta_data = parse_fasta(sys.argv[1])
    alignement_length = get_alignement_length(fasta_data)
    print("Delete gap only postions")
    fasta_data_degap = remove_gap_only_positions(fasta_data, alignement_length)
    print("sort sequences per species")
    sequences_per_species = sort_by_name(fasta_data_degap)
    target_species = sys.argv[2]
    window_size = int(sys.argv[3])
    write_fasta_sorted(sequences_per_species, sys.argv[4] + "sorted.fasta")
    print("Calculate mismatches")
    calculate_all(sequences_per_species, target_species, sys.argv[4], alignement_length, window_size)

if __name__ == "__main__":
    main()
