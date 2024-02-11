# Sliding-window-MM

## Script to find conserved regions among DNA sequences of a target and remaining taxa present in an alignment. It calculates mismatch statistics (min, max, sd) within user defined sliding windows across the alignment. It was developed to guide primer design for species detection using a qPCR approach.

### How to use it:
The script takes four positional parameters:
1) Alignment file in a fasta format
2) Name of the target species (under quotation)
3) Sliding window size
4) Prefix to save output files 

The alignment should be in the fasta format with sequence titles separated by spaces in a way that species names would be in the second and third position:
>OP689439.1 Etheostoma vitreum isolate EvitH cytochrome b (cytb) gene, partial cds; mitochondrial

#### Code example:
   python .\Get_differentiating_regions.py Alignment_file.fasta "Species name" 23 Output_name

### Output files
The script saves two output files:
1) Fasta alignment with sequences sorted per species
2) Tab separated text file where each line is the alignment position and columns have information of the sequence motif in the target species and mismatch statistics per species
