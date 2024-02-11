# Sliding-window-MM

## Script to find conserved regions among DNA sequences of a target and remaining taxa present in an alignement. It calculate mismatch statistics (min, max, sd) within user defined sliding windows accross the alignement. It was developed to gide primer design for species detection using a qPCR aproach.

### How to use it:
The script takes four positional parameters:
1) Alignement file in a fasta format
2) Name of the target species (under cootations)
3) Sliding window size
4) Prefix to save output files 

The alignement should be in the fasta format with sequence titles sperated by spaces in a way that species names would be in the second and third position:
>OP689439.1 Etheostoma vitreum isolate EvitH cytochrome b (cytb) gene, partial cds; mitochondrial

#### Code example:
   python .\Get_differentiating_regions.py Alignement_file.fasta "Species name" 23 Output_name

### Output files
The script saves two output files
