"""
This script will read gene annotation and genome sequence to get upstream PAM position.
And save the result as a text file or tsv.
"""
#importing argparse and utils modules
import argparse
import utils
import re

# regular expressions for PAM sequences in crisp-cas 9 biology
re_plus_pam = re.compile(r"^NGG$")  # NGG
re_minus_pam = re.compile(r"^CCN$")  # CCN
#principal script for parsing command-line arguments, reading annotation file and genomic data-output file
def main():
    parser = argparse.ArgumentParser(description="Find PAM upstream")

    # Defining command-line arguments based on suggestions on pdf for -h
    parser.add_argument("-a", "--annotation", required=True, help="Name of annotation file")
    parser.add_argument("-g", "--genome", required=True, help="Name of genome fasta file")
    parser.add_argument("-o", "--output", required=True, help="Name of output file")
    parser.add_argument("-t", "--output-type", default="name", help="Name or position (default = name)")
    parser.add_argument("-u", "--upstream", default=30, type=int, help="Number of base pairs upstream of the gene (default = 30)")

    args = parser.parse_args()# saving the args

    # # Checking if the output type 
    if args.output_type not in ["name", "position"]:#is valid and in format
        print("Invalid output type!")
        return

    try: #reading annotation and genomic files with function on utils
        annotation = utils.read_annotation(args.annotation)
        genome = utils.read_DNA_fasta(args.genome)
    except FileNotFoundError as e: #exception and assertio if file is not found, print Error: annotation or genome file is not 
        print(f"Error: {e}")
        return
    except AssertionError:
        print("Invalid fasta file!")
        return

    # open the file for outut and write
    with open(args.output, "w") as output_file:
        #the output file is specified by the -o argumentiterating through genes in the annotation data and checks for the presence of a PAM sequence upstream of each gene.
        #utils function is_pam_upstream (given) check if there's a pam sequence upstream of the gene and call get_pam_upstream_position function
        for gene_name, gene_annotation in annotation.items():
            if utils.is_pam_upstream(gene_annotation, genome, args.upstream):
                if args.output_type == "name":
                    output_file.write(gene_name + "\n")
                elif args.output_type == "position":
                    strand, position = utils.get_pam_upstream_position(gene_annotation, genome, args.upstream)
                    output_file.write(f"{strand}\t{position}\n")

if __name__ == "__main__": # I like this at the end of the main() for purposes of similarity to C++ where all the defined functions are at the end
    main()

 ### defined functions |
'''
'''
def is_pam_upstream(gene_annotation, genome_seq, upstream):
    # extracts the chromosome sequence for the gene based on the information provided in 'geneannotation'
    chrom_seq = genome_seq[gene_annotation["chrom"]]
    
    # Calculate the start and end positions of the upstream region based on the gene's strand
    if gene_annotation["strand"] == "+": #if the gene is on the plus strand, calculates the region upstream
        chr_end = gene_annotation["txStart"]
        chr_start = chr_end - upstream
    else:
        chr_start = gene_annotation["txEnd"]
        chr_end = chr_start + upstream

    # Extract the upstream region
    upstream_seq = chrom_seq[chr_start:chr_end]

    # Check for the presence of PAM on both strands
    if utils.get_pam(upstream_seq) is not None or utils.get_pam(utils.reverse_complement(upstream_seq)) is not None:
        return True
    else:
        return False

import utils

def get_pam_upstream_position(gene_annotation, genome_seq, upstream):
    #extracting chrom sequence for the gen
    chrom_name = gene_annotation["chrom"]
    chrom_seq = genome_seq[chrom_name]
#calculating the start and end positions of the upstream region based on the gene's strand using the get_upstrem
    start, end = utils.get_upstream_indices(gene_annotation, upstream)
    upstream_seq = chrom_seq[start:end]

    pam_plus_strand_index = utils.get_pam(upstream_seq) #checking for the presence of a PAM if no pam is found in + it will check - more below
    if pam_plus_strand_index is not None:
        return ("+", start + pam_plus_strand_index)

    upstream_minus_seq = utils.reverse_complement(upstream_seq)
    pam_minus_strand_index = utils.get_pam(upstream_minus_seq)
    if pam_minus_strand_index is not None:
        return ("-", end - pam_minus_strand_index - 3)

    return None

    


    