import argparse
import sys
import utils

def is_pam_upstream(seq: str) -> bool:
    """
    Check if there is a PAM sequence (NGG or CCN) on either strand.

    Return True if there is one, False if not.
    """
    return "NGG" in seq or "CCN" in seq

def get_pam_upstream_position(seq: str) -> tuple[str, int] | None:
    """
    Get the strand and 5' reference chromosome position of the PAM sequence.

    Return a tuple (strand, position) or None if no PAM is found.
    """
    if "NGG" in seq:
        strand = "+"
        position = seq.index("NGG")
    elif "CCN" in seq:
        strand = "-"
        position = seq.index("CCN")
    else:
        return None

    return strand, position

def main():
    parser = argparse.ArgumentParser(description="Find genes with PAM upstream or get the position of upstream PAM.")
    parser.add_argument("-a", "--annotation", required=True, help="Name of annotation file")
    parser.add_argument("-g", "--genome", required=True, help="Name of genome fasta file")
    parser.add_argument("-o", "--output", required=True, help="Name of output file")
    parser.add_argument("-t", "--output-type", default="name", help="Name or position (default = name)")
    parser.add_argument("-u", "--upstream", type=int, default=30, help="Number of base pairs upstream of the gene (default = 30)")

    args = parser.parse_args()

    if args.output_type not in ["name", "position"]:
        print("Invalid output type!")
        sys.exit(1)

    try:
        annotation_data = utils.read_annotation(args.annotation)
        dna_data = utils.read_DNA_fasta(args.genome)

        genes_with_pam_upstream = []

        for gene_name, gene_annotation in annotation_data.items():
            if gene_annotation["strand"] == "+":
                chr_seq = dna_data[gene_annotation["chrom"]]
            else:
                chr_seq = utils.reverse_complement(dna_data[gene_annotation["chrom"]])

            start, end = utils.get_upstream_indices(gene_annotation, args.upstream)
            upstream_seq = chr_seq[start:end]

            if is_pam_upstream(upstream_seq):
                if args.output_type == "name":
                    genes_with_pam_upstream.append(gene_name)
                elif args.output_type == "position":
                    pam_position = get_pam_upstream_position(upstream_seq)
                    if pam_position:
                        print(f"{gene_name} ({pam_position[0]}, {pam_position[1] + start})")

        if args.output_type == "name":
            with open(args.output, "w") as output_file:
                for gene_name in genes_with_pam_upstream:
                    output_file.write(f"{gene_name}\n")

    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except AssertionError:
        print("Invalid fasta file!")
        sys.exit(1)

if __name__ == "__main__":
    main()
