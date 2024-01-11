import argparse
from Bio.Blast import NCBIXML

##function that take a blast xml file as input-- ans parses it to extract the best hit
###it iterates through the blast records finding the lowest e-val so best hit 
def parse_blast_xml(blast_xml_file):
    results = {}
    with open(blast_xml_file, "r") as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        for record in blast_records:
            query_name = record.query
            best_hit = None
            best_evalue = float("inf")
            
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < best_evalue:
                        best_hit = alignment.hit_def
                        best_evalue = hsp.expect
            
            if best_hit is not None:
                results[query_name] = best_hit
    
    return results
#the rguments are two dictionaries  one for each species (in this case 2) it can be extended too more
###compare reciprocal matches as sugested (best hits) between the two species 
def find_reciprocal_best_matches(species1_results, species2_results):
    # if a query seq in spec 1 have a best hit in 2 spec and viceversa, they are considered reciprocal best match 
    reciprocal_best_matches = {}
    for species1_protein, species2_protein in species1_results.items():
        if species2_results.get(species2_protein) == species1_protein:
            reciprocal_best_matches[species1_protein] = species2_protein
    return reciprocal_best_matches 



def main():
    
    parser = argparse.ArgumentParser(description='Identify reciprocal best matches between two species based on BLAST results.')
    parser.add_argument('--species1', required=True, help='Path to the BLAST XML file for species 1')
    parser.add_argument('--species2', required=True, help='Path to the BLAST XML file for species 2')
    parser.add_argument('--output', default='homologs.txt', help='Path to the output file (default: homologs.txt)')

    args = parser.parse_args()

    species1_xml_file = args.species1 #spec1_file = path to the blast xml file for spec 1
    species2_xml_file = args.species2 #spec2_file = path to the blast xml file for spec 2

    species1_results = parse_blast_xml(species1_xml_file) #calling 'parse_blast_xml' function to parse the blast results for spec 1
    species2_results = parse_blast_xml(species2_xml_file)#calling 'parse_blast_xml' function to parse the blast results for spec 2
#calling reciprocal_best_matches
    reciprocal_best_matches = find_reciprocal_best_matches(species1_results, species2_results) 

    with open(args.output, 'w') as output_file:
        for species1_protein, species2_protein in reciprocal_best_matches.items():
            output_file.write(f"{species1_protein}\t{species2_protein}\n")

if __name__ == "__main__":
    main()


