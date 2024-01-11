"""
This script is a sequence and annotation utility script.
"""

# TODO: 1. complete assert_DNA
# TODO: 2. complete read_DNA_fasta
# TODO: 3. raise assertion error if this script is executed from command line.

# raise assertion error if this script is executed from command line.


def assert_DNA(seq: str):
    
    valid_bases = {'A', 'T', 'C', 'G'}# defining the valid bases on a DNA sequence
    
  
    if all(base in valid_bases for base in seq): #checking if all characters in hte sequence are like the defined valid bases as the above. If this is true then we return None as we se below
        return None #the sequence is valid so no error
    else: #else then the sequence is invalid so we raise an assertion error that the sequence is not a DNA sequence 
        raise AssertionError("Input sequence is not from DNA or valid at all")


def read_DNA_fasta(file_name: str) -> dict[str, str]:
    """
    Read the FASTA file and return a dict: {sequence_name: sequence}
    Raise a ValueError or AssertionError if a sequence is not DNA (handled by assert_DNA()).
    You can call my read_fasta function below to read a FASTA file.
    """
    #using the defined function above to check for false DNA sequences
    #same function
    def assert_DNA(seq: str):
        valid_bases = {'A', 'T', 'C', 'G'}
        if all(base in valid_bases for base in seq):
            return None
        else:
            raise AssertionError("Input sequence is not a valid DNA sequence")

    # Use the read_fasta function provided to read the FASTA file
    fasta_data = read_fasta(file_name)

    # Check if each sequence is DNA and store it in a new dictionary
    dna_sequences = {}
    for name, sequence in fasta_data.items():
        try:
            assert_DNA(sequence)
            dna_sequences[name] = sequence
        except AssertionError:
            raise AssertionError(f"Sequence '{name}' in the FASTA file is false or not DNA at all")

    return dna_sequences

#region functions provided by me, you do not need to do anything here.
def assert_3divisible(seq: str):
    """
    Check if this sequence is a multiple of 3.

    Raise a ValueError or AssertionError if seq is not divisible by 3.
    """

    result = len(seq) % 3 == 0

    # raise ValueError if seq is not divisible by 3.
    # if not result:
    #     raise ValueError("Invalid sequence, input sequence can only be multiple of three.")
    
    # raise AssertionError if seq is not divisible by 3.
    assert result, "Invalid sequence, input sequence can only be multiple of three."

def read_annotation(file_name: str) -> dict[str, dict]:
    """
    read the gene annotation and return a nested dict.
    {
        gene_name : {
            "chrom" : str, # chromosome of that gene
            "strand" : str, # + or -,
            "txStart" : int, # Transcription start position (or end position for minus strand item)
            "txEnd" : int,   # Transcription end position (or start position for minus strand item)
            "exonCount" : int, # 1, 2
            "exonStarts" : int tuple, # (465, 4969)
            "exonEnds" : int tuple # (465, 4969)
        }
    }

    If file not exists, it will raise FileNotFoundError.
    """

    file = open(file_name) # this may raise a FileNotFoundError

    # skip the first line (column names)
    next(file)

    sac_genes = {}

    # loop through remaining lines
    for row in file.readlines():
        gene = row[:-1].split("\t")
        sac_genes[gene[1]] = {
            "chrom" : gene[2],
            "strand" : gene[3],
            "txStart" : int(gene[4]),
            "txEnd" : int(gene[5]),
            "exonCount" : int(gene[8]),
            "exonStarts" : tuple(int(i) for i in gene[9].split(",")[:-1]), # use list comprehension to split the 
            # there is an additional comma at the end, [:-1] to remove the last item
            "exonEnds" : tuple(int(i) for i in gene[10].split(",")[:-1])
        }

    file.close()
    return sac_genes

def read_fasta(file_name: str) -> dict[str, str]:
    """
    read the fasta file and return a dict.
    {
        sequence_name : sequence
    }

    If file not exists, it will raise FileNotFoundError.
    """

    file = open(file_name)
    
    contents = file.read()

    # split contents by >, because each sequence starts with >
    entries = contents.split('>')[1:]

    # >name
    # sequence
    # We use partition to separate name (1st element) and sequence (3rd element)
    partitioned_entries = [entry.partition('\n') for entry in entries]
    # use dict comprehension to get a dict of each chromosome
    return { entry[0] : entry[2].replace("\n", "") for entry in partitioned_entries}

def write_fasta(file_name: str, sequence_dict: dict[str, str]):
    """
    save sequence_dict to a fasta file.
    
    > dict_key
    dict_value
    """

    file = open(file_name, "w")

    for name, seq in sequence_dict.items():
        file.write(f">{name}\n")

        # add line break if each line is more than 60 characters
        seq_len =  len(seq)
        for i in range(0, seq_len, 60):
            if i + 60 > seq_len:
                file.write(seq[i:seq_len])
            else:
                file.write(seq[i:i+60])
            file.write("\n")
    
    file.close()

def get_pam(seq: str) -> int | None:
    """
    Get the first PAM index of the DNA sequence.

    Return None if there is no PAM.

    Raise a ValueError or AssertionError if seq is not DNA. (This should be handled by assert_DNA)
    """
    
    assert_DNA(seq)

    # remove the first nt and find GG.
    # so now, the index is the index of N of NGG
    index = seq[1:].find("GG")

    if index == -1:
        # print("No PAM in provided sequence")
        return None
    else:
        return index

def get_upstream_indices(gene_annotation: dict, upstream: int) -> tuple[int, int]:
    """
    Get the upstream chromosome start and end indices using gene_annotation.

    This means you can use chr_seq[start:end] to slice upstream sequence of reference genome.

    If you want to get upstream sequence of gene on minus strand, you still need to use reverse_complement.
    """
    if gene_annotation["strand"] == "+":
        chr_end = gene_annotation["txStart"] # ends at txStart
        chr_start = chr_end - upstream       # starts at txStart - upstream
    else:
        chr_start = gene_annotation["txEnd"] # starts at txEnd
        chr_end = chr_start + upstream       # ends at txEnd + upstream
    
    return (chr_start, chr_end)

_BASE_PAIRS = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} # _ for hidden
"""
constant dict of DNA pairs (should not be changed)
"""

def reverse_complement(seq: str) -> str:
    """
    Get the reverse complement sequence.

    Raise a ValueError or AssertionError if seq is not DNA. (This should be handled by assert_DNA)
    """

    assert_DNA(seq) # It will raise a ValueError or AssertionError if seq is not DNA.

    # get complement nucleotide of each one
    # [::-1] reverse the seq string
    complement_nts = [ _BASE_PAIRS[nt] for nt in seq[::-1] ]
    return ''.join(complement_nts) # join the list to a str
#endregion
