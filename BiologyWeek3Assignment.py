# full codon table
codon_table = {
    'AUG': ('Methionine', 'M'), 'UUU': ('Phenylalanine', 'F'), 'UUC': ('Phenylalanine', 'F'),
    'UUA': ('Leucine', 'L'), 'UUG': ('Leucine', 'L'), 'UCU': ('Serine', 'S'), 'UCC': ('Serine', 'S'),
    'UCA': ('Serine', 'S'), 'UCG': ('Serine', 'S'), 'UAU': ('Tyrosine', 'Y'), 'UAC': ('Tyrosine', 'Y'),
    'UGU': ('Cysteine', 'C'), 'UGC': ('Cysteine', 'C'), 'UGG': ('Tryptophan', 'W'), 'CUU': ('Leucine', 'L'),
    'CUC': ('Leucine', 'L'), 'CUA': ('Leucine', 'L'), 'CUG': ('Leucine', 'L'), 'CCU': ('Proline', 'P'),
    'CCC': ('Proline', 'P'), 'CCA': ('Proline', 'P'), 'CCG': ('Proline', 'P'), 'CAU': ('Histidine', 'H'),
    'CAC': ('Histidine', 'H'), 'CAA': ('Glutamine', 'Q'), 'CAG': ('Glutamine', 'Q'), 'CGU': ('Arginine', 'R'),
    'CGC': ('Arginine', 'R'), 'CGA': ('Arginine', 'R'), 'CGG': ('Arginine', 'R'), 'AUU': ('Isoleucine', 'I'),
    'AUC': ('Isoleucine', 'I'), 'AUA': ('Isoleucine', 'I'), 'GUU': ('Valine', 'V'), 'GUC': ('Valine', 'V'),
    'GUA': ('Valine', 'V'), 'GUG': ('Valine', 'V'), 'GCU': ('Alanine', 'A'), 'GCC': ('Alanine', 'A'),
    'GCA': ('Alanine', 'A'), 'GCG': ('Alanine', 'A'), 'GAU': ('Aspartic acid', 'D'), 'GAC': ('Aspartic acid', 'D'),
    'GAA': ('Glutamic acid', 'E'), 'GAG': ('Glutamic acid', 'E'), 'GGU': ('Glycine', 'G'), 'GGC': ('Glycine', 'G'),
    'GGA': ('Glycine', 'G'), 'GGG': ('Glycine', 'G'), 'UAA': ('Stop', None), 'UAG': ('Stop', None),
    'UGA': ('Stop', None), 'ACU': ('Threonine', 'T'), 'ACC': ('Threonine', 'T'), 'ACA': ('Threonine', 'T'),
    'ACG': ('Threonine', 'T'), 'AAU': ('Asparagine', 'N'), 'AAC': ('Asparagine', 'N'), 'AAA': ('Lysine', 'K'),
    'AAG': ('Lysine', 'K'), 'AGU': ('Serine', 'S'), 'AGC': ('Serine', 'S'), 'AGA': ('Arginine', 'R'),
    'AGG': ('Arginine', 'R')
}


# create a dictionary for the letter of amino acids, and all possible codons for that acid
def create_amino_acid_to_codon_dict():
    amino_acid_to_codon_dict = {}
    for codon, (name, letter) in codon_table.items():
        if name:
            if letter not in amino_acid_to_codon_dict:
                amino_acid_to_codon_dict[letter] = []
            amino_acid_to_codon_dict[letter].append(codon)

    return amino_acid_to_codon_dict


# call the function
amino_acid_to_codon = create_amino_acid_to_codon_dict()


# function for transcribing
def transcribe(dna_bottom_strand):
    transcription_map = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}

    # reverse complement the bottom strand, but U and not T because RNA
    mrna = ''.join([transcription_map[nucleotide] for nucleotide in dna_bottom_strand])

    # following print statement from instructions
    print("Complement:", mrna.replace('U', 'T'))

    # print out mrna
    print("Transcribed mRNA sequence:", mrna)

    return mrna


# function for transcribing
def translate(mRNA):
    # empty array to contain the acids
    protein = []

    # loop through each 3 RNA
    for index in range(0, len(mRNA), 3):
        # get the current three
        codon = mRNA[index:index + 3]
        # find acid
        amino_acid = codon_table.get(codon, None)

        # stop translation
        if amino_acid is None or amino_acid[0] == 'Stop':
            break

        # add acid to array
        if amino_acid:
            protein.append(str((amino_acid[0] + " (" + amino_acid[1] + ") ")))  # append acid

    # return array
    return protein


# function to call the previoujs two together
def transcribe_and_translate(dna_seq):
    # checks if dna is valid
    if not set(dna_seq.upper()).issubset({'A', 'T', 'C', 'G'}):
        print("DNA sequence is not valid. Contains nucleotides that aren't 'A', 'T', 'C', nor 'G'.\n")
        return

    # transcribe
    mrna_seq = transcribe(dna_seq.upper())

    # translate
    acids = translate(mrna_seq)

    # print results
    print("Aminoacid =", '- '.join(acids))
    print("")


# recursive function to find all codon combinations for up to 3 amino acids
def all_codon_combinations(all_amino_acids, index=0, current_combination='', current_frequencies=None):
    # create a new dict of codon frequencies if it doesn't exist yet
    if current_frequencies is None:
        current_frequencies = {}

    # base case, where all amino acids have been checked
    if index == len(all_amino_acids):
        return [current_combination], [current_frequencies]

    # get current amino acid
    current_amino_acid = all_amino_acids[index]

    # check if acid doesn't exist
    if current_amino_acid not in amino_acid_to_codon:
        print("Invalid amino acid:", current_amino_acid)
        return [], []

    # get all possible codons for this amino acid
    possible_codons = amino_acid_to_codon.get(current_amino_acid, [])

    # list to store all combinations
    all_combinations = []
    all_frequencies = []

    # loop through each possible codon
    for codon in possible_codons:
        # create copy of frequency dict so that everything doesn't alter the original dict
        next_frequencies = current_frequencies.copy()

        # if the frequency of the current codon doesn't exist in the current frequency dict, set it to 0
        if codon not in next_frequencies:
            next_frequencies[codon] = 0

        # increment frequency
        next_frequencies[codon] += 1

        # recursive call to get all combinations for the next amino acid
        next_combinations, next_frequencies = all_codon_combinations(all_amino_acids, index + 1,
                                                                     current_combination + codon, next_frequencies)

        # put all entries to the original arrays
        all_combinations.extend(next_combinations)
        all_frequencies.extend(next_frequencies)

    return all_combinations, all_frequencies


# main loop
while True:
    # choosing which functionality to use
    print("1. Transcribe DNA to mRNA, then translate mRNA to amino acid sequence")
    print("2. List all mRNA combinations for amino acid sequence, along with the codon frequencies")
    print("3. Exit")
    choice = input("Enter 1, 2 or 3: ")
    while choice != '1' and choice != '2' and choice != '3':
        print("Invalid response\n")
        choice = input("Enter 1, 2 or 3 ")

    # transcribe + translate
    if choice == '1':
        DNA = input("\nEnter DNA sequence (top / bottom strand): ")
        # make sure dna is multiple of 3
        if len(DNA) % 3 != 0:
            print("Invalid. Should be a multiple of 3.\n")
        else:
            transcribe_and_translate(DNA)

    # codon combinations for amino acids
    elif choice == '2':
        amino_acids = input("\nEnter amino acid sequence (max 3 amino acids):  ")
        # make sure there are at most 3 amino acids
        if len(amino_acids) > 3:
            print("Invalid. There can only be 3 amino acids at most.\n")
        else:
            combinations, frequencies = all_codon_combinations(amino_acids.upper())

            # print out all combinations and the frequency of each codon for each combination
            for i in range(len(combinations)):
                print(str(i + 1) + ".", "mRNA:", combinations[i])
                for used_codon, codon_frequency in frequencies[i].items():
                    print(used_codon, "=", codon_frequency)
                print("")
            print("")
    # exit the program
    else:
        break
