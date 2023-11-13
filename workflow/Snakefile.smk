from Bio import AlignIO
from collections import Counter

configfile: "config/config.yaml"

# the 4 amino acid types

AA_TYPES = {
    "K": "pos",
    "R": "pos",
    "H": "pos",
    "D": "neg",
    "E": "neg",
    "S": "pol",
    "T": "pol",
    "Y": "pol",
    "N": "pol",
    "U": "pol",
    "G": "npol",
    "A": "npol",
    "V": "npol",
    "C": "npol",
    "P": "npol",
    "L": "npol",
    "I": "npol",
    "M": "npol",
    "F": "npol",
    "W": "npol",
    "-": "gap"
}

POSITIONS = list(config["SNP"])


rule all:
    input: expand("results/out{SNP}.txt", SNP = range(0, len(config["SNP"]))) # parallelized

rule make_msa:
    input: config["multifasta"]
    output: config["msa"]
    params: config["clustalo"]
    message: "Making {output} file from {input} using {params[0]}"
    shell: "{params[0]} -i {input} -o {output}"

rule read_msa:
    input: rules.make_msa.output
    output: "resources/position{SNP}.txt"
    run:
        alignment = AlignIO.read(open(input[0]), "fasta")
        with open(output[0], "w") as FILE:
            for rec in alignment:
                print(rec.seq[POSITIONS[int(wildcards.SNP)] -1 ], file=FILE) # write away the amino acid of the snp position



rule evaluate:
    input: rules.read_msa.output
    output: "results/out{SNP}.txt"
    run:
        def parse():
            with open(input[0], "r") as FILE:
                x = [line.strip() for line in FILE]
                y = [AA_TYPES[aa] for aa in x]
                count = Counter(x)
                type = Counter(y)
            return type, count


        def determine_conservation(type, count, size):
            type_conserve = max(type.values())/size
            aa_conserve = max(count.values())/size
            return type_conserve, aa_conserve


        def score(mutated_aa):
            mutated_aa_type = AA_TYPES[mutated_aa]
            try:
                if counts[mutated_aa] == max(counts.values()): # if most prevalent aa in position
                    return 1 # likely not even an SNP
                else:
                    if a_con >= 0.8: # if not most prevalent aa in a conserved position
                        s = 3 + (1 - a_con) * 5 # heavy penalty
                        if t_con >= 0.9: # if not most prevalent type in conserved type position
                            return 3 + (1 - t_con) * 5 + s # heavy penalty
                        else: # if not most prevalent type in non conserved position
                            return s + (1 - t_con) * 5 # weak penalty
                    else: # if not most prevalent aa in a non conserved position
                        s = (1 - a_con) * 5 # weak penalty
                        if t_con >= 0.9: # if not most prevalent type in a type conserved position
                            return 2 + (1 - t_con) * 5 + s # heavy penalty
                        else: # if not most prevalent type in a type non conserved position
                            return s + (1 - t_con) * 5 # weak penalty

            except KeyError: # if aa not in the comparison set then it is likely more deleterious
                try:
                    if types[mutated_aa_type] == max(types.values()): # if type is most prevalent in position
                        return 3 # not same aa penalty
                    else:
                        if t_con >= 0.9: # if not most prevalent type in conserved type position
                            return 5 + (1 - t_con) * 10 # heavy penalty
                        else:
                            return 3 + (1 - t_con) * 10 # weaker penalty
                except KeyError: # if type also not in the comparison set then it is most deleterious
                    return 10

        def write():
            with open(output[0], "w") as FILE:
                print(f"""SNP: {mutated_aa} at position {POSITIONS[int(wildcards.SNP)]}
has been rated as {rating}
Amino acid conservation rate is {a_con}
Amino acid type conservation rate is {t_con}""", file=FILE)


        mutated_aa = config["SNP"][POSITIONS[int(wildcards.SNP)]] # get the snp amino acid
        types, counts = parse()
        n_samples = sum(counts.values())
        t_con, a_con = determine_conservation(types, counts, n_samples)
        rating = score(mutated_aa)
        write()

