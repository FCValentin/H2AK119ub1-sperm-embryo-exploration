-*- coding: utf-8 -*-

def open_file(fichier):
    with open(fichier, 'r') as f:
        lines = f.readlines()
    # No return end of line
    return ''.join(line.strip() for line in lines[1:])

def find_commom_seq(sequence1, sequence2):
    if len(sequence1) != len(sequence2):
        raise ValueError("Sequence length are differents")

    commom_seq = ''
    for letter1, letter2 in zip(sequence1, sequence2):
        # Ignore positions where both letters are N
        if letter1 != letter2 and letter1 != 'N':
            commom_seq += letter1
        elif letter1 != letter2 and letter2 != 'N':
            commom_seq += letter2
        else:
            commom_seq += letter1  # Letters are identical or N

    return commom_seq

if __name__ == "__main__":
    
    
    import os

    # path to directory
    new_directory = 'Documents/These/Analyses/Chen2021'

    # change Path
    os.chdir(new_directory)

    for i in range(1, 23):
        if i == 20:
            value = "X"
        elif i == 21:
            value = "Y"
        elif i == 22:
            value = "MT"
        else:
            value = str(i)      
        # For each chromosome, compute the script on SNP masked genome
        file1 = 'PWK_PhJ_DBA_2J_dual_hybrid.based_on_GRCm39_N-masked/chr'+value+'.N-masked.fa'
        file2 = 'PWK_PhJ_C57BL_6NJ_dual_hybrid.based_on_GRCm39_N-masked/chr'+value+'.N-masked.fa'
        
        # Read sequences
        sequence1 = open_file(file1)
        sequence2 = open_file(file2)

        # Identify a common sequence
        commom_seq = find_commom_seq(sequence1, sequence2)
        filename = 'mask/chr'+valeur+'.N-masked.fa'

        # Export genome ('w')
        with open(filename, 'w') as file:
            file.write(commom_seq)
