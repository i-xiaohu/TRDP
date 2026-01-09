import random
import pandas as pd


def gen_repeat(motif_len: int, repeat_n: int, mutation_rate: float, flank_l: int, flank_r: int, tab: dict):
    alphabet = 'ACGT'
    text = ''
    for i in range(0, flank_l):
        text += alphabet[random.randint(0, 3)]
    # Generate a random motif with the given length
    pattern = list()
    for i in range(motif_len):
        pattern.append(alphabet[random.randint(0, 3)])
    pattern = ''.join(pattern)
    print('Pattern:', pattern)

    mutated_n = 0
    if mutation_rate < 1e-9:  # Avoid division by zero
        number_pool = -1
    else:
        number_pool = round(1 / mutation_rate)
    for i in range(repeat_n):
        result = []
        for j in range(motif_len):
            # If the picked number from the pool is 1, then a mutation happens
            if number_pool == -1 or random.randint(1, number_pool) != 1:
                result.append(pattern[j])
            else:
                mutated_n += 1
                var_type = random.randint(0, 2)
                if var_type == 0:  # Mismatch
                    v = alphabet[random.randint(0, 3)]
                    while v == pattern[j]:
                        v = alphabet[random.randint(0, 3)]
                    result.append(v)
                    # print('  Mismatch at %d: %s -> %s' % (j, pattern[j], v))
                elif var_type == 1:  # Deletion
                    pass
                    # print('  Deletion at %d' % j)
                else:  # Insertion
                    v = alphabet[random.randint(0, 3)]
                    # print('  Insertion at %d: %s' % (j, v))
                    result.append(v)
                    result.append(pattern[j])
        result = ''.join(result)
        # print('  Mutated motif: %s\n' % result)
        text += result
    for i in range(0, flank_r):
        text += alphabet[random.randint(0, 3)]
    tab['ID'].append(len(tab['ID']) + 1)
    tab['motif'].append(pattern)
    tab['periods'].append(repeat_n)
    tab['mutation'].append(mutated_n)
    tab['flank_l'].append(flank_l)
    tab['flank_r'].append(flank_r)
    tab['text'].append(text)
    return text


def gen_vntr(motif_len: int, r1: int, r2: int, var_rate, mut_rate: float, flank_l: int, flank_r: int, tab: dict):
    alphabet = 'ACGT'
    text1, text2 = '', ''
    if var_rate < 1e-9:
        var_pool = -1
    else:
        var_pool = round(1 / var_rate)

    for i in range(0, flank_l):
        c = alphabet[random.randint(0, 3)]
        text1 += c
        if var_pool == -1 or random.randint(1, var_pool) != 1:
            text2 += c
        else:
            var_type = random.randint(0, 2)
            if var_type == 0:  # Mismatch
                v = alphabet[random.randint(0, 3)]
                while v == c:
                    v = alphabet[random.randint(0, 3)]
                text2 += v
            elif var_type == 1:  # Deletion
                pass
            else:  # Insertion
                v = alphabet[random.randint(0, 3)]
                text2 += c
                text2 += v

    # Generate a random motif with the given length
    pattern = list()
    for i in range(motif_len):
        pattern.append(alphabet[random.randint(0, 3)])
    pattern = ''.join(pattern)
    print('Pattern:', pattern)
    if mut_rate < 1e-9:  # Avoid division by zero
        mut_pool = -1
    else:
        mut_pool = round(1 / mut_rate)

    mutated_n1 = 0
    for i in range(r1):
        result = []
        for j in range(motif_len):
            # If the picked number from the pool is 1, then a mutation happens
            if mut_pool == -1 or random.randint(1, mut_pool) != 1:
                result.append(pattern[j])
            else:
                mutated_n1 += 1
                var_type = random.randint(0, 2)
                if var_type == 0:  # Mismatch
                    v = alphabet[random.randint(0, 3)]
                    while v == pattern[j]:
                        v = alphabet[random.randint(0, 3)]
                    result.append(v)
                    # print('  Mismatch at %d: %s -> %s' % (j, pattern[j], v))
                elif var_type == 1:  # Deletion
                    pass
                    # print('  Deletion at %d' % j)
                else:  # Insertion
                    v = alphabet[random.randint(0, 3)]
                    # print('  Insertion at %d: %s' % (j, v))
                    result.append(v)
                    result.append(pattern[j])
        result = ''.join(result)
        # print('  Mutated motif: %s\n' % result)
        text1 += result

    mutated_n2 = 0
    for i in range(r2):
        result = []
        for j in range(motif_len):
            # If the picked number from the pool is 1, then a mutation happens
            if mut_pool == -1 or random.randint(1, mut_pool) != 1:
                result.append(pattern[j])
            else:
                mutated_n2 += 1
                var_type = random.randint(0, 2)
                if var_type == 0:  # Mismatch
                    v = alphabet[random.randint(0, 3)]
                    while v == pattern[j]:
                        v = alphabet[random.randint(0, 3)]
                    result.append(v)
                    # print('  Mismatch at %d: %s -> %s' % (j, pattern[j], v))
                elif var_type == 1:  # Deletion
                    pass
                    # print('  Deletion at %d' % j)
                else:  # Insertion
                    v = alphabet[random.randint(0, 3)]
                    # print('  Insertion at %d: %s' % (j, v))
                    result.append(v)
                    result.append(pattern[j])
        result = ''.join(result)
        # print('  Mutated motif: %s\n' % result)
        text2 += result

    for i in range(0, flank_r):
        c = alphabet[random.randint(0, 3)]
        text1 += c
        if var_pool == -1 or random.randint(1, var_pool) != 1:
            text2 += c
        else:
            var_type = random.randint(0, 2)
            if var_type == 0:  # Mismatch
                v = alphabet[random.randint(0, 3)]
                while v == c:
                    v = alphabet[random.randint(0, 3)]
                text2 += v
            elif var_type == 1:  # Deletion
                pass
            else:  # Insertion
                v = alphabet[random.randint(0, 3)]
                text2 += c
                text2 += v

    tab['ID'].append(len(tab['ID']) + 1)
    tab['motif'].append(pattern)
    tab['periods'].append(r1)
    tab['mutation'].append(mutated_n1)
    tab['flank_l'].append(flank_l)
    tab['flank_r'].append(flank_r)
    tab['text'].append(text1)

    tab['ID'].append(len(tab['ID']) + 1)
    tab['motif'].append(pattern)
    tab['periods'].append(r2)
    tab['mutation'].append(mutated_n1)
    tab['flank_l'].append(flank_l)
    tab['flank_r'].append(flank_r)
    tab['text'].append(text2)


def combine_text(text_array: tuple, tab: pd):
    tab['ID'].append(len(tab['ID']) + 1)
    tab['motif'].append('NA')
    tab['periods'].append('0')
    tab['mutation'].append('0')
    tab['flank_l'].append('0')
    tab['flank_r'].append('0')
    tab['text'].append(''.join(text_array))


if __name__ == '__main__':
    table = {'ID': [], 'motif': [], 'periods': [], 'mutation': [], 'flank_l': [], 'flank_r': [], 'text': []}
    text1 = gen_repeat(motif_len=3, repeat_n=15, mutation_rate=0.0, flank_l=8, flank_r=8, tab=table)
    text2 = gen_repeat(motif_len=5, repeat_n=10, mutation_rate=0.0, flank_l=20, flank_r=20, tab=table)
    combine_text(text_array=(text1, text2), tab=table)
    gen_vntr(motif_len=5, r1=10, r2=15, var_rate=0.02, mut_rate=0.0, flank_l=20, flank_r=25, tab=table)
    gen_repeat(motif_len=10, repeat_n=20, mutation_rate=0.03, flank_l=80, flank_r=80, tab=table)
    gen_repeat(motif_len=50, repeat_n=40, mutation_rate=0.05, flank_l=150, flank_r=125, tab=table)
    gen_vntr(motif_len=20, r1=7, r2=11, var_rate=0.01, mut_rate=0.01, flank_l=120, flank_r=90, tab=table)
    gen_repeat(motif_len=100, repeat_n=100, mutation_rate=0.01, flank_l=1000, flank_r=1000, tab=table)
    df = pd.DataFrame(table)
    df.to_csv('motif.csv', index=False)
