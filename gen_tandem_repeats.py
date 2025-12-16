import random
import pandas as pd


def gen_repeat(motif_len: int, repeat_n: int, mutation_rate: float, tab: dict):
    # Generate a random motif with the given length
    alphabet = 'ACGT'
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
    text = ''
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
    tab['ID'].append(len(tab['ID']) + 1)
    tab['motif'].append(pattern)
    tab['periods'].append(repeat_n)
    tab['mutation'].append(mutated_n)
    tab['text'].append(text)


if __name__ == '__main__':
    table = {'ID': [], 'motif': [], 'periods': [], 'mutation': [], 'text': []}
    gen_repeat(motif_len=3, repeat_n=5, mutation_rate=0.0, tab=table)
    gen_repeat(motif_len=5, repeat_n=8, mutation_rate=0.02, tab=table)
    gen_repeat(motif_len=10, repeat_n=20, mutation_rate=0.05, tab=table)
    gen_repeat(motif_len=50, repeat_n=40, mutation_rate=0.07, tab=table)
    gen_repeat(motif_len=100, repeat_n=100, mutation_rate=0.1, tab=table)
    df = pd.DataFrame(table)
    df.to_csv('motif.csv', index=False)
