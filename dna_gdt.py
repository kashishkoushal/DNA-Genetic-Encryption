"""
DNA Genetic Decryption Technique
"""
from time import time
import ast
import re  # For punctuation removal

import utils
from utils import *

def encrypt_key(data, key):
    if len(data) > len(key):
        factor = int(len(data) / len(key))
        key += key * factor
        return bitxor(data, key)
    return bitxor(data, key)

def reshape(dna_sequence, reshape_info):
    chromosome_length = int(reshape_info[0])
    chromosomes = []
    for i in range(0, len(dna_sequence), chromosome_length):
        chromosomes.append(dna_sequence[i:i + chromosome_length])
    return chromosomes

def reverse_reshape(population):
    return "".join(population)

def rotate_crossover(population, rotate_info):
    new_population = []
    rotation_offset = int(get_pattern(rotation_offset_del, rotate_info)[0])
    rotations = get_pattern(rotation_types_del, rotate_info)[0].split("|")[:-1]

    for i in range(len(population)):
        chromosome = population[i]
        direction = rotations[i]
        if direction == "left":
            part1 = chromosome[:-rotation_offset]
            part2 = chromosome[-rotation_offset:]
            new_population.append(part2 + part1)
        elif direction == "right":
            part1 = chromosome[:rotation_offset]
            part2 = chromosome[rotation_offset:]
            new_population.append(part2 + part1)
    return new_population

def single_point_crossover(population, single_point_info):
    crossover_points = [int(p) for p in single_point_info.split("|") if p]
    new_population = []
    for i in range(0, len(population) - 1, 2):
        c1 = population[i]
        c2 = population[i + 1]
        point = crossover_points[int(i / 2)]
        new_population.append(c2[:point] + c1[point:])
        new_population.append(c1[:point] + c2[point:])
    if len(population) % 2 == 1:
        new_population.append(population[-1])
    return new_population

def crossover(population, crossover_info):
    crossover_type = get_pattern(crossover_type_del, crossover_info)[0]
    if crossover_type == "rotate_crossover":
        rotate_info = get_pattern(rotate_crossover_del, crossover_info)[0]
        return rotate_crossover(population, rotate_info)
    elif crossover_type == "single_point_crossover":
        single_point_info = get_pattern(single_point_crossover_del, crossover_info)[0]
        return single_point_crossover(population, single_point_info)
    elif crossover_type == "both":
        rotate_info = get_pattern(rotate_crossover_del, crossover_info)[0]
        single_point_info = get_pattern(single_point_crossover_del, crossover_info)[0]
        population = single_point_crossover(population, single_point_info)
        return rotate_crossover(population, rotate_info)

def complement(chromosome, point1, point2):
    return ''.join(
        '1' if point1 <= i <= point2 and bit == '0'
        else '0' if point1 <= i <= point2 and bit == '1'
        else bit
        for i, bit in enumerate(chromosome)
    )

def mutation(population, mutation_info):
    alter_dna_table = ast.literal_eval(get_pattern(mutation_table_del, mutation_info[0])[0])
    chromosomes_info = get_pattern(chromosome_del, mutation_info[0])
    new_population = []

    for i in range(len(population)):
        chromosome = population[i]
        chromosome_info = chromosomes_info[i]

        alter_info = get_pattern(alter_mutation_del, chromosome_info)[0]
        point1, point2 = ast.literal_eval(alter_info)
        new_chromosome = ''
        for j in range(len(chromosome)):
            if point1 <= j <= point2:
                new_chromosome += alter_dna_table[chromosome[j]]
            else:
                new_chromosome += chromosome[j]

        two_bases_vector = group_bases(new_chromosome)
        last_two_bits = None
        if len(new_chromosome) % 2 == 1:
            last_two_bits = utils.dna_base_to_two_bits_table[new_chromosome[-1]]
            two_bases_vector = two_bases_vector[:-1]

        bits_seq = dna_to_bits(two_bases_vector, utils.two_dna_base_to_four_bits_table)
        if last_two_bits:
            bits_seq += last_two_bits

        complement_info = get_pattern(complement_mutation_del, chromosome_info)[0]
        point1, point2 = ast.literal_eval(complement_info)
        b_chromosome = complement(bits_seq, point1, point2)
        b_chromosome = group_bits(b_chromosome)
        new_chromosome = bits_to_dna(b_chromosome, utils.two_bits_to_dna_base_table)
        new_population.append(new_chromosome)

    return new_population

def dna_gdt(text, key):
    print("\n🔓 DNA-GDT Decryption Running...\n")

    rounds_no = int(get_pattern(no_rounds_del, key)[0])
    rounds = get_pattern(round_del, key)
    b_data = text

    print("Initial DNA sequence:", b_data)

    while rounds_no > 0:
        round_info = rounds[rounds_no - 1]

        b_data = reshape(b_data, get_pattern(reshape_del, round_info))
        b_data = mutation(b_data, get_pattern(mutation_del, round_info))
        b_data = crossover(b_data, round_info)

        encryption_key = get_pattern(key_del, key)[0]
        b_data = bits_to_dna(
            group_bits(
                encrypt_key(
                    dna_to_bits(reverse_reshape(b_data), utils.dna_base_to_two_bits_table),
                    encryption_key
                )
            ),
            utils.two_bits_to_dna_base_table
        )

        rounds_no -= 1

    decrypted_bin = dna_to_bits(b_data, utils.dna_base_to_two_bits_table)
    return bin2str(decrypted_bin).decode()

def main():
    original_file = open(original_filename, "r")
    encrypted_file = open(encrypted_filename, "r")
    decrypted_file = open(decrypted_filename, "w")
    key_file = open(key_filename, "r")

    original_text = original_file.read()
    encrypted_text = encrypted_file.read()
    key = key_file.read()

    print("Encrypted text:", encrypted_text)
    print("Key:", key)

    generate_pre_processing_tables()
    generate_mutation_tables()

    start = time()
    decrypted_text = dna_gdt(encrypted_text, key)

    # ✅ Remove punctuation: , . ( )
    decrypted_text_cleaned = re.sub(r'[.,()]', '', decrypted_text)

    print("\n🧾 Decrypted text :")
    print("=" * 80)
    print(decrypted_text_cleaned)
    print("=" * 80)

    decrypted_file.write(decrypted_text_cleaned)
    end = time()

    print("\n⏱ Total execution time:", end - start)

    original_clean = re.sub(r'[.,()]', '', original_text)

    if original_clean != decrypted_text_cleaned:
        print("\n❌ Decryption failed.")
    else:
        print("\n✅ Decryption succeeded.")

    original_file.close()
    encrypted_file.close()
    decrypted_file.close()
    key_file.close()

if __name__ == '__main__':
    main()
