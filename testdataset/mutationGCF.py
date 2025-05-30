import os
import random

def mutate_sequence(sequence, mutation_rate):
    bases = ['A', 'C', 'G', 'T']
    mutated_sequence = ""
    
    for base in sequence:
        if base not in bases:
            mutated_sequence += base  # Keep non-DNA characters unchanged
        elif random.random() < mutation_rate:
            possible_bases = [b for b in bases if b != base]
            new_base = random.choice(possible_bases)
            mutated_sequence += new_base
        else:
            mutated_sequence += base
    
    return mutated_sequence

def get_mutation_rate():
    while True:
        try:
            mutation_rate = float(input("Enter the mutation rate (between 0 and 1): "))
            if 0 <= mutation_rate <= 1:
                return mutation_rate
            else:
                print("Mutation rate must be between 0 and 1.")
        except ValueError:
            print("Invalid input. Please enter a number.")

def main():
    input_folder = input("Enter the path to the input folder: ")
    output_folder = input("Enter the path to the output folder: ")

    mutation_rate = get_mutation_rate()

    os.makedirs(output_folder, exist_ok=True)

    for filename in os.listdir(input_folder):
        if filename.startswith('GCF_') and filename.endswith('.fna'):
            input_path = os.path.join(input_folder, filename)
            output_path = os.path.join(output_folder, filename)
            
            with open(input_path, 'r') as file:
                lines = file.readlines()
                header = lines[0] if lines[0].startswith('>') else ""
                original_sequence = ''.join(lines[1:]) if header else ''.join(lines)
                original_sequence = original_sequence.strip()

            mutated_sequence = mutate_sequence(original_sequence, mutation_rate)

            with open(output_path, 'w') as file:
                if header:
                    file.write(header)
                file.write(mutated_sequence)

            print(f"Processed and mutated: {filename}")

    print(f"All sequences have been mutated and written to {output_folder}")

if __name__ == "__main__":
    main()
