def calculate_immonium_ions(peptide):
    """Calculate immonium ions for a peptide sequence."""
    # Get unique amino acids in the peptide
    unique_aas = set(peptide)
    
    # Calculate immonium ions for each unique amino acid
    immonium_ions = {}
    for aa in unique_aas:
        if aa in immonium_masses:
            immonium_ions[f'im_{aa}'] = round(immonium_masses[aa], 4)
    
    return dict(sorted(immonium_ions.items(), key=lambda x: x[1]))

def calculate_internal_fragments(peptide, charge=1):
    """Calculate internal fragment ions for a peptide sequence."""
    proton = 1.007276
    water = 18.010565
    
    internal_ions = {}
    
    # Generate all possible internal fragments of length 2 or more
    for i in range(len(peptide)-1):
        for j in range(i+2, len(peptide)+1):
            fragment = peptide[i:j]
            if len(fragment) >= 2:  # Only consider fragments of length 2 or more
                # Calculate mass of internal fragment
                mass = sum(amino_acid_masses[aa] for aa in fragment)
                # Add proton for charge
                mz = (mass + proton) / charge
                internal_ions[f'int_{i+1}_{j+1}'] = round(mz, 4)
    
    return dict(sorted(internal_ions.items(), key=lambda x: x[1]))

def calculate_peptide_ions(peptide, charge=1):
    """
    Calculate theoretical m/z values for b, y, c, z ions, immonium ions, and internal fragments for a peptide.
    """
    # Get the existing ion calculations
    ions = {
        'b': [],
        'y': [],
        'c': [],
        'z': []
    }
    
    # Constants for ion calculations
    proton = 1.007276
    water = 18.010565
    NH3 = 17.0265

    # Calculate b ions (N-terminal fragments)
    cumulative = 0.0
    for i in range(len(peptide) - 1):
        cumulative += amino_acid_masses[peptide[i]]
        ions['b'].append(round((cumulative + proton) / charge, 4))
    ions['b'].sort()

    # Calculate y ions (C-terminal fragments)
    cumulative = 0.0
    for i in range(len(peptide) - 1, 0, -1):
        cumulative += amino_acid_masses[peptide[i]]
        ions['y'].append(round((cumulative + water + proton) / charge, 4))
    ions['y'].sort()

    # Calculate c ions (N-terminal fragments + NH3)
    cumulative = 0.0
    for i in range(len(peptide) - 1):
        cumulative += amino_acid_masses[peptide[i]]
        ions['c'].append(round((cumulative + NH3 + proton) / charge, 4))
    ions['c'].sort()

    # Calculate z ions (C-terminal fragments - NH3)
    cumulative = 0.0
    for i in range(len(peptide) - 1, 0, -1):
        cumulative += amino_acid_masses[peptide[i]]
        ions['z'].append(round((cumulative + proton - NH3) / charge, 4))
    ions['z'].sort()

    # Add immonium ions
    ions['immonium'] = calculate_immonium_ions(peptide)
    
    # Add internal fragments
    ions['internal'] = calculate_internal_fragments(peptide, charge)

    return ions

def calculate_n_glycopeptide_ions(peptide, glycan_composition, glycan_frag_order=None, charge=1):
    """
    Calculate theoretical m/z values for b, y, c, z, Y, B, oxonium ions, immonium ions, internal ions,
    and Y ions with peptide losses for an N-glycopeptide.
    All ion series are sorted from smallest to largest m/z.

    Parameters:
      peptide (str): The peptide sequence (e.g., "NTSK").
      glycan_composition (str): A string like "HexNAc(5)Hex(5)dHex(1)NeuAc(2)".
      glycan_frag_order (list, optional): A list specifying the sugar loss order.
      charge (int): The charge state (default is 1).

    Returns:
      dict: A dictionary with keys for ion types and their m/z values.
    """
    # Get peptide backbone ions first
    peptide_ions = calculate_peptide_ions(peptide, charge)
    
    # Constants
    proton = 1.007276
    water = 18.010565

    # Find the glycosite position (N-X-S/T motif)
    glycosite_match = re.search(r'N[^P][ST]', peptide)
    if not glycosite_match:
        return peptide_ions  # Return only peptide ions if no glycosite found
    
    glycosite_start = glycosite_match.start()
    glycosite_end = glycosite_match.end()

    # Compute the neutral mass of the peptide (including water)
    peptide_mass = sum(amino_acid_masses[aa] for aa in peptide) + water

    # Parse glycan composition into a dictionary
    glycan_dict = {}
    for part in glycan_composition.split(')'):
        if part:
            sugar, count = part.split('(')
            glycan_dict[sugar] = int(count)

    # --- Calculate immonium ions ---
    immonium_ions = {}
    for aa in peptide:
        if aa in immonium_masses:
            immonium_ions[f'im_{aa}'] = round(immonium_masses[aa] + proton, 4)

    # --- Calculate internal ions ---
    internal_ions = {}
    for i in range(len(peptide) - 2):
        for j in range(i + 2, len(peptide)):
            # Skip if internal fragment would break the glycosite
            if i < glycosite_start and j > glycosite_end:
                continue
            internal_seq = peptide[i:j+1]
            internal_mass = sum(amino_acid_masses[aa] for aa in internal_seq)
            internal_ions[f'internal_{i+1}_{j+1}'] = round((internal_mass + proton) / charge, 4)

    # --- Calculate Y ions (glycan-attached peptide fragments) ---
    Y_ions = {}
    current_mass = peptide_mass
    Y_ions['Y0'] = round((current_mass + proton) / charge, 4)

    # Process glycan additions in order
    remaining_sugars = []
    for sugar, count in glycan_dict.items():
        remaining_sugars.extend([sugar] * count)

    # Sort sugars by mass for consistent ordering
    remaining_sugars.sort(key=lambda x: monosaccharide_library[x]['mass'])

    # Add sugars in order
    for i, sugar in enumerate(remaining_sugars, start=1):
        current_mass += monosaccharide_library[sugar]['mass']
        Y_ions[f'Y{i}'] = round((current_mass + proton) / charge, 4)

    # --- Calculate Y ions with peptide backbone losses ---
    Y_peptide_loss = {}
    
    # Calculate b-ion masses (without charge) - only for positions that retain glycosite
    b_masses = []
    cumulative = 0.0
    for i in range(len(peptide) - 1):
        cumulative += amino_acid_masses[peptide[i]]
        # Only include b-ions that would leave the glycosite intact
        if i < glycosite_start:
            b_masses.append((i+1, cumulative))

    # Calculate y-ion masses (without charge) - only for positions that retain glycosite
    y_masses = []
    cumulative = 0.0
    for i in range(len(peptide) - 1, 0, -1):
        cumulative += amino_acid_masses[peptide[i]]
        # Only include y-ions that would leave the glycosite intact
        if i > glycosite_end:
            y_masses.append((len(peptide) - i, cumulative + water))

    # For each Y ion, subtract each valid b-ion and y-ion mass
    for y_num, y_mass in Y_ions.items():
        y_neutral = (y_mass * charge) - proton
        
        # Subtract b-ions (N-terminal losses)
        for b_idx, b_mass in b_masses:
            fragment_mass = y_neutral - b_mass
            Y_peptide_loss[f'{y_num}-b{b_idx}'] = round((fragment_mass + proton) / charge, 4)
        
        # Subtract y-ions (C-terminal losses)
        for y_idx, y_mass in y_masses:
            fragment_mass = y_neutral - y_mass
            Y_peptide_loss[f'{y_num}-y{y_idx}'] = round((fragment_mass + proton) / charge, 4)

    # --- Calculate B ions (glycan fragment ions) ---
    B_ions = {}
    cumulative = 0.0
    
    # Sort sugars by mass
    sorted_sugars = sorted(glycan_dict.items(), key=lambda x: monosaccharide_library[x[0]]['mass'])
    
    # Calculate B ions in order
    for i, (sugar, count) in enumerate(sorted_sugars, start=1):
        for j in range(count):
            cumulative += monosaccharide_library[sugar]['mass']
            B_ions[f'B{i}_{j+1}'] = round((cumulative + proton) / charge, 4)

    # --- Calculate Oxonium ions ---
    oxonium_ions = {}
    for sugar in sorted(glycan_dict.keys(), key=lambda x: monosaccharide_library[x]['mass']):
        if glycan_dict[sugar] > 0:
            oxonium_ions[f'ox_{sugar}'] = round(monosaccharide_library[sugar]['mass'] + proton, 4)

    return {
        'b': peptide_ions['b'],
        'y': peptide_ions['y'],
        'c': peptide_ions['c'],
        'z': peptide_ions['z'],
        'immonium': dict(sorted(immonium_ions.items(), key=lambda item: item[1])),
        'internal': dict(sorted(internal_ions.items(), key=lambda item: item[1])),
        'Y': dict(sorted(Y_ions.items(), key=lambda item: item[1])),
        'Y_peptide_loss': dict(sorted(Y_peptide_loss.items(), key=lambda item: item[1])),
        'B': dict(sorted(B_ions.items(), key=lambda item: item[1])),
        'oxonium': dict(sorted(oxonium_ions.items(), key=lambda item: item[1]))
    }

## Experimental -- testing
#from itertools import permutations
def generate_all_y_ions(peptide, glycan_composition, charge=1):
    """
    Generate all possible Y ion series considering different fragmentation paths.

    Parameters:
      peptide (str): Peptide sequence (e.g., "NTSK").
      glycan_composition (str): Glycan composition string (e.g., "HexNAc(5)Hex(5)dHex(1)NeuAc(2)").
      charge (int): The charge state (default is 1).

    Returns:
      dict: A dictionary where keys are fragmentation paths and values are Y ion series.
    """
    # Constants
    proton = 1.007276
    water  = 18.010565

    # Compute the neutral mass of the peptide
    peptide_mass = sum(amino_acid_masses[aa] for aa in peptide) + water

    # Parse glycan composition
    glycan_dict = {}
    for part in glycan_composition.split(')'):
        if part:
            sugar, count = part.split('(')
            glycan_dict[sugar] = int(count)

    # Convert glycan composition into a list of sugars
    glycan_list = []
    for sugar, count in glycan_dict.items():
        glycan_list.extend([sugar] * count)  # Expand each sugar into individual occurrences

    # Generate all possible fragmentation paths
    unique_permutations = set(permutations(glycan_list))  # Unique orders only

    # Compute Y ions for each fragmentation path
    all_Y_ions = {}

    for perm in unique_permutations:
        path_name = " -> ".join(perm)  # Name this path

        current_mass = peptide_mass  # Start with peptide only
        y_series = {'Y0': round((current_mass + proton) / charge, 4)}

        for i, sugar in enumerate(perm, start=1):
            current_mass += monosaccharide_library[sugar]['mass']
            y_series[f'Y{i}'] = round((current_mass + proton) / charge, 4)

        all_Y_ions[path_name] = y_series

    return all_Y_ions