from pathlib import Path
import psi4
import numpy as np

psi4.set_num_threads(16)
psi4.core.clean_options()
psi4.set_memory('32 GB')

path = Path('./xyz')
for item in path.iterdir():
    if item.is_file() and str(item)[-4:] == '.xyz':
        name = str(item)
        print(name, end='\t')
        ss = open(name, 'r').read()
        ss = '\n'.join(ss.split('\n')[2:])

        psi4.set_output_file(f'{name}.out', False)
        mol = psi4.geometry(ss)

        # Set computational options
        psi4.set_options({
            'scf_type': 'pk',
            'dft_pruning_scheme': 'none',
            'e_convergence': 1e-8,
            'd_convergence': 1e-8,
            'basis': 'cc-pvdz',
            'reference': 'rks', 
            'dft_radial_scheme': 'becke', 
            'dft_nuclear_scheme': 'becke'
        })

        energy = psi4.energy('blyp') # pbe, etc.
        print(energy)
