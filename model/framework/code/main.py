import os
import sys
import numpy as np
from ersilia_pack_utils.core import read_smiles, write_out

input_file = sys.argv[1]
output_file = sys.argv[2]

root = os.path.dirname(os.path.abspath(__file__))
checkpoints = os.path.abspath(os.path.join(root, "..", "..", "checkpoints"))
sys.path.insert(0, checkpoints)
sys.path.insert(0, root)

from generate_analogues import generate_analogues

N_GENERATE = 100

_, smiles_list = read_smiles(input_file)

outputs = []
for smi in smiles_list:
    analogues = generate_analogues(smi, checkpoints, n=N_GENERATE)
    outputs.append(analogues)

outputs = np.array(outputs, dtype=object)

assert len(smiles_list) == len(outputs)

headers = [f"smi_{str(i).zfill(2)}" for i in range(N_GENERATE)]

write_out(outputs, headers, output_file, str)
