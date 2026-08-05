# PyMolGen Drug-Like Molecule Generation

Generates drug-like molecules from a parent compound by deriving fragment combination rules from a molecular database (ChEMBL). Starting from a core structure, new analogues are randomly sampled according to the source databases fragment-bond frequencies, deduplicated, and checked for RDKit validity. Each output is a distinct generated analogue of the input compound.

This model was incorporated on 2026-08-04.Last packaged on 2026-08-04.

## Information
### Identifiers
- **Ersilia Identifier:** `eos8zvb`
- **Slug:** `pymolgen`

### Domain
- **Task:** `Sampling`
- **Subtask:** `Generation`
- **Biomedical Area:** `Any`
- **Target Organism:** `Any`
- **Tags:** `Compound generation`, `Drug-likeness`, `ChEMBL`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `100`
- **Output Consistency:** `Variable`
- **Interpretation:** Generated drug-like analogues of the input parent compound sampled from ChEMBL fragment combination rules.

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| smi_00 | string |  | Generated drug-like analogue index 0 sampled from ChEMBL fragment combination rules |
| smi_01 | string |  | Generated drug-like analogue index 1 sampled from ChEMBL fragment combination rules |
| smi_02 | string |  | Generated drug-like analogue index 2 sampled from ChEMBL fragment combination rules |
| smi_03 | string |  | Generated drug-like analogue index 3 sampled from ChEMBL fragment combination rules |
| smi_04 | string |  | Generated drug-like analogue index 4 sampled from ChEMBL fragment combination rules |
| smi_05 | string |  | Generated drug-like analogue index 5 sampled from ChEMBL fragment combination rules |
| smi_06 | string |  | Generated drug-like analogue index 6 sampled from ChEMBL fragment combination rules |
| smi_07 | string |  | Generated drug-like analogue index 7 sampled from ChEMBL fragment combination rules |
| smi_08 | string |  | Generated drug-like analogue index 8 sampled from ChEMBL fragment combination rules |
| smi_09 | string |  | Generated drug-like analogue index 9 sampled from ChEMBL fragment combination rules |

_10 of 100 columns are shown_
### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos8zvb](https://hub.docker.com/r/ersiliaos/eos8zvb)
- **Docker Architecture:** `AMD64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos8zvb.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos8zvb.zip)

### Resource Consumption
- **Model Size (Mb):** `42`
- **Environment Size (Mb):** `1030`
- **Image Size (Mb):** `1143.06`

**Computational Performance (seconds):**
- 10 inputs: `38.84`
- 100 inputs: `46.94`
- 10000 inputs: `1448.8`

### References
- **Source Code**: [https://github.com/HirstGroup/PyMolGen](https://github.com/HirstGroup/PyMolGen)
- **Publication**: [https://doi.org/10.1021/acs.jcim.6c00689](https://doi.org/10.1021/acs.jcim.6c00689)
- **Publication Type:** `Peer reviewed`
- **Publication Year:** `2026`
- **Ersilia Contributor:** [arnaucoma24](https://github.com/arnaucoma24)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [MIT](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos8zvb
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos8zvb
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
