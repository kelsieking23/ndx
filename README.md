# ndx
Creates custom index files for GROMACS


**ndx.py** creates custom index files for GROMACS. The script makes use of selection strings that are intuitive, making complex selections easier.
Some notes on selection:
* Can select many GROMACS default index groups: backbone, mainchain, sidechain, protein, nonprotein, solvent, lipids, ions
* Can select by atom type (CA, C, O, N, H, P, S), or by specific atom name
* Can select by residue ID <residue_number><residue_name> (ex:110CYS)
* Can select by chain. Chains are labelled p1, p2, p3... etc
* Can select by residue index. Residue indeces are labelled ri1, ri2, ri3...
* Currently, only logical operators supported are parentheses, "and" and "or" 

## Usage
Positional arguments:

**gro** _(.gro)_ Structure file

**-selections** _(str)_ Selection string(s). Each individual selection must be enclosed in quotations. Run with option --show-examples to see demonstrations.

**-ligands** _(str)_ Residue names for ligand(s) or other non-standard residues in the structure file, if any. If not specified and non-residue ligands are present, will be grouped under 'non protein'

**-o, --output** _(default = index.ndx)_ Output index file

**--show-examples** _(default = False)_ Show examples and quit 

**-h, --help** Show help message and quit

## Examples
* Selecting by chain (Chain A and B in a single index): `python ndx.py structure.gro -selections "p1 or p2" -o index.ndx`
* Selecting by chain (Chain A and B in different indices): `python ndx.py structure.gro -selections "p1" "p2" -o index.ndx`
* Selecting by residue index: `python ndx.py structure.gro -selections "ri1-ri15" "ri16-ri30" -o index.ndx`
* Selecting by residue index, only backbone `python ndx.py structure.gro -selections "(ri1-ri15) and backbone" "(ri16-ri30) and backbone"`
