# ndx.py
# Version 1.0
# Written by Kelsie M. King
# Github: kelsieking23
# Contact for issues and requests: kelsieking23@vt.edu
# Last updated: 10/14/2022

import os
import sys
import argparse
from argparse import RawTextHelpFormatter

class Ndx:

    def __init__(self, gro, peptides=1, ndxt=None, ligands=None):
        self.gro = gro
        self.peptides= peptides
        self.ndxt_groups = ndxt
        if ndxt is not None:
            self.ndxt_groups = self.parseNDXT(ndxt)
        types = self.getTypes(self.ndxt_groups, ligands)
        self.types = {}
        for key in types.keys():
            self.types[key] = Type(key, types[key])
        self.selections = {}

    @property
    def residues(self):
        residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HSD']
        residues = [item.lower() for item in residues]
        return residues
    
    @property
    def caps(self):
        caps = ['ace', 'nh2', 'nme']
        return caps
    
    @property
    def lipids(self):
        lipids = ['POPC', 'CHL1', 'SDPE', 'POPE', 'PSM', 'SOPS', 'POPE', 'POPS', 'SM', 'CHOL', 'DLPG', 'DDPC']
        lipids = [item.lower() for item in lipids]
        return lipids
    
    @property
    def backbone(self):
        backbone = ['ca', 'c', 'n']
        return backbone
    
    @property
    def mainchain(self):
        mainchain = ['ca', 'c', 'o', 'n', 'hn', 'h', 'ha']
        return mainchain
    
    @property
    def ions(self):
        ions = ['k', 'cl', 'na', 'sod', 'cla']
        return ions
    
    @property
    def solvent(self):
        solvent = ['sol', 'tip3p']
        return solvent
    
    @property
    def ri(self):
        return 'ri'
    
    @property
    def chain(self):
        return 'p'
    
    def parseNDXT(self, ndxt):
        ndxt_groups = {}
        with open(ndxt, 'r') as f:
            for line in f:
                name = line.split(':')[0]
                atom_names = line.split(':')[1].split(',')
                atom_names[-1] = atom_names[-1].strip()
                ndxt_groups[name] = atom_names
                
        if ndxt_groups != {}:
            return ndxt_groups

    def select(self, sele):
        if isinstance(sele, str):
            selection = self.parseSelection(sele)
            self.selections[selection] = self.interpret(selection)
            return self.interpret(selection)


    def parseSelection(self, sele):
        parsed = []
        sele_in = sele
        sele_split = sele.split()
        # parse hyphens
        i = 0
        for exp in sele_split:
            if '-' in exp:
                items = exp.split('-')
                items = [item.strip('(').strip(')').strip() for item in items]
                # print(items[0][-3:], print(type(items[0])))
                if items[0] in self.types.keys():
                    if items[0].startswith('ri'):
                        index_1 = int(items[0][2:])
                        index_2 = int(items[-1][2:])+1
                        residue_ids = ['ri{}'.format(index) for index in range(index_1, index_2)]
                        sele_string = ' or '.join(residue_ids)
                        sele_string = '({})'.format(sele_string)
                        sele_split[i] = sele_string
                    if (items[0][-3:] in self.residues):
                        rimap = self.mapRI()
                        ri1 = rimap[items[0]]
                        ri2 = rimap[items[-1]]
                        index_1 = int(ri1[2:])
                        index_2 = int(ri2[2:])+1
                        residue_ids = ['ri{}'.format(index) for index in range(index_1, index_2)]
                        sele_string = ' or '.join(residue_ids)
                        sele_string = '({})'.format(sele_string)
                        sele_split[i] = sele_string
                    if items[0].startswith('p'):
                        index_1 = int(items[0][1:])
                        index_2 = int(items[-1][1:])+1
                        residue_ids = ['p{}'.format(index) for index in range(index_1, index_2)]
                        sele_string = ' or '.join(residue_ids)
                        sele_string = '({})'.format(sele_string)
                        sele_split[i] = sele_string
            i += 1
        sele = ' '.join(sele_split)
        # parse parentheses
        i = 0
        if '(' in sele:
            sele_string = ''
            subgroups = self.parse_parentheses(sele)
            for grp in subgroups:
                if isinstance(grp, list):
                    _sele_string = ''.join(grp)
                    atoms = self.interpret(_sele_string)
                    type_string = sele_in.split()[i]
                    self.types[type_string] = Type(type_string, atoms)
                    sele_string += type_string
                    # parsed.append(_sele_string)

                else:
                    sele_string += grp
            parsed.append(sele_string)
        i += 1
        if parsed == []:
            return sele
        return parsed[0]
                
    def push(self, obj, l, depth):
        while depth:
            l = l[-1]
            depth -= 1

        l.append(obj)

    def parse_parentheses(self, s):
        groups = []
        depth = 0

        try:
            for char in s:
                if char == '(':
                    self.push([], groups, depth)
                    depth += 1
                elif char == ')':
                    depth -= 1
                else:
                    self.push(char, groups, depth)
        except IndexError:
            raise ValueError('Parentheses mismatch')

        if depth > 0:
            raise ValueError('Parentheses mismatch')
        else:
            return groups
    
    def interpret(self, sele):
        sele_split = sele.split()

        atoms = []
        # ## and comes first but we working on or rn
        # i = 0
        sele_split_or = [item.strip('(').strip(')').strip() for item in sele.split('or')]
        sele_split_or = [item.strip('(').strip(')').strip() for item in sele_split_or]
        if len(sele_split_or) > 1:
            for item in sele_split_or:
                item
                atoms += self.getAtoms(item)

        sele_split_and = [item.strip('(').strip(')').strip() for item in sele.split('and')]
        sele_split_and = [item.strip('(').strip(')').strip() for item in sele_split_and]
        if len(sele_split_and) > 1:
            sets = []
            for item in sele_split_and:
                sets.append(set(self.getAtoms(item)))
            if len(sets) > 2:
                atoms = list(sets[0].intersection(*sets[1:]))
            else:
                atoms = list(sets[0].intersection(sets[1]))
        if atoms == []:
            return self.getAtoms(sele)
        return atoms
    
    def getAtoms(self, sele):
        return self.types[sele].items

    def mapRI(self):
        i = 0
        ri = 0
        last_id = None
        rimap = {}
        for line in self.lineGenerator():
            if i < 2:
                i += 1
                continue
            res_id = line[0:10].strip()
            res_name = line[5:10].strip()
            if (res_name in self.residues):
                if res_id != last_id:
                    ri += 1
                rimap[res_id] = 'ri{}'.format(ri)
            last_id = res_id
            i += 1
        return rimap
                
    def lineGenerator(self):
        '''
        For memory saving, generates lines for any given file
        '''
        f = open(self.gro, 'r')
        contents = f.readlines()
        f.close()
        for line in contents:
            yield line
    
    def getTypes(self, ndxt_groups=None, ligands=None):
        '''
        Assigns atoms as/by:
        ** backbone
        ** mainchain,  mainchain_h, mainchain_nocaps
        ** sidechain, sidechain_h
        ** caps
        ** nonprotein (nonprotein really specifies small molecules, lumped together into one index)
        ** ions
        ** solvent
        ** residue index
        ** residue id
        ** residue name
        ** small molecule
        Returns:
        ** dict {str:[str]}
        '''
        residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HSD']
        residues = [item.lower() for item in residues]
        caps = ['ace', 'nh2']
        lipids = ['POPC', 'CHL1', 'SDPE', 'POPE', 'PSM', 'SOPS', 'POPE', 'POPS', 'SM', 'CHOL', 'DLPG', 'DDPC']
        lipids = [item.lower() for item in lipids]
        backbone = ['ca', 'c', 'n']
        mainchain = ['ca', 'c', 'o', 'n', 'hn', 'h', 'ha']
        ions = ['k', 'cl', 'na', 'sod', 'cla']
        solvent = ['sol', 'tip3p']
        headgroup_ids = ['o7', 'p8', 'p9', 'o10', 'o11']
        types = {
            'system':[],
            'protein':[],
            'protein_caps':[],
            'backbone':[],
            'backbone_nocaps':[],
            'mainchain':[],
            'mainchain_nocaps':[],
            'mainchain_h':[],
            'mainchain_h_nocaps':[],
            'sidechain':[],
            'sidechain_h':[],
            'caps':[],
            'nonprotein':[],
            'ions':[],
            'solvent':[],
            'lipids':[],
            'drude_particles':[],
            'O':[],
            'N':[],
            'H':[],
            'C':[],
            'P':[],
            'headgroups':[],
            'headgroups_noh':[]
        }
        i = 0
        k = 0
        if ligands is None:
            ligands = []
        else:
            if isinstance(ligands, str):
                ligands = [ligands]
        # if self.system.ligands is not None:
        #     if isinstance(self.system.ligands, str):
        #         types[self.system.ligands] = []
        #         ligands.append(self.system.ligands)
        #     if isinstance(self.system.ligands, list):
        #         for l in self.system.ligands:
        #             types[l] = []
        #             ligands.append(l)
        if ndxt_groups is not None:
            to_del = []
            for key in ndxt_groups.keys():
                if key not in types.keys():
                    types[key] = []
                else:
                    to_del.append(key)
            # if to_del != []:
            #     for item in to_del:
            #         del ndxt_groups[item]
        restarted = False
        atom_index = 100000
        last_res_id = None
        visited_residues = []
        peptides = 0
        last_atom_name = None
        last_res_name = None
        for line in self.lineGenerator():
            if i < 2:
                i += 1
                continue
            else:
                line_parts = line.split()
                if len(line_parts) <= 3:
                    break
                
                res_id = line[0:10].strip().lower()
                atom_name = line[8:15].strip().lower()
                if len(atom_name.split()) == 2:
                    atom_name = atom_name.split()[-1].strip()
                res_name = line[5:10].strip().lower()
                res_num = line[:5].strip()
                atom_num = line[15:20].strip()
                if atom_num == '0':
                    atom_num = str(atom_index)
                    restarted = True
                    atom_index += 1
                elif restarted is True:
                    atom_num = str(atom_index)
                    atom_index += 1
                else:
                    pass
                types['system'].append(atom_num)
                if res_name in residues:
                    # protein atoms
                    types['protein'].append(atom_num)
                    types['protein_caps'].append(atom_num)

                    # # backbone
                    # if atom_name in backbone:
                    #     types['backbone'].append(atom_num)
                    #     types['backbone_nocaps'].append(atom_num)

                    # # mainchain
                    # if atom_name in mainchain:
                    #     types['mainchain'].append(atom_num)
                    #     types['mainchain_nocaps'].append(atom_num)
                    #     types['mainchain_h'].append(atom_num)
                    #     types['mainchain_h_nocaps'].append(atom_num)

                    # backbone
                    if atom_name in backbone:
                        types['backbone'].append(atom_num)
                        types['backbone_nocaps'].append(atom_num)
                    # mainchain
                    elif atom_name in mainchain:
                        types['mainchain'].append(atom_num)
                        types['mainchain_nocaps'].append(atom_num)
                        types['mainchain_h'].append(atom_num)
                        types['mainchain_h_nocaps'].append(atom_num)
                    else:
                        types['sidechain'].append(atom_num)

                    # atom types
                    if 'O' in atom_name:
                        types['O'].append(atom_num)
                    if 'N' in atom_name:
                        types['N'].append(atom_num)
                    if 'H' in atom_name:
                        types['H'].append(atom_num)
                    if 'C' in atom_name:
                        types['C'].append(atom_num)
                    if 'P' in atom_name:
                        types['P'].append(atom_num)
                    if atom_name not in types.keys():
                        types[atom_name] = [atom_num]
                    else:
                        types[atom_name].append(atom_num)
                    # res id 
                    if res_id not in list(types.keys()):
                        types[res_id] = []
                        types[res_id].append(atom_num)
                    else:
                        types[res_id].append(atom_num)

                    # residue index
                    if last_res_id != res_id:
                        k += 1
                        ri = 'ri{}'.format(str(k))
                        types[ri] = []
                        types[ri].append(atom_num)
                    else:
                        types[ri].append(atom_num)

                    # res name
                    if res_name not in list(types.keys()):
                        types[res_name] = []
                        types[res_name].append(atom_num)
                    else:
                        types[res_name].append(atom_num)

                    # peptide
                    if self.peptides is not None:
                        if last_res_id is None:
                                peptides += 1
                                pi = 'p{}'.format(str(peptides))
                                types[pi] = []
                                types[pi].append(atom_num)
                        else:
                            if last_res_id != res_id:
                                if int(res_num) != int(last_res_num) + 1:
                                    peptides += 1
                                    pi = 'p{}'.format(str(peptides))
                                    types[pi] = []
                                    types[pi].append(atom_num)
                                else:
                                    types[pi].append(atom_num)
                            else:
                                types[pi].append(atom_num)
                            ### This commended out block works for peptides with repeating residue numbers.
                            ### for example the abeta systems 
                            # if res_id == visited_residues[0]:
                            #     if res_id != last_res_id:
                            #         peptides += 1
                            #         pi = 'p{}'.format(str(peptides))
                            #         types[pi] = []
                            #         types[pi].append(atom_num)
                            #     else:
                            #         types[pi].append(atom_num)
                            # else:
                            #     types[pi].append(atom_num)
                    if res_id not in visited_residues:
                        visited_residues.append(res_id)
                    last_res_id = res_id
                    last_res_num = res_num
                # caps
                elif res_name in caps:
                    types['protein_caps'].append(atom_num)
                    types['caps'].append(atom_num)
                    if atom_name in backbone:
                        types['backbone'].append(atom_num)
                    if atom_name in mainchain:
                        types['mainchain'].append(atom_num)
                    if 'H' in atom_name:
                        types['mainchain_h'].append(atom_num)
                    if res_id not in types.keys():
                        types[res_id] = []
                        types[res_id].append(atom_num)
                    else:
                        types[res_id].append(atom_num)
                # lipids
                elif res_name in lipids:
                    types['lipids'].append(atom_num)
                    if res_name not in types.keys():
                        types[res_name] = [atom_num]
                    else:
                        types[res_name].append(atom_num)
                # ions
                elif res_name in ions:
                    types['ions'].append(atom_num)
                # solvent
                elif res_name in solvent:
                    types['solvent'].append(atom_num)
                    # break
                # ligands
                elif res_name in ligands:
                    types[res_name].append(atom_num)
                    if res_id not in list(types.keys()):
                        types[res_id] = []
                        types[res_id].append(atom_num)
                    else:
                        types[res_id].append(atom_num)
                # other non protein
                else:
                    types['nonprotein'].append(atom_num)
                    if res_name not in list(types.keys()):
                        types[res_name] = []
                        types[res_name].append(atom_num)
                    else:
                        types[res_name].append(atom_num)
                # ndxt groups
                if ndxt_groups is not None:
                    for key in ndxt_groups.keys():
                        if atom_name in ndxt_groups[key]:
                            if atom_num not in types[key]:
                                types[key].append(atom_num)
                        if atom_num in ndxt_groups[key]:
                            if atom_num not in types[key]:
                                types[key].append(atom_num)
                # drude particles
                if (atom_name.startswith('D')) or (atom_name.startswith('LP')):
                    types['drude_particles'].append(atom_num)

                # lipid headgroups:
                if (atom_name.startswith('P')) and (res_name in lipids):
                    types['headgroups'].append(atom_num)
                    types['headgroups'].append(str(int(atom_num) - 1))
                    types['headgroups'].append(str(int(atom_num) + 1))
                    types['headgroups'].append(str(int(atom_num) + 2))
                    types['headgroups'].append(str(int(atom_num) + 3))
                    types['headgroups_noh'].append(atom_num)
                    types['headgroups_noh'].append(str(int(atom_num) - 1))
                    types['headgroups_noh'].append(str(int(atom_num) + 1))
                    types['headgroups_noh'].append(str(int(atom_num) + 2))
                    types['headgroups_noh'].append(str(int(atom_num) + 3))
                if (res_name == 'CHOL'):
                    if ('O' in atom_name):
                        types['headgroups'].append(atom_num)
                        types['headgroups'].append(str(int(atom_num) + 1))
                        types['headgroups_noh'].append(atom_num)


        # if set_types is True:
        #     self.types = types 
        self.types = types      
        return types
    
    def makeNDX(self, filename):
        lines = self.writeLines(self.selections)
        f = open(filename, 'w')
        for line in lines:
            f.write(line)
        f.close()
        return lines

    def chunkGenerator(self, lis, n):
        '''
        Internal function to split lists into groups of 15
        '''
        for i in range(0, len(lis), n):
            try:
                yield lis[i:i+n]
            except:
                yield lis[i:]
    
    def writeLines(self, indeces):
        '''
        Internal function to generate lines for NDX file.
        '''
        strings = []
        for key in indeces.keys():
            ints = [int(num) for num in indeces[key]]
            sorted_ints = sorted(ints)
            indeces[key] = [str(i) for i in sorted_ints]
            value = indeces[key]
            if len(value) > 0:
                string = '[ {} ]\n'.format(key)
                strings.append(string)
                string = ''
                for line in self.chunkGenerator(value, 15):
                    length = 0
                    for item in line:
                        if len(str(item)) > length:
                            length = len(str(item))
                    if length < 4:
                        string = ''.join(['%4s' % i for i in line])
                    else:
                        string = ''.join(['%4s ' % i for i in line])
                    string = string + '\n'
                    strings.append(string)
        return strings


class Type:

    def __init__(self, name, items):
        self.name = name
        self.items = items
        if len(self.items) > 1:
            self.index = 0
            self.end = len(self.items) - 1
        else:
            self.index = None
            self.end = None

    def __len__(self):
        return len(self.items)

    def __add__(self, other):
        return self.items + other.items

    def __iter__(self):
        return self
    
    def __next__(self):
        if (self.index is None):
            raise StopIteration
        elif (self.index > self.end):
            raise StopIteration
        else:
            self.index += 1
            return self.items[(self.index - 1)]

def fileChoicesInput(filename, parser):
    choices = ['.gro', '.pdb', '.pdbqt']
    ext = os.path.splitext(filename)[-1]
    if ext not in choices:
        parser.print_help()
        raise argparse.ArgumentTypeError('\nStructure file must be a .gro, .pdb, or .pdbqt ({} was given)'.format(ext))

def fileChoicesOutput(filename, parser):
    choices = ['.ndx']
    if filename is None:
        return 'index.ndx'
    ext = os.path.splitext(filename)[-1]
    if ext not in choices:
        parser.print_help()
        raise argparse.ArgumentTypeError('\nOutput file must be .ndx ({} was given)'.format(ext))

def commandline():
    description = '''
             ndx.py v1.0 
Creates custom index files for GROMACS
**************************************
Usage:'''
    epilog = '''             ndx.py v1.0 
Creates custom index files for GROMACS
**************************************
Notes on selection:
* Can select many GROMACS default index groups: backbone, mainchain, sidechain, protein, nonprotein, solvent, lipids, ions
* Can select by residue ID <residue_number><residue_name> (ex:110CYS)
* Can select by chain. Chains are labelled p1, p2, p3... etc
* Can select by residue index. Residue indeces are labelled ri1, ri2, ri3...
* Currently, only logical operators supported are parentheses, "and" and "or" 
Selection Examples:
* Selecting by chain (Chain A and B in a single index): python ndx.py structure.gro -selections "p1 or p2"' -o index.ndx
* Selecting by chain (Chain A and B in different indices): python ndx.py structure.gro -selections "p1" "p2" -o index.ndx
* Selecting by residue index: python ndx.py structure.gro -selections "ri1-ri15" "ri16-ri30" -o index.ndx
* Selecting by residue index, only backbone python ndx.py structure.gro -selections "(ri1-ri15) and backbone" "(ri16-ri30) and backbone"'''
    parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

    
    parser.add_argument('gro', type=str, help='[.gro/ .pdb/ .pdbqt] (required)\n  Structure file')
    parser.add_argument('-selections', nargs='+', type=str, help='  Selection string(s). Each individual selection must be enclosed in quotations. Run with option --show-examples to see demonstrations.')
    parser.add_argument('-ligands', nargs='*', type=str, help="  Residue names for ligand(s) or other non-standard residues in the structure file, if any. If not specified and non-residue ligands are present, will be grouped under 'non protein'")
    parser.add_argument('-o', '--output', nargs='?', default='index.ndx', help='[.ndx] (default: index.ndx)\n  Output index file')
    parser.add_argument('--show-examples', action='store_true', help='Show examples and quit')

    if '--show-examples' in sys.argv:
        print(epilog)
        sys.exit(0)
    # show help if ran with no arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    errorHandler(parser, args)
    return args

def errorHandler(parser, args):
    if not os.path.exists(args.gro):
        parser.print_help()
        raise FileNotFoundError('Input file {} either does not exist or is not accessible'.format(args.gro))
    if args.selections is None:
        parser.print_help()
        raise ValueError('No selections specified. See help (python ndx.py -h)')
    fileChoicesOutput(args.output, parser)
    fileChoicesInput(args.gro, parser)


def main():
    # get & check args
    args = commandline()
    # make ndx
    n = Ndx(args.gro, ligands=args.ligands)
    for selection in args.selections:
        n.select(selection.lower())
    n.makeNDX(args.output)

if __name__ == '__main__':
    main()


        