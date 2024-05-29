import pandas as pd
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import Descriptors
import re


def extract_rows_with_rf(Dataframe: pd.DataFrame, column_name: str):
    """Function that extracts rows with Rf values from a dataframe and returns a new dataframe containing only these rows. 

    Args:
        Dataframe (_type_): Dataframe containing the extracted data from the US patents
    """
    # copy the dataframe to leave old dataframe unchanged
    df = Dataframe.copy()
    
    # Define the regex patterns
    Rf_check = r'( ?R[fF][ :=(]?)'
    
    # List to store indices of rows without Rf values
    rows_to_drop = []
    
    # Search for rows with Rf values in the paragraphText column
    for index, row in df.iterrows():
        checkRf = re.findall(Rf_check, row[column_name])

        if not checkRf:
            rows_to_drop.append(index)
               
    # Drop rows without Rf values
    df = df.drop(rows_to_drop).reset_index(drop=True)
            
    return df



def canonicalise_smiles(Smiles: str):
    '''
       Converts Smile to a Mol file and back to a Smiles again to create
       a consistent Smiles string.

       Args: Smiles string
    '''
    mol = Chem.MolFromSmiles(Smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {Smiles}")
    can_Smiles = Chem.MolToSmiles(mol)
    return can_Smiles



def get_solvents(Dataframe: pd.DataFrame):
    """Function that iterates over all solvent entries in the dataframe and checks, if the solvent mentioned is already in the list. If not, it adds the solvent to the list.

    Args:
        Dataframe (_type_): dataframe with the extracted data from the US patents, preprocessed with the get_values function.
    """
    solvents = []
    
    for index, row in Dataframe.iterrows():
        solvent_a = row['Solvent_A']
        solvent_b = row['Solvent_B']
        
        if solvent_a not in solvents and solvent_a is not None:
            solvents.append(solvent_a)
        if solvent_b not in solvents and solvent_b is not None:
            solvents.append(solvent_b)
        else:
            continue
        
    return solvents



def convert_solvents(Dataframe: pd.DataFrame, column_solvent_A: str, column_solvent_B: str):
    """Converts the solvent names in the dataframe to a standardized form.

    Args:
        Dataframe (pd.DataFrame): Dataframe containing the extracted data from the US patents, 
                                preprocessed with the get_values and clean_up function.
                                
    Returns:
        Dataframe (pd.DataFrame): Dataframe with standardized solvent names and SMILES.
        df_sorted_out (pd.DataFrame): Dataframe containing rows with solvents that could not be converted.
    """
    Dataframe = Dataframe.copy()
    df_sorted_out = pd.DataFrame()
    indices_to_drop = []
    size_pre_conversion = Dataframe.shape[0]  # get the size of the dataframe
    
    for index, row in Dataframe.iterrows():
        for solvent_column in [column_solvent_A, column_solvent_B]: # iterate through the solvent columns
            solvent_name = row[solvent_column].lower()
            
            if '/' in solvent_name or '%' in solvent_name or '(' in solvent_name or ')' in solvent_name: # if there is a mixture of solvents, take out the entire row
                df_sorted_out = pd.concat([df_sorted_out, pd.DataFrame(row).transpose()], ignore_index=True)
                indices_to_drop.append(index) # add the index to the list of indices to drop not to change size of Dataframe with every iteration
                
            elif 'dcm' in solvent_name or 'dichlorometh' in solvent_name or 'methylene chloride' in solvent_name or 'ch2cl2' in solvent_name and 'chloroform' not in solvent_name and 'trichloromethane' not in solvent_name: # exclude chloroform and trichloromethane for substring search
                Dataframe.at[index, solvent_column] = 'DCM'     
                Dataframe.at[index, solvent_column + '_Smiles'] = 'ClCCl' # add SMILES of the solvent A 
                
            elif 'meoh' in solvent_name or 'methanol' in solvent_name or 'ch3oh' in solvent_name:
                Dataframe.at[index, solvent_column] = 'MeOH'
                Dataframe.at[index, solvent_column + '_Smiles'] = 'CO'
                
            elif 'acetonitrile' in solvent_name or 'acetonitril' in solvent_name or 'mecn' in solvent_name or 'ch3cn' in solvent_name or 'acn' in solvent_name:
                Dataframe.at[index, solvent_column] = 'MeCN'
                Dataframe.at[index, solvent_column + '_Smiles'] = 'CC#N'
                
            elif 'toluene' in solvent_name:
                Dataframe.at[index, solvent_column] = 'Toluene'
                Dataframe.at[index, solvent_column + '_Smiles'] = 'Cc1ccccc1'
                
            elif 'hexane' in solvent_name or 'hex' in solvent_name or 'n-hex' in solvent_name or 'hexanes' in solvent_name and 'cycl' not in solvent_name: # exclude cyclohexane
                Dataframe.at[index, solvent_column] = 'Hexane'
                Dataframe.at[index, solvent_column + '_Smiles'] = 'CCCCCC'
                
            elif 'chloroform' in solvent_name or 'trichloromethane' in solvent_name or 'chcl3' in solvent_name:
                Dataframe.at[index, solvent_column] = 'Chloroform'
                Dataframe.at[index, solvent_column + '_Smiles'] = 'ClC(Cl)Cl'
                
            elif 'acetone' in solvent_name or 'me2co' in solvent_name or 'ch3coch3' in solvent_name or '2-propanone' in solvent_name or 'dimethyl ketone' in solvent_name:
                Dataframe.at[index, solvent_column] = 'Acetone'
                Dataframe.at[index, solvent_column + '_Smiles'] = 'CC(=O)C'
                
            elif 'etoh' in solvent_name or 'ethanol' in solvent_name or 'ch3ch2oh' in solvent_name or 'alcohol' in solvent_name:
                Dataframe.at[index, solvent_column] = 'EtOH'
                Dataframe.at[index, solvent_column + '_Smiles'] = 'CCO'
                
            elif 'ether' in solvent_name or 'diethylether' in solvent_name or 'et2o' in solvent_name or 'diethyl ether' in solvent_name:
                Dataframe.at[index, solvent_column] = 'diethyl ether'
                Dataframe.at[index, solvent_column + '_Smiles'] = 'CCOCC'
                
            elif 'hep' in solvent_name or 'heptane' in solvent_name or 'n-heptane' in solvent_name or 'heptanes' in solvent_name:
                Dataframe.at[index, solvent_column] = 'heptane'
                Dataframe.at[index, solvent_column + '_Smiles'] = 'CCCCCCC'
                
            elif 'petroleum ether' in solvent_name or 'light petroleum' in solvent_name or 'pe' in solvent_name or 'pet ether' in solvent_name or 'pet. ether' in solvent_name or 'pet.ether' in solvent_name or 'petroleum ether (bp 40-60)' in solvent_name or 'petroleum ether (bp 60-80)' in solvent_name or 'petroleum ether (bp 80-100)' in solvent_name or 'petroleum ether (bp 100-140)' in solvent_name or 'petroleum ether (bp 140-180)' in solvent_name or 'petroleum ether (bp 180-210)' in solvent_name or 'petroleum ether (bp 210-240)' in solvent_name or 'petroleum ether (bp 240-270)' in solvent_name or 'petroleum ether (bp 270-300)' in solvent_name or 'petroleum ether (bp 300-400)' in solvent_name or 'petroleum ether (bp 400-600)' in solvent_name or 'petroleum ether (bp 600-800)' in solvent_name or 'petroleum ether (bp 800-1000)' in solvent_name or 'petroleum ether (bp 1000-1200)' in solvent_name or 'petroleum ether (bp 1200-1400)' in solvent_name or 'petroleum ether (bp 1400-1600)' in solvent_name or 'petroleum ether (bp 1600-1800)' in solvent_name or 'petroleum ether (bp 1800-2000)' in solvent_name or 'petroleum ether (bp 2000-2200)' in solvent_name or 'petroleum ether (bp 2200-2400)' in solvent_name or 'petroleum ether (bp 2400-2600)' in solvent_name or 'petroleum ether (bp 2600-2800)' in solvent_name or 'petroleum ether (bp 2800-3000)' in solvent_name or 'petroleum ether (bp 3000-3200)' in solvent_name or 'petroleum ether (bp 3200-3400)' in solvent_name or 'petroleum ether (bp 3400-3600)' in solvent_name or 'petroleum ether (bp 3600-3800)' in solvent_name or 'petroleum ether (bp 3800-4000)' in solvent_name or 'petroleum ether (bp 4000-4200)' in solvent_name:
                Dataframe.at[index, solvent_column] = 'petroleum ether (2-methylpentane)'
                Dataframe.at[index, solvent_column + '_Smiles'] = 'CCCC(C)C'
                
            elif 'triethylamine' in solvent_name or 'et3n' in solvent_name or 'tea' in solvent_name:
                Dataframe.at[index, solvent_column] = 'triethylamine'
                Dataframe.at[index, solvent_column + '_Smiles'] = 'CCN(CC)CC'
                
            elif 'etoac' in solvent_name or 'ethyl acetate' in solvent_name or 'ethylacetate' in solvent_name or 'ethylacetat' in solvent_name or 'ea' in solvent_name or 'ethanol acetate' in solvent_name:
                Dataframe.at[index, solvent_column] = 'EtOAc'
                Dataframe.at[index, solvent_column + '_Smiles'] = 'O=C(OCC)C' 
                
            elif 'thf' in solvent_name or 'tetrahydrofuran' in solvent_name:
                Dataframe.at[index, solvent_column] = 'THF'
                Dataframe.at[index, solvent_column + '_Smiles'] = 'C1CCOC1'
                
            elif 'none' in solvent_name:
                Dataframe.at[index, solvent_column + '_Smiles'] = None
                
            else: # if the solvent is not in the list, add it to the another dataframe but still sort it out. We don't want things to get too complex.
                df_sorted_out = pd.concat([df_sorted_out, pd.DataFrame(row).transpose()], ignore_index=True)
                indices_to_drop.append(index) # add the index to the list of indices to drop not to change size of Dataframe with every iteration

    Dataframe = Dataframe.drop(indices_to_drop)
    Dataframe.reset_index(drop=True, inplace=True)
    df_sorted_out.reset_index(drop=True, inplace=True)
    
    size_post_conversion = Dataframe.shape[0]  # get the size of the dataframe after conversion to see how many rows were dropped
    print(f"Size of the dataframe before conversion: {size_pre_conversion}")
    print(f"Size of the dataframe after conversion: {size_post_conversion}")
    print(f"Number of rows dropped: {size_pre_conversion - size_post_conversion}")
    print(f"Percentage of rows dropped: {(size_pre_conversion - size_post_conversion) / size_pre_conversion * 100}%")
        
    return Dataframe, df_sorted_out



def clean_smiles(df: pd.DataFrame, column_name: str):
    """Cleanes SMILES as apparently there are quotes around the SMILES which makes it impossible for rDKit to parse.

    Args:
        df (pd.DataFrame): dataframe cleaned with all previous functions
        column_name (str): productSmiles column

    Returns:
        df: pd.DataFrame: cleaned dataframe with SMILES that can be parsed
    """
    # Remove any leading or trailing whitespace
    df[column_name] = df[column_name].str.strip()
    
    # Remove any single or double quotes
    df[column_name] = df[column_name].str.replace("'", "").replace('"', '')
    
    return df



def are_enantiomers(smi1: str, smi2: str):
    """
    Function that checks if two molecules are enantiomers. Credits to: https://github.com/rdkit/rdkit/discussions/7169 for the inspiration with the swapping of the "@" and "@@" in the SMILES strings.

    Args:
        smi1 (str): SMILES string of molecule 1
        smi2 (str): SMILES string of molecule 2

    Returns:
        tuple (Bool, int): (True, index) if the molecules are enantiomers, (False, None) otherwise
    """
    # Convert the SMILES strings to RDKit molecule objects
    mol1 = Chem.MolFromSmiles(smi1)
    mol2 = Chem.MolFromSmiles(smi2)

    # If either SMILES string is invalid, return False
    if mol1 is None or mol2 is None:
        return False, None

    # Convert the molecule objects back to canonical SMILES strings
    can_smi1 = Chem.MolToSmiles(mol1)
    can_smi2 = Chem.MolToSmiles(mol2)

    # If either SMILES string does not contain "@", return False
    if "@" not in can_smi1 or "@" not in can_smi2:
        return False, None

    # Swap "@" and "@@" in can_smi1
    swapped_can_smi1 = can_smi1.replace("@@", "__DOUBLE_AT__").replace("@", "@@").replace("__DOUBLE_AT__", "@")

    # If the swapped can_smi1 is equal to can_smi2, return True
    if swapped_can_smi1 == can_smi2:
        # Find the index of the chiral atom
        for atom in mol1.GetAtoms():
            if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                return True, atom.GetIdx()
    
    # If none of the above conditions were met, return False
    return False, None


def get_maccs(smiles: str):
    """
    Generate MACCS keys for a molecule from its SMILES string.
    
    Args:
        smiles (str): The SMILES string of the molecule.
        
    Returns:
        np.array: A numpy array of the MACCS keys.
    """
    maccs = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(smiles))
    return np.array([int(x) for x in list(maccs.ToBitString())])  # Convert MACCS keys to numpy array


def get_solvent_features(solvent_A: str, solvent_B: str, percent_A: float, solvents: list):
    """
    Generate solvent features for a given solvent combination.
    
    Args:
        solvent_A (str): The first solvent.
        solvent_B (str): The second solvent.
        percent_A (float): The percentage of solvent A.
        solvents (list): The list of all possible solvents.
        
    Returns:
        np.array: A numpy array of the solvent features.
    """
    percent_B = 100 - percent_A
    solvent_feature = np.zeros(len(solvents))
    if solvent_A in solvents:
        solvent_A_index = solvents.index(solvent_A)
        solvent_feature[solvent_A_index] = percent_A
    if solvent_B in solvents:
        solvent_B_index = solvents.index(solvent_B)
        solvent_feature[solvent_B_index] = percent_B
    return solvent_feature

def get_rdkit_descriptors(smiles: str):
    """
    Generate RDKit descriptors for a molecule from its SMILES string.
    
    Args:
        smiles (str): The SMILES string of the molecule.
        
    Returns:
        np.array: A numpy array of the RDKit descriptors.
    """
    mol = Chem.MolFromSmiles(smiles)
    return np.array([Descriptors.MolWt(mol), Descriptors.MolLogP(mol),
                     Descriptors.NumHDonors(mol), Descriptors.NumHAcceptors(mol)])

    
def process_input(smiles: str, solvent_A: str, solvent_B: str, percent_A: float):
    """
    Process the input data into a feature matrix.
    
    Args:
        smiles (str): The SMILES string of the molecule.
        solvent_A (str): The first solvent.
        solvent_B (str): The second solvent.
        percent_A (float): The percentage of solvent A.
        
    Returns:
        np.array: The feature matrix.
    """
    solvents = ['DCM', 'MeOH', 'MeCN', 'Toluene', 'Hexane', 'Chloroform', 'Acetone', 'EtOH', 'diethyl ether', 'heptane', 'petroleum ether (2-methylpentane)', 'triethylamine', 'EtOAc', 'THF']
    maccs = get_maccs(smiles)
    solvent_features = get_solvent_features(solvent_A, solvent_B, percent_A, solvents)
    rdkit_descriptors = get_rdkit_descriptors(smiles)
    X = np.concatenate([maccs, rdkit_descriptors, solvent_features])
    return X