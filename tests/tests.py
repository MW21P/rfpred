from typing import Literal
import pytest
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import Descriptors

from rfpred.functions import (
    extract_rows_with_rf,
    canonicalise_smiles,
    get_solvents,
    convert_solvents,
    clean_smiles,
    are_enantiomers,
    get_maccs,
    get_solvent_features,
    get_rdkit_descriptors,
    process_input
)


def test_extract_rows_with_rf():
    """
    Test the extract_rows_with_rf function for filtering paragraphs containing Rf values.

    This function tests the extract_rows_with_rf function, which filters rows
    in a DataFrame that contain specified Rf values in the text. The test checks
    if the function correctly identifies and extracts rows containing valid Rf
    values from the given text.

    Test Case:
        - Input DataFrame contains various paragraphs, some of which contain Rf values 
          in different formats, and others do not contain any Rf values.

    Expected Output:
        - DataFrame containing only the rows with paragraphs that have valid Rf values.

    Asserts:
        - That the output DataFrame from the extract_rows_with_rf function is equal to
          the expected DataFrame containing only the relevant rows with Rf values.
    """
    
    # Test texts with Rf and without Rf values
    test = {'paragraphText': [
        'The value of the compound is', 
        'The compound has an RF value of 0.5', 
        'The measured Rf:0.5 of the product...', 
        'Dont show this f:0.5 value',
        'An RF(0.5) was measured',
        'Dont show this R: 0.5 value',
        'The Rf value of the compound is Rf=0.4',
        'The valueRf of the cpmound is 0.5',
        'The measured Rf : 0.5 of the product...',
        'Sometimes theRfvalue is hidden in the text',
        'And somtimes the Rf~0.6 is noted with a ~.'
    ]}
    
    df_test = pd.DataFrame(test)

    # Expected output of the function
    expected = {'paragraphText': [ 
        'The compound has an RF value of 0.5', 
        'The measured Rf:0.5 of the product...', 
        'An RF(0.5) was measured',
        'The Rf value of the compound is Rf=0.4',
        'The valueRf of the cpmound is 0.5',
        'The measured Rf : 0.5 of the product...',
        'Sometimes theRfvalue is hidden in the text',
        'And somtimes the Rf~0.6 is noted with a ~.'

    ]}

    df_expected = pd.DataFrame(expected)

    #running the function
    result = extract_rows_with_rf(df_test, 'paragraphText')

    #comparing the result of the function with the expected output
    pd.testing.assert_frame_equal(result, df_expected)

 
@pytest.mark.parametrize("smiles, expected", [
    ("C1=CC=CC=C1", "c1ccccc1"),  # Benzene
    ("C(C(=O)O)N", "NCC(=O)O"),  # Alanine
    ("O=C(O)CC(=O)O", "O=C(O)CC(=O)O"),  # Succinic acid
    ("CC(C)CC1=CC=CC=C1", "CC(C)Cc1ccccc1"),  # Isopropylbenzene
    ("N[C@@H](C)C(=O)O", "C[C@H](N)C(=O)O"),  # L-alanine
    ("N[C@H](C)C(=O)O", "C[C@@H](N)C(=O)O")  # D-alanine
])
def test_canonical_smiles(smiles: Literal['C1=CC=CC=C1'] | Literal['C(C(=O)O)N'] | Literal['O=C(O)CC(=O)O'] | Literal['CC(C)CC1=CC=CC=C1'] | Literal['N[C@@H](C)C(=O)O'] | Literal['N[C@H](C)C(=O)O'], expected: Literal['c1ccccc1'] | Literal['NCC(=O)O'] | Literal['O=C(O)CC(=O)O'] | Literal['CC(C)Cc1ccccc1'] | Literal['C[C@H](N)C(=O)O'] | Literal['C[C@@H](N)C(=O)O']):
    """
    Test the canonicalization of SMILES strings.

    This function checks if the given SMILES string is correctly converted to its
    canonical form.

    Args:
        smiles (str): The SMILES string to be canonicalized.
        expected (str): The expected canonical SMILES string.

    Asserts:
        That the function canonicalise_smiles returns the expected canonical SMILES string.
    """
    
    assert canonicalise_smiles(smiles) == expected


def test_empty_smiles():
    """
    Test the canonicalization of an empty SMILES string.

    This function checks if the canonicalization of an empty SMILES string
    correctly results in an empty string.

    Asserts:
        That the function canonicalise_smiles returns an empty string when given an empty string.
    """

    empty_smiles = ""
    assert canonicalise_smiles(empty_smiles) == ""


def test_invalide_smiles():
    """
    Test handling of an invalid SMILES string.

    This function checks if an invalid SMILES string correctly raises a
    ValueError exception.

    Asserts:
        That the function canonicalise_smiles raises a ValueError exception
        when given an invalid SMILES string.
        That the error message contains the expected message.
    """

    invalid_smiles = "invalid_smiles"
    
    with pytest.raises(ValueError) as e:
        canonicalise_smiles(invalid_smiles)
    assert str(e.value).startswith(f"Invalid SMILES string: {invalid_smiles}")


def test_get_solvents():
    """
    Test the extraction of unique solvents from a DataFrame.

    This function checks if the get_solvents function correctly identifies and
    extracts unique solvents from the provided DataFrame.

    The DataFrame contains various columns such as 'productSmiles', 'Rf', 
    'Solvent_A', 'Solvent_B', 'Percent_A', 'Percent_B', 'Additive_C', 
    and 'Percent_C'. The function should identify all unique solvents listed 
    in the 'Solvent_A' and 'Solvent_B' columns.

    Test Data:
        - 'productSmiles': List of product SMILES strings.
        - 'Rf': List of Rf values.
        - 'Solvent_A': Primary solvents.
        - 'Solvent_B': Secondary solvents.
        - 'Percent_A': Percentage of primary solvent.
        - 'Percent_B': Percentage of secondary solvent.
        - 'Additive_C': Additives used.
        - 'Percent_C': Percentage of additives.

    Expected Output:
        - A sorted list of unique solvents.

    Asserts:
        That the function get_solvents returns the expected list of unique solvents,
        sorted in alphabetical order.
    """

    test = {
            'productSmiles': ["['CCO']", "['CC(=O)O']", "['CCN(CC)CC']", "['CCCCCCCC']", "['CCCCCCC']", "['CCCCCC']", "['CCCCC']", "['CCCC']", "['CCC']", "['CC']", "['C']", "['C1CCCCC1']", "['CC(C)C(=O)O']"],
            'Rf': ['0.5', '0.77', '0.1', '1.3', '0.6', '0.8', '0.9', '0.4', '0.7', '0.3', '0.2', None, '0.33'],
            'Solvent_A': ['ethyl acetate', None, 'hexane', 'hexane', 'Methanol', 'Hexane', 'Ethyl acetate', 'Ethyl acetate', 'hexane', 'Ethyl acetate', 'hexane', 'DCM', None],
            'Solvent_B': ['hexane', 'hexane', 'ethyl acetate', None, 'Ethanol', 'Ethyl acetate', 'Tetrahydrofuran', 'hexane', 'DCM', 'hexane', 'DCM', 'MeOH', None],
            'Percent_A': ['50', None, '20', '20', '70', '10', '100', '40', '66.666666', '60', '50', '80', '20'],
            'Percent_B': ['50', '30', None, '80', '30', '25', None, '80', '33.333333', '40', '45', '20', '80'],
            'Additive_C': [None, None, None, None, 'TEA', None, None, None, None, 'TEA', None, None, None],
            'Percent_C': [None, None, None, None, '5', None, None, None, None, None, '5', None, None]
        }
    df_test = pd.DataFrame(test)

    expected = ['DCM', 'Ethanol', 'Ethyl acetate', 'ethyl acetate', 'Hexane', 'hexane', 'MeOH', 'Methanol', 'Tetrahydrofuran']

    result = get_solvents(df_test)

    assert sorted(result) == sorted(expected)


def test_convert_solvents():
    """
    Test the conversion and sorting of solvents in a DataFrame.

    This function checks if the `convert_solvents` function correctly converts
    solvent names in columns 'Solvent_A' and 'Solvent_B' to their standardized
    forms and sorts out the rows with non-standard solvent names.

    Test Data:
        - 'Solvent_A': List of solvents in non-standard names.
        - 'Solvent_B': List of solvents in non-standard names.

    Expected Output:
        - A DataFrame with standardized solvent names and their SMILES strings.
        - A DataFrame with rows that could not be converted.

    Asserts:
        - That the function `convert_solvents` returns the expected DataFrame 
          with standardized solvent names and SMILES strings.
        - That the function `convert_solvents` returns the expected DataFrame 
          with sorted out rows containing non-standard solvent names.
    """

    test = {
        'Solvent_A': ['ethyl acetate', 'Hexane', 'methylene chloride', 'mecn', 'n-hex', 'acetone', 'ethanol', 'n-heptane', 'et3n', 'none', '20%Etahnol', 'none', 'EA/HE'],
        'Solvent_B': ['hexane', 'EA', 'acetonitril', 'ch3oh', 'toluene', 'trichloromethane', 'et2o', 'pe', 'tetrahydrofuran', 'EA', 'et2o', '(TEA)', 'THF']
    }

    df_test = pd.DataFrame(test)

    expected = {
        'Solvent_A': ['EtOAc', 'Hexane', 'DCM', 'MeCN', 'Hexane', 'Acetone', 'EtOH', 'heptane', 'triethylamine', 'none'],
        'Solvent_B': ['Hexane', 'EtOAc', 'MeCN', 'MeOH', 'Toluene', 'Chloroform', 'diethyl ether', 'petroleum ether (2-methylpentane)', 'THF', 'EtOAc' ],
        'Solvent_A_Smiles': ['O=C(OCC)C', 'CCCCCC', 'ClCCl', 'CC#N', 'CCCCCC', 'CC(=O)C', 'CCO', 'CCCCCCC', 'CCN(CC)CC', None],
        'Solvent_B_Smiles': ['CCCCCC', 'O=C(OCC)C', 'CC#N', 'CO', 'Cc1ccccc1', 'ClC(Cl)Cl', 'CCOCC', 'CCCC(C)C', 'C1CCOC1', 'O=C(OCC)C']
    }

    df_expected = pd.DataFrame(expected)

    sorted_out = {
        'Solvent_A': ['20%Etahnol', 'none', 'EA/HE'],
        'Solvent_B': ['et2o', '(TEA)', 'THF']
    }   

    df_sorted_out_expected = pd.DataFrame(sorted_out)

    result, df_sorted_out = convert_solvents(df_test, 'Solvent_A', 'Solvent_B')

    pd.testing.assert_frame_equal(df_expected, result)

    pd.testing.assert_frame_equal(df_sorted_out_expected, df_sorted_out)


def clean_smiles():
    """
    Test the cleaning of SMILES strings in a DataFrame.

    This function checks if the `clean_smiles` function correctly removes leading
    and trailing whitespace, and extraneous quotes from SMILES strings in a given
    column of a DataFrame.

    Test Data:
        - 'Smiles': List of SMILES strings with leading/trailing whitespace and quotes.

    Expected Output:
        - A DataFrame with cleaned SMILES strings, where leading/trailing whitespace
          and quotes have been removed.

    Asserts:
        - That the function `clean_smiles` returns the expected DataFrame with cleaned
          SMILES strings.
    """
    
    test = {
        'Smiles': ['  COO', 'CCC  ', '"COC"', 'CC', ' COO"', '" CCCC', '   C   ']
    }

    df_test = pd.DataFrame(test)

    expected = {
        'Smiles': ['COO', 'CCC', 'COC', 'CC', 'COO', 'CCCC', 'C']
    }

    df_expected = pd.DataFrame(expected)

    result = clean_smiles(df_test, 'Smiles')

    pd.testing.assert_frame_equal(df_expected, result)


def test_are_enantiomers():
    """
    Test the are_enantiomers function for various cases.

    This function tests the are_enantiomers function, which checks if two
    molecules represented by SMILES strings are enantiomers and returns
    a tuple (True, index) if they are, or (False, None) if they are not.

    Tests:
        1. Known enantiomers:
            - L-alanine and D-alanine should be recognized as enantiomers.
        2. Known non-enantiomers:
            - Ethanol and Butane should not be recognized as enantiomers.
        3. Identical molecules:
            - Two identical Ethanol molecules should not be recognized as enantiomers.
        4. Invalid SMILES strings:
            - Any comparison involving invalid SMILES strings should return (False, None).
        5. Additional enantiomer pair:
            - Additional test case for another known pair of enantiomers.

    Asserts:
        - That the are_enantiomers function returns the correct tuple for each test case.
    """

    # Test with known enantiomers
    smiles_A = 'N[C@@H](C)C(=O)O'  # L-alanine
    smiles_B = 'N[C@H](C)C(=O)O'  # D-alanine
    assert are_enantiomers(smiles_A, smiles_B) == (True, 1)

    # Test with known non-enantiomers
    smiles_A = 'CCO'  # Ethanol
    smiles_B = 'CCCC'  # Butane
    assert are_enantiomers(smiles_A, smiles_B) == (False, None)

    # Test with the same molecule
    smiles_A = 'CCO'  # Ethanol
    smiles_B = 'CCO'  # Ethanol
    assert are_enantiomers(smiles_A, smiles_B) == (False, None)

    # Test with invalid SMILES
    smiles_A = 'invalid_smiles'
    smiles_B = 'CCO'  # Ethanol
    assert are_enantiomers(smiles_A, smiles_B) == (False, None)

    smiles_A = 'CCO'  # Ethanol
    smiles_B = 'invalid_smiles'
    assert are_enantiomers(smiles_A, smiles_B) == (False, None)

    smiles_A = 'invalid_smiles'
    smiles_B = 'invalid_smiles'
    assert are_enantiomers(smiles_A, smiles_B) == (False, None)

    smiles_A = 'CC(O)[C@@H](N)C(=O)O'
    smiles_B = 'CC(O)[C@H](N)C(=O)O'
    assert are_enantiomers(smiles_A, smiles_B) == (True, 3)


def test_get_maccs():
    """
    Test the get_maccs function for correctness.

    This function tests the get_maccs function, which generates the MACCS keys
    fingerprint for a given molecule represented by a SMILES string. The test
    compares the output of the get_maccs function with the expected MACCS keys
    fingerprint generated using RDKit.

    Tests:
        1. Known molecule (Ethanol):
            - SMILES string: 'CCO'
            - Expected MACCS keys fingerprint generated using RDKit.
        2. Another known molecule (Butane):
            - SMILES string: 'CCCC'
            - Expected MACCS keys fingerprint generated using RDKit.

    Asserts:
        - That the get_maccs function returns the correct MACCS keys fingerprint
          array for each test case.
    """

    # Test with a known molecule
    smiles = 'CCO'  # Ethanol
    expected_maccs = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(smiles))
    expected_array = np.array([int(x) for x in list(expected_maccs.ToBitString())])
    result = get_maccs(smiles)
    np.testing.assert_array_equal(result, expected_array)

    # Test with another known molecule
    smiles = 'CCCC'  # Butane
    expected_maccs = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(smiles))
    expected_array = np.array([int(x) for x in list(expected_maccs.ToBitString())])
    result = get_maccs(smiles)
    np.testing.assert_array_equal(result, expected_array)


def test_get_solvent_features():
    """
    Test the get_solvent_features function for correctness.

    This function tests the get_solvent_features function, which generates a feature
    array representing the composition of two solvents in a mixture, given their names
    and the percentage of the first solvent. The test covers various scenarios with
    different solvents and percentages.

    Tests:
        1. Known solvents and percentages:
            - Solvent A: 'water'
            - Solvent B: 'ethanol'
            - Percent A: 60.0
            - Expected feature array: [60.0, 40.0, 0.0, 0.0]
        
        2. Different combination:
            - Solvent A: 'methanol'
            - Solvent B: 'acetone'
            - Percent A: 30.0
            - Expected feature array: [0.0, 0.0, 30.0, 70.0]
        
        3. Solvent A not in the list:
            - Solvent A: 'hexane'
            - Solvent B: 'ethanol'
            - Percent A: 50.0
            - Expected feature array: [0.0, 50.0, 0.0, 0.0]
        
        4. Solvent B not in the list:
            - Solvent A: 'water'
            - Solvent B: 'hexane'
            - Percent A: 80.0
            - Expected feature array: [80.0, 0.0, 0.0, 0.0]
        
        5. Both solvents not in the list:
            - Solvent A: 'hexane'
            - Solvent B: 'benzene'
            - Percent A: 20.0
            - Expected feature array: [0.0, 0.0, 0.0, 0.0]
        
        6. Percent A as 0:
            - Solvent A: 'ethanol'
            - Solvent B: 'acetone'
            - Percent A: 0.0
            - Expected feature array: [0.0, 0.0, 0.0, 100.0]
        
        7. Percent A as 100:
            - Solvent A: 'water'
            - Solvent B: 'methanol'
            - Percent A: 100.0
            - Expected feature array: [100.0, 0.0, 0.0, 0.0]

    Asserts:
        - That the get_solvent_features function returns the correct feature array for
          each test case.
    """

    solvents = ['water', 'ethanol', 'methanol', 'acetone']

    # Test with known solvents and percentages
    solvent_A = 'water'
    solvent_B = 'ethanol'
    percent_A = 60.0
    expected_feature = np.array([60.0, 40.0, 0.0, 0.0])
    result = get_solvent_features(solvent_A, solvent_B, percent_A, solvents)
    np.testing.assert_array_equal(result, expected_feature)

    # Test with a different combination
    solvent_A = 'methanol'
    solvent_B = 'acetone'
    percent_A = 30.0
    expected_feature = np.array([0.0, 0.0, 30.0, 70.0])
    result = get_solvent_features(solvent_A, solvent_B, percent_A, solvents)
    np.testing.assert_array_equal(result, expected_feature)

    # Test with solvent_A not in the list
    solvent_A = 'hexane'
    solvent_B = 'ethanol'
    percent_A = 50.0
    expected_feature = np.array([0.0, 50.0, 0.0, 0.0])
    result = get_solvent_features(solvent_A, solvent_B, percent_A, solvents)
    np.testing.assert_array_equal(result, expected_feature)

    # Test with solvent_B not in the list
    solvent_A = 'water'
    solvent_B = 'hexane'
    percent_A = 80.0
    expected_feature = np.array([80.0, 0.0, 0.0, 0.0])
    result = get_solvent_features(solvent_A, solvent_B, percent_A, solvents)
    np.testing.assert_array_equal(result, expected_feature)

    # Test with both solvents not in the list
    solvent_A = 'hexane'
    solvent_B = 'benzene'
    percent_A = 20.0
    expected_feature = np.array([0.0, 0.0, 0.0, 0.0])
    result = get_solvent_features(solvent_A, solvent_B, percent_A, solvents)
    np.testing.assert_array_equal(result, expected_feature)

    # Test with percent_A as 0
    solvent_A = 'ethanol'
    solvent_B = 'acetone'
    percent_A = 0.0
    expected_feature = np.array([0.0, 0.0, 0.0, 100.0])
    result = get_solvent_features(solvent_A, solvent_B, percent_A, solvents)
    np.testing.assert_array_equal(result, expected_feature)

    # Test with percent_A as 100
    solvent_A = 'water'
    solvent_B = 'methanol'
    percent_A = 100.0
    expected_feature = np.array([100.0, 0.0, 0.0, 0.0])
    result = get_solvent_features(solvent_A, solvent_B, percent_A, solvents)
    np.testing.assert_array_equal(result, expected_feature)

def test_get_rdkit_descriptors():
    """
    Test the get_rdkit_descriptors function for correctness.

    This function tests the get_rdkit_descriptors function, which calculates a set of
    molecular descriptors using RDKit for a given SMILES string. The test covers 
    several known molecules and compares the calculated descriptors to expected values.

    Test Cases:
        1. Wasser (H2O):
            - SMILES: 'O'
            - Expected descriptors:
                - Molecular Weight
                - LogP
                - Number of Hydrogen Donors
                - Number of Hydrogen Acceptors
        
        2. Benzene (C6H6):
            - SMILES: 'c1ccccc1'
            - Expected descriptors:
                - Molecular Weight
                - LogP
                - Number of Hydrogen Donors
                - Number of Hydrogen Acceptors
        
        3. Ethanol (C2H6O):
            - SMILES: 'CCO'
            - Expected descriptors:
                - Molecular Weight
                - LogP
                - Number of Hydrogen Donors
                - Number of Hydrogen Acceptors

    Asserts:
        - That the get_rdkit_descriptors function returns the correct descriptor array 
          for each test case, with values almost equal to the expected ones (up to 6
          decimal places).
    """

    # Testcase 1: Wasser (H2O)
    smiles = 'O'
    expected = np.array([
        Descriptors.MolWt(Chem.MolFromSmiles(smiles)),
        Descriptors.MolLogP(Chem.MolFromSmiles(smiles)),
        Descriptors.NumHDonors(Chem.MolFromSmiles(smiles)),
        Descriptors.NumHAcceptors(Chem.MolFromSmiles(smiles))
    ])
    result = get_rdkit_descriptors(smiles)
    np.testing.assert_array_almost_equal(result, expected, decimal=6)

    # Testcase 2: Benzene (C6H6)
    smiles = 'c1ccccc1'
    expected = np.array([
        Descriptors.MolWt(Chem.MolFromSmiles(smiles)),
        Descriptors.MolLogP(Chem.MolFromSmiles(smiles)),
        Descriptors.NumHDonors(Chem.MolFromSmiles(smiles)),
        Descriptors.NumHAcceptors(Chem.MolFromSmiles(smiles))
    ])
    result = get_rdkit_descriptors(smiles)
    np.testing.assert_array_almost_equal(result, expected, decimal=6)

    # Testcase 3: Ethanol (C2H6O)
    smiles = 'CCO'
    expected = np.array([
        Descriptors.MolWt(Chem.MolFromSmiles(smiles)),
        Descriptors.MolLogP(Chem.MolFromSmiles(smiles)),
        Descriptors.NumHDonors(Chem.MolFromSmiles(smiles)),
        Descriptors.NumHAcceptors(Chem.MolFromSmiles(smiles))
    ])
    result = get_rdkit_descriptors(smiles)
    np.testing.assert_array_almost_equal(result, expected, decimal=6)

def test_process_input():
    """
    Test the process_input function for correctness.

    This function tests the process_input function, which processes the input 
    consisting of a SMILES string, two solvents, and a percentage of solvent A. 
    The test checks if the function correctly combines MACCS keys, RDKit 
    descriptors, and solvent features into a single feature array.

    Test Case:
        - SMILES: 'CCO' (Ethanol)
        - Solvent A: 'DCM' (Dichloromethane)
        - Solvent B: 'MeOH' (Methanol)
        - Percentage of Solvent A: 60.0

    Expected Output:
        - A concatenated array of MACCS keys, RDKit descriptors, and solvent 
          features for the given input.

    Asserts:
        - That the output of the process_input function is almost equal (up to 6 
          decimal places) to the expected concatenated feature array.
    """

# Example input
    smiles = 'CCO'  # Ethanol
    solvent_A = 'DCM'
    solvent_B = 'MeOH'
    percent_A = 60.0

    # Expected output
    solvents = ['DCM', 'MeOH', 'MeCN', 'Toluene', 'Hexane', 'Chloroform', 'Acetone', 'EtOH', 'diethyl ether', 'heptane', 'petroleum ether (2-methylpentane)', 'triethylamine', 'EtOAc', 'THF']
    maccs = get_maccs(smiles)
    solvent_features = get_solvent_features(solvent_A, solvent_B, percent_A, solvents)
    rdkit_descriptors = get_rdkit_descriptors(smiles)
    expected = np.concatenate([maccs, rdkit_descriptors, solvent_features])

    # output of the function
    result = process_input(smiles, solvent_A, solvent_B, percent_A)

    # Test if the output is as expected
    np.testing.assert_array_almost_equal(result, expected, decimal=6)
