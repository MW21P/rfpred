import pytest
from rfpred.models import InputProcessing
from rfpred.rfpred_prediction import Prediction
import numpy as np

# Test InputProcessing class

def test_get_solvent_features():
    processor = InputProcessing()
    solvents = ['DCM', 'MeOH', 'MeCN', 'Toluene', 'Hexane', 'Chloroform', 'Acetone', 'EtOH', 'diethyl ether', 'heptane', 'petroleum ether (2-methylpentane)', 'triethylamine', 'EtOAc', 'THF']
    solvent_features = processor.get_solvent_features('DCM', 'MeOH', 50, solvents)
    assert len(solvent_features) == len(solvents), "Solvent features length should match solvents list length"

def test_get_solvent_features_invalid_solvents():
    processor = InputProcessing()
    solvents = ['DCM', 'MeOH', 'MeCN', 'Toluene']
    with pytest.raises(ValueError):
        processor.get_solvent_features('InvalidSolvent', 'MeOH', 50, solvents)

def test_get_rdkit_descriptors():
    processor = InputProcessing()
    smiles = "CCO"
    descriptors = processor.get_rdkit_descriptors(smiles)
    assert len(descriptors) == 4, "There should be 4 RDKit descriptors"

def test_process_input():
    processor = InputProcessing()
    smiles = "CCO"
    solvent_A = "DCM"
    solvent_B = "MeOH"
    percent_A = 50
    feature_matrix = processor.process_input(smiles, solvent_A, solvent_B, percent_A)
    assert feature_matrix.shape[1] == 182, "Feature matrix should have 182 features"

# Test Prediction class
def test_prediction():
    predictor = Prediction()
    smiles = "CCO"
    solvent_A = "DCM"
    solvent_B = "MeOH"
    percent_A = 50
    prediction = predictor.predict(smiles, solvent_A, solvent_B, percent_A)
    assert isinstance(prediction, float), "Prediction should be a float"