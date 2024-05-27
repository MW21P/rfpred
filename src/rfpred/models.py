import os
from rdkit import Chem
import numpy as np
from rdkit.Chem import Descriptors, MACCSkeys
import lightgbm as lgb
    
class LightGBM_model:
    def __init__(self):
        """
        Sets parameter of lgb modal class.
        Loads the trained model from a txt file.
        
        Args:
            self: instance of the class
        """
        self.lgb_params = {
            'objective': 'regression',
            'num_leaves': 100,
            'learning_rate': 0.05,
            'n_estimators': 1000,
            'random_state': 1,
            'n_jobs': 1
        }
        self.model = None
        self.load_model()

    def load_model(self):
        """
        Loads the trained model from a txt file.
        
        Args:
            self: instance of the class
            
        Returns:
            lgb.model: The loaded lgb model
        """
        # Model class must be defined somewhere
        print('loading model...')
        dir_path = os.path.dirname(os.path.realpath(__file__))
        path = f'{dir_path}\\saved_models\\lgb_model_maccs_best_params.txt'
        # Load the model
        self.model= lgb.Booster(model_file=path)
        print(f'model loaded, path: {path}')
        return self.model
    
    def predict(self,input: np.array):
        """
        The trained model predicts the output of a given input.
        
        Args:
            self: instance of the class
            input (np.array): The feature matrix (input for model).
            
        Returns:
            lgb.model: None / predicted value 
        """
        if self.model != None:
            result = self.model.predict(input)
            if len(result) > 0:
                return round(float(result[0]), 2)
        return None
        



class InputProcessing:
    def __init__(self) -> None:
        pass
    def get_maccs(self,smiles):
        """
        Generate MACCS keys for a molecule from its SMILES string.
        
        Args:
            smiles (str): The SMILES string of the molecule.
            
        Returns:
            np.array: A numpy array of the MACCS keys.
        """
        maccs = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(smiles))
        return np.array([int(x) for x in list(maccs.ToBitString())])  # Convert MACCS keys to numpy array

    def get_solvent_features(self,solvent_A, solvent_B, percent_A, solvents):
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

    def get_rdkit_descriptors(self,smiles):
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

    def process_input(self,smiles, solvent_A, solvent_B, percent_A):
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
        maccs = self.get_maccs(smiles)
        solvent_features = self.get_solvent_features(solvent_A, solvent_B, percent_A, solvents)
        rdkit_descriptors = self.get_rdkit_descriptors(smiles)
        X = np.concatenate([maccs, rdkit_descriptors, solvent_features])
        x_2d = X[np.newaxis, :]
        return x_2d
