import os
import rfpred.models

class App:
    """A class to run a graphical interface."""
    def __init__(self) -> None:
        pass
    def run():
        """
        Starts the streamlit server.
        """
        #get current file path
        dir_path = os.path.dirname(os.path.realpath(__file__))
        #start the streamlit server runtime
        os.system(f"streamlit run {dir_path}\\streamlit_gui.py")

class Prediction:
    """A class to get a prediction."""
    
    def __init__(self):
        self.model = models.LightGBM_model()
        self.input_processor = models.InputProcessing()
        self.solvents = ['DCM', 'MeOH', 'MeCN', 'Toluene', 'Hexane', 'Chloroform', 'Acetone', 'EtOH', 'diethyl ether', 'heptane', 'petroleum ether (2-methylpentane)', 'triethylamine', 'EtOAc', 'THF'] 
    
    def predict(self,compound_smile, solvent_a, solvent_b, percent_solvent_a):
        """
        Process the input data into the model prediction. Retuns the predicted rf value as a float.
        
        Args:
            smiles (str): The SMILES string of the molecule.
            solvent_A (str): The first solvent.
            solvent_B (str): The second solvent.
            percent_A (float): The percentage of solvent A.
            
        Returns:
            float: The predicted rf value
        """
        arr = self.input_processor.process_input(smiles=compound_smile,
                                                 solvent_A=solvent_a,
                                                 solvent_B=solvent_b,
                                                 percent_A=percent_solvent_a)
        return float(self.model.predict(arr))

       
    
    