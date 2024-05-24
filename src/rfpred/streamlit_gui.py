import streamlit as st
from streamlit_ketcher import st_ketcher
import matplotlib.pyplot as plt

from models import InputProcessing
from models import LightGBM_model

input_processor = InputProcessing()
model = LightGBM_model()

DEFAULT_MOL = (
    r"C"
)

solvent_options = {
    'DCM': 'ClCCl',
    'MeOH': 'CO',
    'MeCN': 'CC#N',
    'Toluene': 'Cc1ccccc1',
    'Hexane': 'CCCCCC',
    'Chloroform': 'ClC(Cl)Cl',
    'Acetone': 'CC(=O)C',
    'EtOH': 'CCO',
    'diethyl ether': 'CCOCC',
    'heptane': 'CCCCCCC',
    'petroleum ether (2-methylpentane)': 'CCCC(C)C',
    'triethylamine': 'CCN(CC)CC',
    'EtOAc': 'O=C(OCC)C',
    'THF': 'C1CCOC1'
}
solvents = ['DCM', 'MeOH', 'MeCN', 'Toluene', 'Hexane', 'Chloroform', 'Acetone', 'EtOH', 'diethyl ether', 'heptane', 'petroleum ether (2-methylpentane)', 'triethylamine', 'EtOAc', 'THF']


st.header('TLC Prediction', divider='rainbow')

with st.container(border=True):
    #streamlit_ketcher usage for crafting molecules
    molecule = st.text_input("Molecule", DEFAULT_MOL)
    smile_code = st_ketcher(molecule)
    st.markdown(f"current Smile code: ``{smile_code}``")

    with st.form("ped_form",border=False):
        st.write("choose solvents")
        
        #solvent selection boxes
        select_box_A = st.selectbox('Solvent A', list(solvents))
        select_box_B = st.selectbox('Solvent B', list(solvents))

        #slider
        select_proportion= st.slider('Pick a proportion A to B', 0, 100,50)

        #submit button for form
        submit_form = st.form_submit_button('get answer')

if submit_form:
    print("Input  -------------------------------------------------------")
    print("selected Molecule smile: " + smile_code)
    print("selected Solvent A: " + str(select_box_A))
    print("selected Solvent B: " + str(select_box_B))
    print(f"proportions: Solvent A {select_proportion}% to Solvent B {100-select_proportion}%")
    print("--------------------------------------------------------------")
    with st.container(border=True):
        with st.spinner('Wait for it...'):
            #nn model result ------------
            arr = input_processor.process_input(smiles=str(smile_code),
                                                solvent_A=select_box_A,
                                                solvent_B=select_box_B,
                                                percent_A=select_proportion)
            #arr = model_converter.convert_gui_input_to_nparray(str(smile_code),str(solvent_options[select_box_A]),str(solvent_options[select_box_B]),select_proportion,100-select_proportion)
            rf_value = model.predict(arr)
            
            
            st.success('Done!')
            #Values for Plotting
            x = [0]  
            y = [rf_value] 

            # create the Plot
            fig, ax = plt.subplots()
            
            #configure the plot
            ax.scatter(x, y, color='green', marker='o', s=100)  # draw spot

            #configure plot axis (hide x axis ticks)
            ax.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
            ax.set_ylim(0, 1)  # scales y-axis from0 to 100%

            #show plot
            st.pyplot(fig=fig, clear_figure=None, use_container_width=True)

            st.markdown(f"### Rf: ``{rf_value}``")

            

                        