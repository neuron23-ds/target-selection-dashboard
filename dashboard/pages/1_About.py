import streamlit as st
from streamlit_utils import *


st.set_page_config(
  page_title='Target Selection Dashboard',
  page_icon=n23_icon,
  layout='wide'
  )

add_logo()

# if not check_password():
#     st.stop()  # Do not continue if check_password is not True.

st.markdown(
    """
    # What information is available?
    ### Omics Assessment
    N23 data science's assessment of the omic data of the protein as\
    it pertains to the disease.
    
    ### Function
    As specified from UniProt.

    ### Location
    Tissue specificity and subcellular localization of the protein as\
    specificed by UniProt. 

    ### Protein Structure
    Does the protein have a 3D structure available via PDB? If so,\
    PDB ID is provided. If not, a BLAST has been run and top 5 results\
    are shown. 

    ### Existing drugs
    Does the protein have any existing drugs that target it? As specified\
    from ChEMBL.
    """,
    unsafe_allow_html=True,
)