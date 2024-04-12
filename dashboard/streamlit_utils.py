from pathlib import Path
import streamlit as st
import hmac
import numpy as np


def check_password():
    """Returns `True` if the user had the correct password."""

    def password_entered():
        """Checks whether a password entered by the user is correct."""
        if hmac.compare_digest(st.session_state["password"], st.secrets["password"]):
            st.session_state["password_correct"] = True
            del st.session_state["password"]  # Don't store the password.
        else:
            st.session_state["password_correct"] = False

    # Return True if the password is validated.
    if st.session_state.get("password_correct", False):
        return True

    # Show input for password.
    st.text_input(
        "Enter institutional login password", type="password", on_change=password_entered, key="password"
    )
    if "password_correct" in st.session_state:
        st.error("😕 Password incorrect")
    return False

from PIL import Image
n23_icon = Image.open('dashboard/n23_icon.png')


# Thanks to https://discuss.streamlit.io/t/put-logo-and-title-above-on-top-of-page-navigation-in-sidebar-of-multipage-app/28213/5
def add_logo():
    st.markdown(
        """
        <style>
            [data-testid="stSidebarNav"] {
                background-image: url(https://neuron23.com/wp-content/themes/neuron23/assets/images/logo.svg);
                background-repeat: no-repeat;
                padding-top: 0px;
                background-position: 20px 20px;
                margin-top: 20px;
            }
            [data-testid="stSidebarNav"]::before {
                content: "Target Selection Dashboard";
                font-size: 25px;
                margin-left: 20px;
                top: 80px;
                position: relative;
            }
        </style>
        """,
        unsafe_allow_html=True,
    )

# Thanks to https://discuss.streamlit.io/t/increase-expanders-label-text-size/36227/4
def local_css(file_name):
    with open(file_name) as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

# Dictionary of EC codes for classification
ec_codes = {
    1:'Oxioreductase',
    2:'Transferase',
    3:'Hydrolase',
    4:'Lysase',
    5:'Isomerase',
    6:'Ligase',
    7:'Translocase'
}

# Dictionary of nodes for filtering protein classification
def get_other_subclasses(overview_df):
    other_subclasses = np.unique(overview_df.loc[(overview_df['Protein Family'] == 'Other') & ~(overview_df['Protein Subfamily'].apply(lambda x: None in x)), 'Protein Subfamily'].sum())
    other_subclasses = [{"label":x, "value":x} for x in other_subclasses if x not in ['Oxidoreductase','Transferase','Hydrolase','Lyase','Isomerase','Ligase']]
    return other_subclasses

def get_protein_classification_nodes(overview_df):
    return [
        {
            'label':'Enzyme',
            'value':'Enzyme',
            'children':[
                {"label": "Oxidoreductase", "value": "Oxidoreductase"},
                {"label": "Transferase", "value": "Transferase"},
                {"label": "Hydrolase", "value": "Hydrolase"},
                {"label": "Lyase", "value": "Lyase"},
                {"label": "Isomerase", "value": "Isomerase"},
                {"label": "Ligase", "value": "Ligase"},
            ]
        },
        {'label':'GPCR','value':'GPCR'},
        {'label':'Ion Channel','value':'Ion Channel'},
        {
            'label':'Other',
            'value':'Other',
            'children':get_other_subclasses(overview_df)
        },
    ]

# Dictionary containing details of what each omic result means
omic_reassignment = {
    'risk_gwas':'Risk GWAS',
    'SMR_omicsynth':'SMR OmicSynth',
    'coloc_gtex':'Coloc GTEx',
    'tissue_specificity_score':'GTEx Tissue Specificity',
    'SMR_UKB':'SMR Plasma',
    'coloc_UKB':'Coloc Plasma',
    'SMR_Brain':'SMR Brain',
    'coloc_brain':'Coloc Brain',
    'SMR_CSF':'SMR CSF',
    'coloc_CSF':'Coloc CSF',
    'sc_diff_exp':'Single Cell'
}

omic_assessment_descriptions = {
    'risk_gwas':'At least one genetic risk loci physically maps to this gene.',
    'SMR_omicsynth':'There is a causal relationship between at least one omic measurement (eQTLs, caQTLs, mQTLs) of this gene and disease risk.',
    'coloc_gtex':'Changes in bulk tissue RNA expression for this gene and disease risk are both associated with the same (cis) genetic variant(s) in at least one tissue.',
    'tissue_specificity_score':'Closer to 1 the more exclusive our findings suggesting differential RNA expression are to disease specific tissues',
    'SMR_UKB':"There is a causal relationship between this gene's plasma protein levels and disease risk.", 
    'coloc_UKB':'Changes in plasma protein levels for this gene and disease risk are both associated with the same genetic variant(s)', 
    'SMR_Brain':"There is a causal relationship between this gene's brain protein levels and disease risk.", 
    'coloc_Brain':'Changes in brain protein levels for this gene and disease risk are both associated with the same genetic variant(s).', 
    'SMR_CSF':"There is a causal relationship between this gene's CSF protein levels and disease risk.",
    'coloc_CSF':'Changes in CSF protein levels for this gene and disease risk are both associated with the same genetic variant(s).',
    'sc_diff_exp':'RNA for this gene is differentially expressed between patients and controls in at least one single cell type from the human primary motor cortex.'
}

