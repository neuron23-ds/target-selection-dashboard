import streamlit as st
from streamlit_utils import n23_icon, add_logo, check_password

st.set_page_config(
  page_title='Target Selection Dashboard',
  page_icon=n23_icon,
  # layout='wide'
  )

add_logo()

if not check_password():
    st.stop()  # Do not continue if check_password is not True.

st.markdown("# Target Selection Dashboard")

st.markdown("## How to use this dashboard")
with st.expander('**Overview**', expanded=False):
  st.markdown("""
  Select an indication from the sidebar to load results for that indication.
  At the top, you will find the main table, which contains information for all targets for the selected indication.
  Targets in the main table have been included based a risk omics analysis and progression omics analyses as performed internally by the Data Science team.
  See [this slide deck](https://neuron23com.sharepoint.com/:p:/s/Neuron23Inc/EUt6Ioiz5BtKtc2i_bqCDS8BtGYSnRwm6VAllE73hqCppw?e=MXgjto) for more details on the omics anlyses
  """) 
with st.expander('**Using the main table**'):
  st.markdown("""
  The main table interactive, allowing you to sort, filter, and search columns. You can also select a row to view detailed information about that target.
  
  **Selecting a row:**
  Click on a row to view detailed information about that target.
  
  **Sorting:**
  Click on a column header to sort Click again to reverse the sort order. Click again to remove the sort.\n
  
  **Filtering:**
  There are two way to filter: use the "Filters" tab on the right side of table or click three lines (≡) in a column header.            

  **Displaying columns:**
  Customize columns to display by clicking the "Columns" tab on the right side of the table. 
  There are 30 columns available to display in the main table. 
  The column presets are a good place to start if you're not sure which columns to display.
  
  **Exporting data:**
  Export your customized main table to a CSV file clicking the "Export" button at the bottom of the table.
              
  **Column descritpions:**
  Hover over a column header to view the description of that column. The same descriptions are provided below:
  """)
  from load_data import column_descriptions
  import pandas as pd
  df = pd.DataFrame.from_dict(column_descriptions, orient='index').reset_index().rename(columns={'index':'Column', 0:'Description'})
  st.dataframe(df, hide_index=True, use_container_width=True)

with st.expander("**Detailed information**"):
  st.markdown("""
      You can click on a row in the main table to view detailed information about that target.        
      It is mostly from the main table, but with more details.
      
      The detailed has been divided into 5 sections:
      - **Function**: General information about the function of the target. Sourced from UniProt, PANTHER, GO and Reactome.
      - **Expression**: Three levels of expression data.    
        - **Tissue specific**: Internal coloc and SMR results. Uniprot information.
        - **Cell specific**: Internal single cell expression results. Uniprot information.
        - **Subcellular location**: Uniprot information.
      - **Network & Pathway**: 
      - **Druggability**:
      - **Disease involvement**: 
  """)

st.markdown("## Data Sources") 
with st.expander("**GWAS Summary Statistics**", expanded=False):
  st.markdown("* ALS Risk: [Van Rheenen 2022](https://doi.org/10.1038/s41588-022-01020-3)")

with st.expander("**Omics Resources**", expanded=False):
  st.markdown("""
              
    * UK Biobank Plasma Proteomics
    * Yang Brain and CSF pQTL
    * MetaBrain eQTL
    * GTEX eQTL
    * NIH CARDD Omicsynth Browser
              
  """)
