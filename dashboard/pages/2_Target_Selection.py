import pandas as pd
import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder, ColumnsAutoSizeMode

from streamlit_utils import n23_icon, add_logo, check_password


st.set_page_config(
  page_title='Target Selection Dashboard',
  page_icon=n23_icon,
  layout='wide'
 )

add_logo()

if not check_password():
    st.stop()  # Do not continue if check_password is not True.


# ===========
#   Sidebar 
# ===========

with st.sidebar:
    # Select an indication
    indications = ['ALS'] # PSP, FTD
    indication = st.selectbox("Select an indication", indications)

    # Select column presets
    column_presets = ['Default','Omics Analysis Risk','Omics Analysis Progression','Biology','Chemistry']
    column_preset = st.radio('Select a column preset', column_presets, index=0)


### TEMPORARY GWAS, need to come up with a way to do GWAS
if indication == 'ALS':
   gwas = 'ALS_GWAS_vanRheenen' 

# Load data
from load_data import data_loader
data = data_loader(indication, gwas)

# ======================
#   Overview Dataframe
# ======================

st.title(f'Results for {indication}')

# Main dataframe
main_df = data.load_main_df()

builder = GridOptionsBuilder.from_dataframe(main_df)

builder.configure_selection(selection_mode='single', pre_selected_rows=[0], use_checkbox=False, suppressRowDeselection=True)
builder.configure_pagination(enabled=True, paginationPageSize=10000, paginationAutoPageSize=False)
builder.configure_side_bar()
builder.configure_grid_options(tooltipInteraction=True, tooltipShowDelay=1000, alwaysShowHorizontalScroll=True)

from load_data import column_descriptions
for col, desc in column_descriptions.items():
    builder.configure_column(col, headerTooltip=desc)

builder.configure_default_column(filterable=True, sorteable=True, groupable=False, editable=False,
                                 enableCellTextSelection=True, ensureDomOrder=True, enableValue=True, 
                                 hide=True, filter='agSetColumnFilter')

builder.configure_column('symbol', hide=False)
builder.configure_column('name', hide=False)
builder.configure_column('score_risk', hide=False, type=['numericColumn','numberColumnFilter','customNumericFormat'], precision=2)
builder.configure_column('score_progression', hide=False, type=['numericColumn','numberColumnFilter','customNumericFormat'], precision=0)

if column_preset == 'Default':
    cols_to_show = []
elif column_preset == 'Omics Analysis Risk':
    cols_to_show = ['score_risk','score_risk_n_sources', 'gwas_hit_risk', 'protein_smr_risk','additional_smr_risk', 'protein_coloc_risk', 'expression_coloc_risk','single_cell_expression_risk']
elif column_preset == 'Omics Analysis Progression':
    cols_to_show = ['score_progression', 'gwas_hit_progression','publication_progression', 'single_cell_expression_progression', 'single_cell_protein_progression']
elif column_preset == 'Biology':
    cols_to_show = ['protein_class','molecular_function','biological_process','pathway','cellular_component','uniprot_id']
elif column_preset == 'Chemistry':
    cols_to_show = ['Protein class','3D structure', 'Avg BLAST identity', 'Chemical matter','chembl_id']

for col in cols_to_show:
    builder.configure_column(col, hide=False)


grid_options = builder.build()
grid_return = AgGrid(main_df, 
                     height=500, 
                     gridOptions=grid_options,
                     allow_unsafe_jscode=True,
                     fit_columns_on_grid_load=True,
                     update_on=['sortChanged', 'filterChanged', 'filterModified','columnMoved','columnVisible'])


# Add download button
def convert_df(main_df, grid_return):
    visible_cols = [x['colId'] for x in grid_return['columns_state'] if not x['hide']]
    col_order = grid_return['rows_id_after_sort_and_filter']
    df_for_download = main_df.iloc[col_order][visible_cols]
    return df_for_download.to_csv(index=False).encode('utf-8')
try:
    st.download_button(label='Download main table', data=convert_df(main_df, grid_return), file_name='main_table.csv', mime='text/csv')
except TypeError:
    pass

# Jump to sections
st.markdown("Jump to: [Function](#function), [Expression & Localization](#expression-localization), [Network & Pathway](#network-pathway), [Druggability](#druggability), [Disease Involvement](#disease-involvement)") 

st.markdown('***')

# =========================
#   Detailed Information
# =========================

try:
    selected = grid_return['selected_rows'].iloc[0]
except AttributeError:
    st.write('Please select a row from the table to view detailed information.')
    st.stop()

name = selected['name']
symbol = selected['symbol']
uniprot_id = selected['uniprot_id'][0] # Only displays first uniprot_id, need to come up with a way to display all
ensembl_id = selected['ensembl_id'][0]
entrez_id = int(selected['entrez_id'])
chembl_id = selected['chembl_id']


# st.write(f"# {symbol}")

st.write(f"# {symbol} ({name})")


with st.container():
    st.write(f"## FUNCTION")

    description = data.get_uniprot_comments(uniprot_id, 'FUNCTION')
    if any(description):
        st.write("**Description:**", description[0])
    else: 
        st.write("**Description:** *No description available*")

    st.write("**Molecular functions**:")
    molecular_functions = data.get_molecular_functions(uniprot_id)
    st.dataframe(molecular_functions, hide_index=True)

    st.write("**Biological processes**:")
    biological_processes = data.get_biological_processes(uniprot_id)
    st.dataframe(biological_processes, hide_index=True)


st.markdown('***')


with st.container():
    st.write(f"## EXPRESSION & LOCALIZATION")

    with st.container(border=True):
        st.write(f"### {indication}-Specific")
        
        st.write("#### Tissue-level")
        
        st.write("**SMR Results**")
        st.write(f"Tests whether the effect size of a SNP on {indication} is mediated by {symbol} protein expression.")
        st.dataframe(data.get_smr_results(symbol), hide_index=True)
        
        st.write("**Colocalization Results**")
        st.write(f"Identifies genes that contain a SNP which introduces phenotypic change in both {symbol} protein expression and {indication} risk.")
        st.dataframe(data.get_coloc_pqtl_results(symbol), hide_index=True)
        st.write(f"Identifies genes that contain a SNP which introduces phenotypic change in both {symbol} RNA expression and {indication} risk.")
        st.dataframe(data.get_coloc_eqtl_tissues(ensembl_id), hide_index=True)


        st.write("#### Cell-level")
        st.write("**Single Cell Differential Expression Results**")
        st.write(f"Cells below show differential {symbol} expression in  RNA (sc_dge) or protein (sc_dpe) with respect to the {indication} phenotype indicated.")
        st.dataframe(data.get_single_cell_diffex(ensembl_id, uniprot_id), hide_index=True)

    with st.container():
        st.write(f"### General (non indication-specific)")
        st.write("#### Tissue-level")
        uniprot_tissues = data.get_uniprot_comments(uniprot_id, 'TISSUE SPECIFICITY')
        if any(uniprot_tissues):
            st.write(f"**From UniProt**: {uniprot_tissues[0]}")
        else:
            st.write("*No hits or no data available*")

        st.write("#### Subcellular Location")
        uniprot_subcell_location = data.get_uniprot_subcell_location(uniprot_id)
        if any(uniprot_subcell_location):
            st.write(f"**From UniProt**: {', '.join(uniprot_subcell_location)}")


st.markdown('***')


with st.container():
    st.write(f"## NETWORK & PATHWAY")
    st.write(f"{symbol} is involved in the following pathways:")
    st.dataframe(data.get_pathways(uniprot_id), hide_index=True)

    uniprot_interactors = data.get_uniprot_comments(uniprot_id, 'INTERACTION')
    if any(uniprot_interactors):
        st.write(f"**Interactions via uniprot**: {''.join(uniprot_interactors)}")


st.markdown('***')


with st.container():
    st.write(f"## DRUGGABILITY")

    regulation = data.get_uniprot_comments(uniprot_id, 'ACTIVITY REGULATION')
    if any(regulation):
        st.write(f"**Activity regulation**: {', '.join(regulation)}")

    st.write(f"### Structure")
    pdb_ids = data.get_uniprot_pdb(uniprot_id)
    if any(pdb_ids):
        st.write(f"**Known strcuture PDB ID(s)**: {', '.join(pdb_ids)}")
    else:
        blast_result = data.get_uniprot_blast(uniprot_id)
        if not blast_result.empty:
            st.dataframe(blast_result, hide_index=True)

    st.write(f"### Chemical Matter")
    molecules = data.get_chembl_molecules(chembl_id)
    if isinstance(molecules, pd.DataFrame):
        st.write(f"**Molecules targeting this protein**:")
        st.dataframe(molecules, hide_index=True)

    activities = data.get_chembl_activities(molecules)
    if isinstance(activities, pd.DataFrame):
        if not activities.empty:
            st.write(f"**Activities of molecules targeting this protein**:")
            st.dataframe(activities, hide_index=True)


st.markdown('***')


with st.container():
    st.write(f"## DISEASE INVOLVEMENT")

    uniprot_disease = data.get_uniprot_keywords(uniprot_id, 'Disease')
    if any(uniprot_disease):
        st.write(f"**Disease(s) associated with this gene via UniProt**: {', '.join(uniprot_disease)}")

    hpo_disease = data.get_hpo_disease(entrez_id)
    if not hpo_disease.empty:
        st.write('**Disease(s) associated with this gene via HPO**:')
        st.dataframe(hpo_disease, hide_index=True)

    hpo_terms = data.get_hpo_terms(entrez_id)
    if not hpo_terms.empty:
        st.write('**Symptoms associated with this gene via HPO**:')
        st.dataframe(hpo_terms, hide_index=True)    


