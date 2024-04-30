import pandas as pd
import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder

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

# Select an indication
indications = ['ALS'] # PSP, FTD
with st.sidebar:
   indication = st.selectbox("Select an indication", indications)


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
builder.configure_pagination(paginationPageSize=10000, paginationAutoPageSize=False)
builder.configure_side_bar()
builder.configure_grid_options(tooltipInteraction=True, tooltipShowDelay=1000, alwaysShowHorizontalScroll=True)

from load_data import column_descriptions
for col, desc in column_descriptions.items():
    builder.configure_column(col, headerTooltip=desc)

builder.configure_default_column(filterable=True, sorteable=True, groupable=False, editable=False,
                                 enableCellTextSelection=True, ensureDomOrder=True, enableValue=True, 
                                 hide=True, filter='agSetColumnFilter')

builder.configure_column('index', hide=False, headerTooltip='index')
builder.configure_column('symbol', hide=False)
builder.configure_column('name', hide=False)
builder.configure_column('score_risk', hide=False, type=['numericColumn','numberColumnFilter','customNumericFormat'], precision=2)
builder.configure_column('score_progression', hide=False, type=['numericColumn','numberColumnFilter','customNumericFormat'], precision=0)
builder.configure_column('Protein class', hide=False, filter='agSetColumnFilter')


grid_options = builder.build()
grid_return = AgGrid(main_df, 
                     height=500, 
                     gridOptions=grid_options,
                     allow_unsafe_jscode=True,
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


st.write(f"## Function")

description = data.get_uniprot_comments(uniprot_id, 'FUNCTION')
if any(description):
   st.write(description[0])
else: 
   st.write("*No description available*")

molecular_function = data.get_uniprot_keywords(uniprot_id, 'Molecular function')
if any(molecular_function):
    st.write(f"**Molecular function**: {', '.join(molecular_function)}")

biological_proceses = data.get_uniprot_keywords(uniprot_id, 'Biological process')
if any(biological_proceses):
    st.write(f"**Biological process**: {', '.join(biological_proceses)}")


st.markdown('***')


st.write(f"## Expression")

st.write(f"### Tissue specific")

st.write(f"**From internal {indication}-specific colocalization analysis**:")
coloc_pqtl_tissues = data.get_coloc_pqtl_results(symbol)
if len(coloc_pqtl_tissues) == 0:
    st.write(f'*No protein data available*')
elif any(coloc_pqtl_tissues):
    st.dataframe(coloc_pqtl_tissues, hide_index=True)
else:
    st.write('No protein hits')


coloc_eqtl_tissues = data.get_coloc_eqtl_tissues(ensembl_id)
if len(coloc_eqtl_tissues) == 0:
    st.write(f'*No gene expression data available*')
elif any(coloc_eqtl_tissues):
    st.dataframe(coloc_eqtl_tissues, hide_index=True)
else:
    st.write('No gene expression hits')

st.text("")

st.write(f"**From internal {indication}-specific SMR analysis**:")
get_smr_results = data.get_smr_results(symbol)
if not get_smr_results.empty:
    st.dataframe(get_smr_results, hide_index=True)
else:
    st.write("*No data available*")

st.text("")

st.write("**From UniProt**:")
uniprot_tissues = data.get_uniprot_comments(uniprot_id, 'TISSUE SPECIFICITY')
if any(uniprot_tissues):
    st.write(f"{uniprot_tissues[0]}")
else:
    st.write("*No hits or no data available*")

st.text("")

st.write(f"### Cell specific")
single_cell_diffex = data.get_single_cell_diffex(ensembl_id, uniprot_id)
if not single_cell_diffex.empty:
    st.dataframe(single_cell_diffex, hide_index=True)
else:
    st.write("*No data available*")

st.text("")

st.write(f"### Subcellular location")
uniprot_subcell_location = data.get_uniprot_subcell_location(uniprot_id)
if any(uniprot_subcell_location):
    st.write(f"**From UniProt**: {', '.join(uniprot_subcell_location)}")


st.markdown('***')


st.write(f"## Network & Pathway")

uniprot_pathways = data.get_uniprot_comments(uniprot_id, 'PATHWAY')
if any(uniprot_pathways):
    st.write(f"**Pathways via uniprot**: {''.join(uniprot_pathways)}")

uniprot_interactors = data.get_uniprot_comments(uniprot_id, 'INTERACTION')
if any(uniprot_interactors):
    st.write(f"**Interactions via uniprot**: {''.join(uniprot_interactors)}")


st.markdown('***')


st.write(f"## Druggability")

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


st.write(f"## Disease Involvement")

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


