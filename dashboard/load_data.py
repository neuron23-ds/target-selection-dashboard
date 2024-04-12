import streamlit as st
import pandas as pd
import numpy as np
from ast import literal_eval

from st_files_connection import FilesConnection
conn = st.connection('gcs', type=FilesConnection)

genes = conn.read('target-selection-pipeline/db/genes.csv')

risk_score = conn.read('target-selection-pipeline/db/genes_omic_risk_score.csv').rename({'available':'score_risk_n_sources'}, axis=1)
progression_score = conn.read('target-selection-pipeline/db/genes_omic_progression_score.csv')
druggability_score = conn.read('target-selection-pipeline/db/genes_druggability_score.csv')

studies = conn.read('target-selection-pipeline/db/omics_studies.csv')
gwas_hits = conn.read('target-selection-pipeline/db/omics_gwas_hits.csv')
smr_omicsynth = conn.read('target-selection-pipeline/db/omics_smr_omicsynth.csv')
smr_pqtl = conn.read('target-selection-pipeline/db/omics_smr_pqtl.csv')
coloc_pqtl = conn.read('target-selection-pipeline/db/omics_coloc_pqtl.csv')
coloc_eqtl_gtex = conn.read('target-selection-pipeline/db/omics_coloc_eqtl_gtex.csv')
coloc_eqtl_metabrain = conn.read('target-selection-pipeline/db/omics_coloc_eqtl_metabrain.csv')
single_cell_expression = conn.read('target-selection-pipeline/db/omics_single_cell_results.csv')
single_cell_proteomics = conn.read('target-selection-pipeline/db/omics_single_cell_results_uniprot.csv')

uniprot_blast = conn.read('target-selection-pipeline/db/uniprot_blast.csv')

uniprot_comments = conn.read('target-selection-pipeline/db/uniprot_comments.csv')
uniprot_cross_references = conn.read('target-selection-pipeline/db/uniprot_cross_references.csv')
uniprot_entries = conn.read('target-selection-pipeline/db/uniprot_entries.csv')
uniprot_keywords = conn.read('target-selection-pipeline/db/uniprot_keywords.csv')

chembl_mechanisms = conn.read('target-selection-pipeline/db/chembl_mechanisms.csv')
chembl_molecules = conn.read('target-selection-pipeline/db/chembl_molecules.csv')
chembl_activities = conn.read('target-selection-pipeline/db/chembl_activities.csv')

hpo_disease = conn.read('target-selection-pipeline/db/disease_hpo_diseases.csv')
hpo_terms = conn.read('target-selection-pipeline/db/disease_hpo_terms.csv')
rare_source = conn.read('target-selection-pipeline/db/disease_rare_source.csv')


column_descriptions = {
    'index':'Rank',
    'symbol':'Gene sybmol',
    'name':'Full name',
    'entrez_id':'Entrez ID',
    'uniprot_id':'UniProt ID',
    'ec':'Enzuyme Commission number',
    'ensembl_id':'Ensembl ID', 
    'score_risk':'Reflects the relative strength of association between a gene and indication. This is an internally determined score based on an analysis designed by the N23 DS team.', 
    'score_risk_n_sources':'The number data sources available for the risk score. Some genes may have data avaiable all sources, while others may have data available from only a few sources. You can take this number into account when interpreting the risk score.', 
    'score_progression':'Omics progression score',
    'Protein class':'Protein class', 
    'Molecular function':'Molecular functions, as annotated by UniProt', 
    'Biological process':'Biological process, as annotated by UniProt',
    '3D structure':'Is there a 3D structure available via PDB?', 
    'Avg BLAST identity':'Average identity to human proteins of top 5 BLAST results', 
    'chembl_id':'ChEMBL ID', 
    'Chemical matter':'Chemical matter associated with this gene',
    'gwas_hit_risk':'GWAS hit for risk', 
    'Protein SMR':'Protein SMR', 
    'Expression Coloc':'Expression Coloc', 
    'Additional SMR':'SMR results from NIH omicsynth browser',
    'sc_exp_risk':'Single cell expression risk',
    'gwas_hit_progression':'indicates whether a published genetic association between this gene and ALS survival or progression has been confirmed in the Answer ALS cohort', 
    'publication':'indicates whether a genetic association between this gene and ALS survival or progression has been published',
    'sc_exp_prog_aals':'cell types that show differential gene expression with fast progressing ALS as compared to slow progressing in Answer ALS based on monthly change in ALSFRS',
    'sc_prot_prog_aals':'cell types that show differential protein with fast progressing ALS as compared to slow progressing in Answer ALS based on monthly change in ALSFRS'
}


def string_list_to_list(x):
    return literal_eval(x.replace("' '", "','").replace('\n',','))


class data_loader():
    def __init__(self, indication, gwas):
        self.indication = indication
        self.gwas = gwas
        self.studies = studies[studies.indication == indication]

        self.smr_pqtl = pd.merge(studies, smr_pqtl, left_on='study_id', right_on='gwas', how='right')
        self.coloc_pqtl = pd.merge(studies, coloc_pqtl, left_on='study_id', right_on='gwas', how='right')
        self.coloc_eqtl_gtex = pd.merge(studies, coloc_eqtl_gtex, left_on='study_id', right_on='gwas', how='right')
        self.coloc_eqtl_metabrain = pd.merge(studies, coloc_eqtl_metabrain, left_on='study_id', right_on='gwas', how='right')
        self.coloc_eqtl = pd.concat([self.coloc_eqtl_gtex, self.coloc_eqtl_metabrain])

        self.single_cell_expression = pd.merge(studies, single_cell_expression, on='study_id', how='right')
        self.smr_omicsynth = smr_omicsynth[smr_omicsynth.indication == self.indication]
        self.single_cell_proteomics = pd.merge(studies, single_cell_proteomics, on='study_id', how='right')
    
    @st.cache_data
    def load_main_df(_self):
        print('Loading main dataframe')
        
        # Initialize main dataframe
        main_df = genes.copy()
        
        # Merge with druggability info
        main_df = pd.merge(main_df, druggability_score, on=['entrez_id', 'name', 'symbol', 'uniprot_id', 'ec', 'ensembl_id'])
        main_df['Molecular function'] = main_df['Molecular function'].map(string_list_to_list, na_action='ignore')
        main_df['Biological process'] = main_df['Biological process'].map(string_list_to_list, na_action='ignore')

        # Merge with scores
        main_df = pd.merge(main_df, risk_score[['entrez_id', 'name', 'symbol', 'uniprot_id', 'ec', 'ensembl_id', 'score','score_risk_n_sources','gwas_hit']], on=['entrez_id', 'name', 'symbol', 'uniprot_id', 'ec', 'ensembl_id'])
        main_df = pd.merge(main_df, progression_score[['entrez_id', 'name', 'symbol', 'uniprot_id', 'ec', 'ensembl_id', 'score','gwas_hit','publication']], on=['entrez_id', 'name', 'symbol', 'uniprot_id', 'ec', 'ensembl_id'], suffixes=('_risk','_progression'))
        
        # Merge with gwas-based results from core analysis
        smr_pqtl_for_main_df = _self.smr_pqtl.groupby('symbol').apply(lambda x: list(x.omic[x.smr_hit == 1])).rename('Protein SMR').reset_index()
        main_df = main_df.merge(smr_pqtl_for_main_df, how='left')
        
        coloc_eqtl_for_main_df = _self.coloc_eqtl.groupby('ensembl_id').apply(lambda x: list(x.tissue[(x.n_coloc > 0)|(x.pp_h4_abf>0.8)])).rename('Expression Coloc').reset_index()
        main_df = main_df.merge(coloc_eqtl_for_main_df, how='left')

        # Merge with supplementary omics data
        smr_omicsynth_for_main_df = _self.smr_omicsynth.groupby('symbol').apply(lambda x: list(x.omic[x.smr_hit == 1])).rename('Additional SMR').reset_index()
        main_df = main_df.merge(smr_omicsynth_for_main_df, how='left')

        coloc_pqtl_for_main_df = _self.coloc_pqtl.groupby('symbol').apply(lambda x: list(x.omic[(x.n_coloc > 0) & (x.n_coloc != x.n_hla)])).rename('Protein Coloc').reset_index()
        main_df = main_df.merge(coloc_pqtl_for_main_df, how='left')

        single_cell_expression_for_main_df = _self.single_cell_expression.groupby(['study_id','ensembl_id']).apply(lambda x: list(x.cell[x.padj <= 0.05])).rename('cells').reset_index()
        single_cell_expression_for_main_df = single_cell_expression_for_main_df.replace({'GSE174332':'sc_exp_risk','aals_progression_expression':'sc_exp_prog_aals'})
        single_cell_expression_for_main_df = single_cell_expression_for_main_df.pivot(index='ensembl_id', columns='study_id', values='cells').reset_index()
        main_df = main_df.merge(single_cell_expression_for_main_df, how='left')

        single_cell_proteomics_for_main_df = _self.single_cell_proteomics.groupby(['study_id','uniprot_id']).apply(lambda x: list(x.cell[x.padj <= 0.05])).rename('cells').reset_index()
        single_cell_proteomics_for_main_df = single_cell_proteomics_for_main_df.replace({'aals_progression_proteomics':'sc_prot_prog_aals'})
        single_cell_proteomics_for_main_df = single_cell_proteomics_for_main_df.pivot(index='uniprot_id', columns='study_id', values='cells').reset_index()
        main_df = main_df.merge(single_cell_proteomics_for_main_df, how='left')

        # Collapse genes with multiple uniprot, ensembl, and ec id's into lists
        new_index_cols = ['entrez_id','name','symbol']
        collapse_cols = ['uniprot_id','ensembl_id','ec']
        new_index = main_df.groupby(new_index_cols, dropna=False)[collapse_cols].agg(lambda x: list(set(x))).reset_index()
        main_df = new_index.merge(main_df.drop(columns=collapse_cols), on=new_index_cols).drop_duplicates(subset=new_index_cols)

        # Sort values by score
        main_df.sort_values(by=['score_risk', 'score_progression'], ascending=[False, False], inplace=True)

        # Configure column order
        main_df = main_df[['symbol','name','entrez_id','uniprot_id', 'ec', 'ensembl_id',
                            'score_risk','score_risk_n_sources','score_progression', 
                            'Protein class', 'Molecular function', 'Biological process',
                            '3D structure', 'Avg BLAST identity', 'chembl_id', 'Chemical matter',
                            'gwas_hit_risk', 'Protein SMR','Additional SMR', 'Protein Coloc', 'Expression Coloc',  
                            'gwas_hit_progression','publication','sc_exp_risk', 'sc_exp_prog_aals', 'sc_prot_prog_aals']]
        
        main_df.reset_index(drop=True, inplace=True)
        main_df.index = main_df.index+1
        main_df.reset_index(inplace=True) 
        
        print(main_df)
        return main_df

    # === Omics

    def get_coloc_pqtl_results(self, symbol):
        result_coloc_pqtl = self.coloc_pqtl.loc[(self.coloc_pqtl.symbol == symbol) & (self.coloc_pqtl.n_coloc > 0), ['omic', 'n_coloc','cis_trans']]
        result_coloc_pqtl.replace({'cis_trans':{1:'cis',2:'trans',3:'cis & trans'}}, inplace=True)
        return result_coloc_pqtl
    
    def get_coloc_eqtl_tissues(self, ensembl_id):
        coloc_eqtl_hits = self.coloc_eqtl[(self.coloc_eqtl.n_coloc > 0 )| (self.coloc_eqtl.pp_h4_abf > 0.8)]
        return coloc_eqtl_hits.loc[coloc_eqtl_hits.ensembl_id==ensembl_id, ['omic','tissue']].drop_duplicates()

    def get_smr_results(self, symbol):
        internal_result = smr_pqtl.loc[(smr_pqtl.gwas==self.gwas) & (smr_pqtl.symbol==symbol), ['omic','b_smr','p_smr']]
        omicsynth_result = smr_omicsynth.loc[(smr_omicsynth.indication==self.indication) & (smr_omicsynth.symbol==symbol), ['omic','b_smr','p_smr']]
        smr_result = pd.concat([internal_result,omicsynth_result])
        smr_result['direction'] = np.select([(smr_result.b_smr < 0) & (smr_result.p_smr <= 0.05), 
                                            (smr_result.b_smr > 0) & (smr_result.p_smr <= 0.05)], ['down', 'up'], default='no change')
        return smr_result

    def get_single_cell_diffex(self, ensembl_id, uniprot_id):
        sigle_cell_expression = self.single_cell_expression.loc[(self.single_cell_expression['ensembl_id'] == ensembl_id)].dropna()
        single_cell_proteomics = self.single_cell_proteomics.loc[(self.single_cell_proteomics['uniprot_id'] == uniprot_id)].dropna()
        
        cols = ['cell','phenotype','study_type','log2_change','pvalue']
        result = pd.concat([sigle_cell_expression[cols], single_cell_proteomics[cols]])
        result['direction'] = np.where(result.log2_change > 0, 'up', 'down')
        return result
    

    # === UniProt

    def get_uniprot_comments(self, uniprot_id, comment_type):
        return uniprot_comments.loc[(uniprot_comments.uniprot_id == uniprot_id) & (uniprot_comments.comment_type == comment_type), 'comment'].values
        
    def get_uniprot_keywords(self, uniprot_id, category):
        return uniprot_keywords.loc[(uniprot_keywords.uniprot_id == uniprot_id) & (uniprot_keywords.category == category), 'name'].values

    def get_uniprot_subcell_location(self, uniprot_id):
        l1 = self.get_uniprot_comments(uniprot_id, 'SUBCELLULAR LOCATION')
        l2 = self.get_uniprot_keywords(uniprot_id, 'Cellular component')
        return set(l1) | set(l2)
    
    def get_uniprot_pdb(self, uniprot_id):
        return uniprot_cross_references.loc[(uniprot_cross_references.uniprot_id == uniprot_id) & (uniprot_cross_references.database == 'PDB'), 'id'].values 

    def get_uniprot_blast(self, uniprot_id):
        return uniprot_blast.loc[uniprot_blast.uniprot_id == uniprot_id]

    # === ChEMBL
    def get_chembl_molecules(self, chembl_id):
        result_chembl_mechanisms = chembl_mechanisms[chembl_mechanisms.target_chembl_id == chembl_id]
        if not result_chembl_mechanisms.empty:
            result_chembl_molecules = chembl_molecules.merge(result_chembl_mechanisms, on='molecule_chembl_id')
            return result_chembl_molecules.drop(['parent_molecule_chembl_id','target_chembl_id'], axis=1)

    def get_chembl_activities(self, result_chembl_molecules):
        if isinstance(result_chembl_molecules, pd.DataFrame):
            return chembl_activities.merge(result_chembl_molecules, on='molecule_chembl_id')[['molecule_chembl_id','assay_chembl_id', 'assay_description', 'type','units', 'value']]

    # === Disease
    def get_hpo_disease(self, entrez_id):
        return hpo_disease[hpo_disease['entrez_id'] == entrez_id]
    
    def get_hpo_terms(self, entrez_id):
        return hpo_terms[hpo_terms['entrez_id'] == entrez_id]