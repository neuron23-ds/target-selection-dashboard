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

uniprot_entries = conn.read('target-selection-pipeline/db/uniprot_entries.csv')
uniprot_cross_references = conn.read('target-selection-pipeline/db/uniprot_cross_references.csv')
uniprot_comments = conn.read('target-selection-pipeline/db/uniprot_comments.csv')
uniprot_keywords = conn.read('target-selection-pipeline/db/uniprot_keywords.csv')
uniprot_blast = conn.read('target-selection-pipeline/db/uniprot_blast.csv')

panther_gene_info = conn.read('target-selection-pipeline/db/panther_gene_info.csv')

chembl_mechanisms = conn.read('target-selection-pipeline/db/chembl_mechanisms.csv')
chembl_molecules = conn.read('target-selection-pipeline/db/chembl_molecules.csv')
chembl_activities = conn.read('target-selection-pipeline/db/chembl_activities.csv')

hpo_disease = conn.read('target-selection-pipeline/db/disease_hpo_diseases.csv')
hpo_terms = conn.read('target-selection-pipeline/db/disease_hpo_terms.csv')
rare_source = conn.read('target-selection-pipeline/db/disease_rare_source.csv')

publications = conn.read('target-selection-pipeline/db/omics_publications_genetics.csv')

column_descriptions = {
    'index':'Rank',
    # ID's
    'symbol':'Gene sybmol',
    'name':'Full name',
    'entrez_id':'Entrez ID',
    'uniprot_id':'UniProt ID',
    'ec':'Enzuyme Commission number',
    'ensembl_id':'Ensembl ID', 
    'chembl_id':'ChEMBL ID', 
    
    # Risk score
    'score_risk':'Reflects the relative strength of association between a gene and indication. This is an internally determined score based on an analysis designed by the N23 DS team.', 
    'score_risk_n_sources':'The number data sources available for the risk score. Some genes may have data avaiable all sources, while others may have data available from only a few sources. You can take this number into account when interpreting the risk score.', 
    'gwas_hit_risk':'GWAS hit for risk', 
    'protein_smr_risk':'Tissues with positive result from protein SMR analysis. Sources include plasma, brain, and CSF.', 
    'expression_coloc_risk':'Tissues with positive result from gene expression colocalization analysis. Sources include GTEx (47 tissues) and MetaBrain (5 brain-specific tissues).', 
    'protein_coloc_risk':'Tiissues with positive result from protein expression colocalization analysis. Sources include plasma, brain, and CSF.', 
    'additional_smr_risk':'SMR results from NIH omicsynth browser',
    'single_cell_expression_risk':'Cell types that show differential gene expression risk between diseased and control subjects.',
    
    # Progression sore
    'score_progression':'Omics progression score',
    'single_cell_expression_progression':'Cell types that show differential gene expression with fast progressing ALS as compared to slow progressing in Answer ALS based on monthly change in ALSFRS',
    'single_cell_protein_progression':'Cell types that show differential protein with fast progressing ALS as compared to slow progressing in Answer ALS based on monthly change in ALSFRS',
    'gwas_hit_progression':'Indicates whether a published genetic association between this gene and ALS survival or progression has been confirmed in the Answer ALS cohort', 
    'publication_progression':'Indicates whether a genetic association between this gene and ALS survival or progression has been published',
    
    # Biological terms
    'protein_class':'Protein classifications of this protein sourced from PANTHER. See About page for more info on the sources of this data.',
    'molecular_function':'Molecular functions of this protein. Sourced from UniProt and GO See About page for more info on the sources of this data.',
    'biological_process':'Biological processes of this protein. Sourced from UniProt and GO. See About page for more info on the sources of this data.',
    'pathway':'Pathways this protein is involved with. Sourced from PANTHER and Reactome. See About page for more info on the sources of this data.',
    'cellular_component':'Cellular components of this protein. Sourced from UniProt and GO. See About page for more info on the sources of this data.',
    
    # Chemistry info
    'Protein class':'Protein classificaion as built internally. Manually classifies proteins into Enzyme, GPCR, Ion channel, Receptor, or Other.', 
    '3D structure':'Is there a 3D structure available via PDB?', 
    'Avg BLAST identity':'Average identity to human proteins of top 5 BLAST results', 
    'Chemical matter':'Chemical matter associated with this gene',
}

panther_annotation_dataset_meaning = {
    'PANTHER GO-Slim Molecular Function':'molecular_function',
    'PANTHER GO-Slim Biological Process':'biological_process',
    'PANTHER GO-Slim Cellular Component':'cellular_component',
    'PANTHER Protein Class':'protein_class',
    'PANTHER Pathways':'pathway',
    'GO molecular function complete':'molecular_function',
    'GO biological process complete':'biological_process',
    'GO cellular component complete':'cellular_component',
    'Reactome pathways':'pathway'
}

uniprot_keyword_category_meaning = {
    'Molecular function':'molecular_function',
    'Biological process':'biological_process',
    'Cellular component':'cellular_component',
}

class data_loader():
    def __init__(self, indication, gwas):
        self.indication = indication
        self.gwas = gwas
        self.studies = studies[studies.indication == indication]

        self.gwas_hits = pd.merge(studies, gwas_hits, left_on='study_id', right_on='gwas', how='right')
        self.publications = publications[publications.indication == indication]

        self.smr_pqtl = pd.merge(self.studies, smr_pqtl, left_on='study_id', right_on='gwas', how='right')
        self.coloc_pqtl = pd.merge(self.studies, coloc_pqtl, left_on='study_id', right_on='gwas', how='right')
        self.coloc_eqtl_gtex = pd.merge(self.studies, coloc_eqtl_gtex, left_on='study_id', right_on='gwas', how='right')
        self.coloc_eqtl_metabrain = pd.merge(self.studies, coloc_eqtl_metabrain, left_on='study_id', right_on='gwas', how='right')
        self.coloc_eqtl = pd.concat([self.coloc_eqtl_gtex, self.coloc_eqtl_metabrain])
        
        self.smr_omicsynth = smr_omicsynth[smr_omicsynth.indication == self.indication]
        self.single_cell_expression_risk = pd.merge(self.studies[self.studies.phenotype == 'risk'], single_cell_expression, on='study_id', how='inner')
        self.single_cell_expression_prog = pd.merge(self.studies[self.studies.phenotype == 'progression'], single_cell_expression, on='study_id', how='inner')
        
        self.single_cell_proteomics = pd.merge(self.studies, single_cell_proteomics, on='study_id', how='right')
    
    @st.cache_data
    def load_main_df(_self):
        print('Loading main dataframe')
        
        # Initialize main dataframe
        main_df = genes.copy()
        
        # Merge with druggability info
        druggability_score_for_main_df = druggability_score[['entrez_id', 'name', 'symbol', 'uniprot_id', 'ec', 'ensembl_id', 'Protein class','3D structure', 'Avg BLAST identity', 'chembl_id', 'Chemical matter']]
        main_df = pd.merge(main_df, druggability_score_for_main_df, on=['entrez_id', 'name', 'symbol', 'uniprot_id', 'ec', 'ensembl_id'])

        # Merge with scores
        main_df = pd.merge(main_df, risk_score[['entrez_id', 'name', 'symbol', 'uniprot_id', 'ec', 'ensembl_id', 'score','score_risk_n_sources']], on=['entrez_id', 'name', 'symbol', 'uniprot_id', 'ec', 'ensembl_id'])
        main_df = pd.merge(main_df, progression_score[['entrez_id', 'name', 'symbol', 'uniprot_id', 'ec', 'ensembl_id', 'score']], on=['entrez_id', 'name', 'symbol', 'uniprot_id', 'ec', 'ensembl_id'], suffixes=('_risk','_progression'))
        
        # Merge with gwas-based results from core analysis
        gwas_hits_risk_for_main_df = _self.gwas_hits[_self.gwas_hits.phenotype == 'risk'].groupby('symbol').snp.apply(list).rename('gwas_hit_risk').reset_index()
        main_df = main_df.merge(gwas_hits_risk_for_main_df, how='left', on='symbol')

        smr_pqtl_for_main_df = _self.smr_pqtl.groupby('symbol').apply(lambda x: list(x.omic[x.smr_hit == 1])).rename('protein_smr_risk').reset_index()
        main_df = main_df.merge(smr_pqtl_for_main_df, how='left')
        
        coloc_eqtl_for_main_df = _self.coloc_eqtl.groupby('ensembl_id').apply(lambda x: list(x.tissue[(x.n_coloc > 0)|(x.pp_h4_abf>0.8)])).rename('expression_coloc_risk').reset_index()
        main_df = main_df.merge(coloc_eqtl_for_main_df, how='left')

        # Merge with supplementary omics data
        smr_omicsynth_for_main_df = _self.smr_omicsynth.groupby('symbol').apply(lambda x: list(x.omic[x.smr_hit == 1])).rename('additional_smr_risk').reset_index()
        main_df = main_df.merge(smr_omicsynth_for_main_df, how='left')

        coloc_pqtl_for_main_df = _self.coloc_pqtl.groupby('symbol').apply(lambda x: list(x.omic[(x.n_coloc > 0) & (x.n_coloc != x.n_hla)])).rename('protein_coloc_risk').reset_index()
        main_df = main_df.merge(coloc_pqtl_for_main_df, how='left')

        single_cell_expression_risk_for_main_df = _self.single_cell_expression_risk.groupby('ensembl_id').apply(lambda x: list(x.cell[x.padj <= 0.05])).rename('single_cell_expression_risk').reset_index()
        main_df = main_df.merge(single_cell_expression_risk_for_main_df, how='left')

        single_cell_expression_prog_for_main_df = _self.single_cell_expression_prog.groupby('ensembl_id').apply(lambda x: list(x.cell[x.pvalue <= 0.05])).rename('single_cell_expression_progression').reset_index()
        main_df = main_df.merge(single_cell_expression_prog_for_main_df, how='left')

        single_cell_proteomics_for_main_df = _self.single_cell_proteomics.groupby(['study_id','uniprot_id']).apply(lambda x: list(x.cell[x.pvalue <= 0.05])).rename('cells').reset_index()
        single_cell_proteomics_for_main_df = single_cell_proteomics_for_main_df.replace({'aals_progression_proteomics':'single_cell_protein_progression'})
        single_cell_proteomics_for_main_df = single_cell_proteomics_for_main_df.pivot(index='uniprot_id', columns='study_id', values='cells').reset_index()
        main_df = main_df.merge(single_cell_proteomics_for_main_df, how='left')

        # Temp AALS progression stuff
        gwas_hits_progression_for_main_df = _self.gwas_hits[(_self.gwas_hits.phenotype == 'progression') & (_self.gwas_hits.pvalue <= 0.05)].groupby('symbol').snp.apply(list).rename('gwas_hit_progression').reset_index()
        main_df = main_df.merge(gwas_hits_progression_for_main_df, how='left', on='symbol')

        pubications_progression_for_main_df = _self.publications[_self.publications.phenotype == 'Answer ALS progression'].groupby('symbol').reference.apply(list).rename('publication_progression').reset_index()
        main_df = main_df.merge(pubications_progression_for_main_df, how='left', on='symbol')

        # Merge with biological terms
        panther_gene_info_for_main_df = panther_gene_info.copy()
        panther_gene_info_for_main_df.replace({'annotation_dataset':panther_annotation_dataset_meaning}, inplace=True)
        panther_gene_info_for_main_df.rename({'annotation_dataset':'header','annotation_name':'term'}, axis=1, inplace=True)

        uniprot_keywords_for_main_df = uniprot_keywords.copy()
        uniprot_keywords_for_main_df = uniprot_keywords_for_main_df[uniprot_keywords_for_main_df.category.isin(uniprot_keyword_category_meaning.keys())]
        uniprot_keywords_for_main_df.replace({'category':uniprot_keyword_category_meaning}, inplace=True)
        uniprot_keywords_for_main_df.rename({'category':'header','name':'term'}, axis=1, inplace=True)

        terms = pd.concat([panther_gene_info_for_main_df, uniprot_keywords_for_main_df])[['uniprot_id','header','term']].drop_duplicates()
        terms = terms.groupby(['uniprot_id','header']).term.apply(list).reset_index()
        terms = terms.pivot(index='uniprot_id', columns='header', values='term').reset_index()
        main_df = main_df.merge(terms, how='left', on='uniprot_id')

        # Collapse genes with multiple uniprot, ensembl, and ec id's into lists
        new_index_cols = ['entrez_id','name','symbol']
        collapse_cols = ['uniprot_id','ensembl_id','ec']
        new_index = main_df.groupby(new_index_cols, dropna=False)[collapse_cols].agg(lambda x: list(set(x))).reset_index()
        main_df = new_index.merge(main_df.drop(columns=collapse_cols), on=new_index_cols).drop_duplicates(subset=new_index_cols)

        # Sort values by score
        main_df.sort_values(by=['score_risk', 'score_progression'], ascending=[False, False], inplace=True)

        # Configure column order
        main_df = main_df[[
            'symbol','name','entrez_id','uniprot_id', 'ec', 'ensembl_id', 'chembl_id', # ids
            'score_risk','score_risk_n_sources', 'gwas_hit_risk', 'protein_smr_risk','additional_smr_risk', 'protein_coloc_risk', 'expression_coloc_risk','single_cell_expression_risk', # risk score
            'score_progression', 'gwas_hit_progression','publication_progression', 'single_cell_expression_progression', 'single_cell_protein_progression', # progression score
            'protein_class','molecular_function','biological_process','pathway','cellular_component',# biological terms
            'Protein class','3D structure', 'Avg BLAST identity', 'Chemical matter', # chemistry info
        ]]
        
        main_df.reset_index(drop=True, inplace=True)
        main_df.index = main_df.index+1
        main_df.reset_index(inplace=True) 
    
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
        single_cell_expression_risk = self.single_cell_expression_risk.loc[(self.single_cell_expression_risk['ensembl_id'] == ensembl_id) & (self.single_cell_expression_risk.padj <= 0.05)].dropna()
        sigle_cell_expression_prog = self.single_cell_expression_prog.loc[(self.single_cell_expression_prog['ensembl_id'] == ensembl_id)].dropna()
        single_cell_proteomics = self.single_cell_proteomics.loc[(self.single_cell_proteomics['uniprot_id'] == uniprot_id)].dropna()
        
        cols = ['cell','phenotype','study_type','log2_change','pvalue']
        result = pd.concat([single_cell_expression_risk[cols], sigle_cell_expression_prog[cols], single_cell_proteomics[cols]])
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

    # === Function
    def get_molecular_functions(self, uniprot_id):
        mf1 = uniprot_keywords.loc[(uniprot_keywords.uniprot_id == uniprot_id) & (uniprot_keywords.category == "Molecular function"), ['name', 'category']].replace("Molecular function", "Uniprot molecular function keywords").rename({'name':'Molecular Function', 'category':'Annotation Source'}, axis=1)
        mf2 = panther_gene_info.loc[(panther_gene_info.uniprot_id == uniprot_id) & (panther_gene_info.annotation_dataset == "PANTHER GO-Slim Molecular Function"), ['annotation_name', 'annotation_dataset']].rename({'annotation_name':'Molecular Function', 'annotation_dataset':'Annotation Source'}, axis=1)
        mf3 = panther_gene_info.loc[(panther_gene_info.uniprot_id == uniprot_id) & (panther_gene_info.annotation_dataset == "GO molecular function complete"), ['annotation_name', 'annotation_dataset']].rename({'annotation_name':'Molecular Function', 'annotation_dataset':'Annotation Source'}, axis=1)
        mf = pd.concat([mf1, mf2, mf3]).drop_duplicates(subset='Molecular Function', keep='first')
        return mf

    def get_biological_processes(self, uniprot_id):
        bp1 = uniprot_keywords.loc[(uniprot_keywords.uniprot_id == uniprot_id) & (uniprot_keywords.category == "Biological process"), ['name', 'category']].replace("Biological process", "Uniprot biological process keywords").rename({'name':'Biological process', 'category':'Annotation Source'}, axis=1)
        bp2 = panther_gene_info.loc[(panther_gene_info.uniprot_id == uniprot_id) & (panther_gene_info.annotation_dataset == "PANTHER GO-Slim Biological Process"), ['annotation_name', 'annotation_dataset']].rename({'annotation_name':'Biological process', 'annotation_dataset':'Annotation Source'}, axis=1)
        bp3 = panther_gene_info.loc[(panther_gene_info.uniprot_id == uniprot_id) & (panther_gene_info.annotation_dataset == "GO biological process complete"), ['annotation_name', 'annotation_dataset']].rename({'annotation_name':'Biological process', 'annotation_dataset':'Annotation Source'}, axis=1)
        bp = pd.concat([bp1, bp2, bp3]).drop_duplicates(subset='Biological process', keep='first')
        return bp

    # === Network and Pathway
    def get_pathways(self, uniprot_id):
        p1 = panther_gene_info.loc[(panther_gene_info.uniprot_id == uniprot_id) & (panther_gene_info.annotation_dataset == "PANTHER Pathways"), ['annotation_name', 'annotation_dataset']].rename({'annotation_name':'Pathway', 'annotation_dataset':'Annotation Source'}, axis=1)
        p2 = panther_gene_info.loc[(panther_gene_info.uniprot_id == uniprot_id) & (panther_gene_info.annotation_dataset == "Reactome pathways"), ['annotation_name', 'annotation_dataset']].rename({'annotation_name':'Pathway', 'annotation_dataset':'Annotation Source'}, axis=1)
        p3 = [x.strip() for x in ';'.join(uniprot_comments.loc[(uniprot_comments.uniprot_id == uniprot_id) & (uniprot_comments.comment_type == 'PATHWAY'), 'comment'].values).split(';')]
        try:
            p3.remove('')
        except:
            pass
        p3 = pd.DataFrame({'Pathway':p3, 'Annotation Source':'UniProt'})
        p = pd.concat([p1, p2, p3]).drop_duplicates(subset='Pathway', keep='first')
        return p


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