import pandas as pd
import os
import argparse
import numpy as np
import time
import logging
import sys

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler("log.txt", mode="a", encoding="utf-8"),
        logging.StreamHandler(sys.stdout)
    ]
)
parser = argparse.ArgumentParser(description="Screening variants based on gene name from an Excel file.")
parser.add_argument("-I", "--input",required=True, type=str, help="Path to the input Excel file.")
parser.add_argument("--gene_name", required=True, nargs="+", type=str, help="Gene name to filter variants.")
parser.add_argument("-O", "--output", required=True, type=str, help="Path to the output Excel file.")
args = parser.parse_args()

INPUT_DIR = args.input
GENE_NAME = args.gene_name
OUTPUT = args.output

genes = []
for g in GENE_NAME:
    genes += g.split()
genes = list(set(genes))
columns = [
    #General
    "CHROM","POS","REF","ALT","DP","AD","QUAL","MQ","Zygosity","FILTER","Effect",
    "Putative_Impact","Gene_Name","Feature_Type","Feature_ID","Transcript_BioType",
    "Rank/Total","HGVS.c","HGVS.p","REF_AA","ALT_AA","cDNA_pos","cDNA_length","CDS_pos",
    "CDS_length","AA_pos","AA_length","Distance",
    #dbSNP138 annotation
    "dbSNP138_ID","dbSNP156_ID",
    #1000 genomes phase 3 annotation
    "p3_1000G_AF","p3_1000G_AFR_AF","p3_1000G_AMR_AF","p3_1000G_EAS_AF","p3_1000G_EUR_AF","p3_1000G_SAS_AF",
    #EVS annotation
    "ESP6500_MAF_EA","ESP6500_MAF_AA","ESP6500_MAF_ALL",
    #Clinvar annotation
    "CLINVAR_CLNSIG","CLINVAR_CLNDISDB","CLINVAR_CLNDN","CLINVAR_CLNREVSTAT",
    "ACMG_SF_v3.2","REF_AA_dbnsfp","ALT_AA_dbnsfp","hg38_chr","hg38_pos(1-based)","cds_strand","refcodon","codonpos",
    "codon_degeneracy","SIFT_score","SIFT_converted_rankscore","SIFT_pred","LRT_score",
    "LRT_converted_rankscore","LRT_pred","LRT_Omega","MutationTaster_score",
    "MutationTaster_converted_rankscore","MutationTaster_pred","MutationTaster_model",
    "MutationTaster_AAE","MutationAssessor_score","MutationAssessor_rankscore",
    "MutationAssessor_pred","FATHMM_score","FATHMM_converted_rankscore","FATHMM_pred",
    "PROVEAN_score","PROVEAN_converted_rankscore","PROVEAN_pred","MetaSVM_score",
    "MetaSVM_rankscore","MetaSVM_pred","MetaLR_score","MetaLR_rankscore","MetaLR_pred",
    "Reliability_index","M-CAP_score","M-CAP_rankscore","M-CAP_pred","MutPred_score",
    "MutPred_rankscore","MutPred_protID","MutPred_AAchange","MutPred_Top5features",
    "fathmm-MKL_coding_score","fathmm-MKL_coding_rankscore","fathmm-MKL_coding_pred",
    "fathmm-MKL_coding_group","Eigen-raw_coding","Eigen-phred_coding",
    "Eigen-PC-raw_coding","Eigen-PC-phred_coding","Eigen-PC-raw_coding_rankscore",
    "integrated_fitCons_score","integrated_fitCons_rankscore","integrated_confidence_value",
    "GERP++_NR","GERP++_RS","GERP++_RS_rankscore","gnomAD_exomes_AC","gnomAD_exomes_AN",
    "gnomAD_exomes_AF","gnomAD_exomes_AFR_AC","gnomAD_exomes_AFR_AN","gnomAD_exomes_AFR_AF",
    "gnomAD_exomes_AMR_AC","gnomAD_exomes_AMR_AN","gnomAD_exomes_AMR_AF",
    "gnomAD_exomes_ASJ_AC","gnomAD_exomes_ASJ_AN","gnomAD_exomes_ASJ_AF",
    "gnomAD_exomes_EAS_AC","gnomAD_exomes_EAS_AN","gnomAD_exomes_EAS_AF",
    "gnomAD_exomes_FIN_AC","gnomAD_exomes_FIN_AN","gnomAD_exomes_FIN_AF",
    "gnomAD_exomes_NFE_AC","gnomAD_exomes_NFE_AN","gnomAD_exomes_NFE_AF",
    "gnomAD_exomes_SAS_AC","gnomAD_exomes_SAS_AN","gnomAD_exomes_SAS_AF",
    "gnomAD_genomes_AC","gnomAD_genomes_AN","gnomAD_genomes_AF","gnomAD_genomes_AFR_AC",
    "gnomAD_genomes_AFR_AN","gnomAD_genomes_AFR_AF","gnomAD_genomes_AMR_AC",
    "gnomAD_genomes_AMR_AN","gnomAD_genomes_AMR_AF","gnomAD_genomes_ASJ_AC",
    "gnomAD_genomes_ASJ_AN","gnomAD_genomes_ASJ_AF","gnomAD_genomes_EAS_AC",
    "gnomAD_genomes_EAS_AN","gnomAD_genomes_EAS_AF","gnomAD_genomes_FIN_AC",
    "gnomAD_genomes_FIN_AN","gnomAD_genomes_FIN_AF","gnomAD_genomes_NFE_AC",
    "gnomAD_genomes_NFE_AN","gnomAD_genomes_NFE_AF","Interpro_domain","GTEx_V8_gene",
    "GTEx_V8_tissue","MIM_id","Gene_old_names","Gene_full_name","Pathway(Uniprot)",
    "Pathway(BioCarta)_short","Pathway(BioCarta)_full","Pathway(ConsensusPathDB)",
    "Pathway(KEGG)_id","Pathway(KEGG)_full","Function_description","Disease_description",
    "MIM_phenotype_id","MIM_disease","Trait_association(GWAS)","GO_biological_process",
    "GO_cellular_component","GO_molecular_function","Tissue_specificity(Uniprot)",
    "Expression(egenetics)","Expression(GNF/Atlas)","Interactions(IntAct)",
    "Interactions(BioGRID)","Interactions(ConsensusPathDB)","P(HI)","P(rec)",
    "Known_rec_info","RVIS_EVS","RVIS_percentile_EVS","LoF-FDR_ExAC","RVIS_ExAC",
    "RVIS_percentile_ExAC","GHIS","GDI","GDI-Phred",
    "Gene_damage_prediction(all_disease-causing_genes)",
    "Gene_damage_prediction(all_Mendelian_disease-causing_genes)",
    "Gene_damage_prediction(Mendelian_AD_disease-causing_genes)",
    "Gene_damage_prediction(Mendelian_AR_disease-causing_genes)",
    "Gene_damage_prediction(all_PID_disease-causing_genes)",
    "Gene_damage_prediction(PID_AD_disease-causing_genes)",
    "Gene_damage_prediction(PID_AR_disease-causing_genes)",
    "Gene_damage_prediction(all_cancer_disease-causing_genes)",
    "Gene_damage_prediction(cancer_recessive_disease-causing_genes)",
    "Gene_damage_prediction(cancer_dominant_disease-causing_genes)"
]
output_frame=pd.DataFrame(columns=columns)

INPUT = [    
    os.path.join(INPUT_DIR, f) for f in os.listdir(INPUT_DIR)
    if os.path.isfile(os.path.join(INPUT_DIR, f)) and f.lower().endswith((".xls", ".xlsx"))
]

def sheetnamechecker (file):
    if "SNP_Indel_ANNO" in pd.ExcelFile(file).sheet_names:
        return "SNP_Indel_ANNO"
    else:
        return None
logging.info("The processing has been started.....")

for i, file in enumerate(INPUT, start=1): 
    start_time = time.time()
    logging.info(f"[{i}/{len(INPUT)}] Processing file: {os.path.basename(file)}")
    if not os.path.exists(file):
        logging.error(f"File not found: {file}") #logging
        continue
    try:
        sheetname=sheetnamechecker(file)
        if not sheetname:
            logging.error(f"No valid sheet found in {file}")
            continue
        sample = pd.read_excel(file, sheet_name=sheetname, usecols=["Gene_Name"])
        selected_rows = sample["Gene_Name"].isin(genes)

        # Đọc lại toàn bộ dữ liệu, nhưng chỉ những hàng cần thiết
        raw_data = pd.read_excel(file, sheet_name=sheetname, skiprows=lambda x: x > 0 and not selected_rows.iloc[x-1])
        raw_data = raw_data.iloc[1:]
        gene_filter = raw_data
        header_row = pd.DataFrame(
            [[file] + [""] * (len(gene_filter.columns) - 1)],
            columns=gene_filter.columns
        )
        if gene_filter.empty:
            logging.info(f"No variants found for gene '{GENE_NAME}' in {os.path.basename(file)}")           
        else:
            combined = pd.concat([header_row, gene_filter], ignore_index=True)

        if not combined.empty:
            output_frame = pd.concat([output_frame, combined], ignore_index=True) 
        elapsed = time.time() - start_time
        logging.info(f"Processed {os.path.basename(file)} in {elapsed:.2f} seconds")
    except Exception as e:
        logging.error(f"Error while processing {os.path.basename(file)}: {e}")

output_frame = output_frame.map(
    lambda x: ", ".join(map(str, x)) if isinstance(x, (list, tuple, np.ndarray)) else x
)
logging.info("Writing results to Excel...")

with pd.ExcelWriter(f"{OUTPUT}", engine="xlsxwriter") as writer:
    output_frame = output_frame.fillna(".")
    output_frame.to_excel(writer, sheet_name="Variant_screening", index=False)
    workbook = writer.book
    worksheet = writer.sheets["Variant_screening"]
    worksheet.autofilter(0, 0, 0, len(output_frame.columns)-1)

    merge_format = workbook.add_format({
        'align': 'left',
        'valign': 'vcenter',
        'font_name': 'Arial',
        'font_size': 10,
        'border': 1,
        'bg_color': "#C0DBEC"
    })

    header_format = workbook.add_format({
        'bold': True,
        'font_color': 'black',
        'font_size': 9,  
        'bg_color': "#ADCAE6",
        'border': 1,
        'font_name': 'Arial',
        'align': 'center',
        'valign': 'vcenter',
        'text_wrap': True,
    })

    data_format = workbook.add_format({
        'font_name': 'Arial',
        'font_color': 'black',
        'font_size': 9,  
        'bold': False,
        'text_wrap': False, 
        'align': 'general' 
    })
    # Apply formatting to header row
    for col_num, value in enumerate(output_frame.columns):
        worksheet.write(0, col_num, value, header_format)
    row_height = 7 * 15
    worksheet.set_row(0, row_height)
    # Set equal column widths
    # equal_width = 20
    # worksheet.set_column(0, len(output_frame.columns) - 1, equal_width)

    # Apply formatting to data rows and path like rows
    row_offset = 1 
    for i in range(len(output_frame)):
        first_col = str(output_frame.iloc[i, 0])
        other_cols = output_frame.iloc[i, 1:]
        row_idx = i + row_offset

        first_col_str = str(first_col).strip()
        is_path_like = (
            first_col_str.startswith("/") or
            first_col_str.startswith("./") or
            ":" in first_col_str[:5] or
            first_col_str.endswith((".xlsx", ".xls", ".csv"))
        )
        is_others_empty = all(
            (pd.isna(x) or str(x).strip() in ["", "."]) for x in other_cols
        )

        if is_path_like and is_others_empty:
            worksheet.set_row(row_idx, 15)
        else:
            worksheet.set_row(row_idx, None, data_format)

logging.info(f"Results saved to: {OUTPUT}")
