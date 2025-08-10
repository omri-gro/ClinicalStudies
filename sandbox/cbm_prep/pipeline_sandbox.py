import sandbox as sb

if __name__ == "__main__":
    metadata = sb.MetadataBundle('config.yaml')
    cpg_test_df_raw = sb.raw_to_df('CPG_example.xlsx', 'CPG', 'CBM')
    cpg_omr_df_raw = sb.raw_to_df('CPG_OMR.csv', 'CPG', 'OMR')

    cpg_test_df_curated = sb.curate_df(cpg_test_df_raw, metadata)
    cpg_omr_df_curated = sb.curate_df(cpg_omr_df_raw, metadata)

    print(cpg_omr_df_raw)
