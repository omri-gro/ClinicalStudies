import sandbox as sb
import os


if __name__ == "__main__":
    analysis_name = 'bwh100_investigator_training'
    site = 'BWH'

    raw_dir = os.path.abspath(os.path.dirname(__file__))
    raw_dir = os.path.join(raw_dir, r'raw', analysis_name)

    metadata = sb.MetadataBundle('config.yaml')
    clv_df_raw = sb.raw_to_df('allCRFsClV.csv', site, 'ClV', dir=raw_dir)
    mnl_df_raw = sb.raw_to_df('allCRFsManual.csv', site, 'Manual', dir=raw_dir)
    sb3163_df_raw = sb.raw_to_df('sb3163csv.csv', 'sb3163', 'CBM', dir=raw_dir)
    sb3209_df_raw = sb.raw_to_df('sb3209csv.csv', 'sb3209', 'CBM', dir=raw_dir)

    clv_df = sb.curate_df(clv_df_raw, metadata, 'ClV')
    mnl_df = sb.curate_df(mnl_df_raw, metadata, 'Manual', wbcs_as_counts=True)
    sb3163_df = sb.curate_df(sb3163_df_raw, metadata, 'Scopio')
    sb3209_df = sb.curate_df(sb3209_df_raw, metadata, 'Scopio')


    print('')

