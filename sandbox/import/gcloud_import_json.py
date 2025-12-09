import subprocess
import shlex
from pandas import read_csv
from google.cloud import storage
from pathlib import Path
import os


def run_cmd(cmd):
    result = subprocess.run(
        shlex.split(cmd),
        capture_output=True,
        text=True,
        shell=False  # safer than shell=True
    )
    if result.returncode != 0:
        print("ERROR:", result.stderr)

def gs_cp(bucket, scan_id, src_suffix, dst_file):
    GSUTIL = r'"C:\Users\omrig\AppData\Local\Google\Cloud SDK\google-cloud-sdk\bin\gsutil.cmd"'
    src = f"gs://{bucket}/{scan_id}/PBS/*/results/{src_suffix}"
    cmd = f'{GSUTIL} cp "{src}" "{dst_file}"'
    run_cmd(cmd)


def download_blob(bucket_name, source_blob_name, destination_file_name):
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(source_blob_name)

    blob.download_to_filename(destination_file_name)

    print(f"Blob {source_blob_name} downloaded to {destination_file_name}.")


def download_single_file(bucket_name, scan_id, suffix, out_file):
    """
        Finds and downloads exactly one file for the given analysis_id whose blob
        name ends with the given suffix.
    """
    client = storage.Client()
    bucket = client.bucket(bucket_name)

    prefix = f"{scan_id}/PBS/"  # known portion before the unknown datetime folder
    full_suffix = f"/results/{suffix}"

    # Search for a blob ending with the right suffix
    match = None
    for blob in bucket.list_blobs(prefix=prefix):
        if blob.name.endswith(full_suffix):
            match = blob
            break

    if match is None:
        print(f"No file found for {scan_id} ending with '{full_suffix}'")
        return None

    # Make sure target directory exists
    Path(out_file).parent.mkdir(parents=True, exist_ok=True)

    match.download_to_filename(out_file)
    print(f"Downloaded {match.name} → {out_file}")
    return out_file

def dig_pbs_rslts_bckt(bucket_name, scan_id, save_dir):
    """
    Assumes bucket structure
    scopio_raw_inference_prod/<scan_id>/PBS/<datetime>/results/results.json
    scopio_raw_inference_prod/<scan_id>/PBS/<datetime>/results/pbs.log
    which will be saved into
    <save_dir>/pbs/<scan_id>.json
    <save_dir>/logs/<scan_id>.log
    """
    download_single_file(
        bucket_name,
        scan_id,
        suffix="results.json",
        out_file=f"{save_dir}/jsons/{scan_id}.json"
    )

    # download pbs.log → logs/<id>.log
    download_single_file(
        bucket_name,
        scan_id,
        suffix="pbs.log",
        out_file=f"{save_dir}/logs/{scan_id}.log"
    )



if __name__ == "__main__":
    """ currently not working, use command line instead"""
    bucket = 'scopio_raw_inference_prod'
    ids_csv = r'C:\Users\omrig\Downloads\SYN_scan_ids.csv'
    save_dir = r'C:\Users\omrig\PycharmProjects\pythonProject\CBM_verification\new_ssh\importing\SYN_from_cloud'

    ids_df = read_csv(ids_csv)
    for uuid in ids_df['Scan IDs total in disk']:
        uuid = uuid.strip()
        gs_cp(bucket, uuid, "results.json", fr"{save_dir}/pbs/{uuid}.json")
        gs_cp(bucket, uuid, "pbs.log", fr"{save_dir}/logs/{uuid}.log")



