import os
import json
import csv
from pathlib import Path


def map_spherocytes_to_csv(base_directory, output_csv_filename):
    """
    Iterates through subdirectories, reads results.json, extracts the
    spherocytes percentage, and saves the mapping to a CSV.
    """
    base_path = Path(base_directory)
    extracted_data = []

    # Iterate through all items in the base directory
    for subdir in base_path.iterdir():
        if subdir.is_dir():
            json_path = subdir / "results.json"

            # Check if results.json exists in this subdirectory
            if json_path.exists():
                try:
                    with open(json_path, 'r', encoding='utf-8') as file:
                        data = json.load(file)

                    spherocytes_percentage = "Not Found"

                    # Navigate through distributions to find observations
                    distributions = data.get("distributions", [])
                    for dist in distributions:
                        observations = dist.get("observations", [])

                        # Find the spherocytes observation
                        for obs in observations:
                            if obs.get("name") == "spherocytes":
                                spherocytes_percentage = obs.get("percentage", "Not Found")
                                break

                    # Append the mapping to our list
                    extracted_data.append([subdir.name, spherocytes_percentage])

                except json.JSONDecodeError:
                    print(f"Skipping {json_path}: Invalid JSON file.")
                except Exception as e:
                    print(f"Error processing {json_path}: {e}")

    # Write the extracted data to a CSV file
    try:
        with open(output_csv_filename, 'w', newline='', encoding='utf-8') as csv_file:
            writer = csv.writer(csv_file)
            # Write the header
            writer.writerow(["Subdirectory_Name", "Spherocytes_Percentage"])
            # Write the data rows
            writer.writerows(extracted_data)
        print(f"Successfully saved results to {output_csv_filename}")
    except Exception as e:
        print(f"Failed to write to CSV: {e}")


if __name__ == "__main__":
    # Replace these variables with your actual directory path and desired output filename
    TARGET_DIRECTORY = r"S:\talm\cbm_clinical_trial\Spherocyte_3.13\roi_results_jsons"
    OUTPUT_CSV = "results/spherocytes_313_results.csv"

    map_spherocytes_to_csv(TARGET_DIRECTORY, OUTPUT_CSV)
