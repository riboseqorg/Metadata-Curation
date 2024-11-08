import csv
import argparse
from Bio import Entrez
import os
import time
from collections import OrderedDict

def read_bioprojects(input_file):
    with open(input_file, 'r') as f:
        reader = csv.reader(f)
        return [row[0] for row in reader]

def get_sra_run_info(bioproject_id):
    # Set your email here. This is required by NCBI to track usage of their API
    Entrez.email = "riboseq@gmail.com"
    
    # Search for SRA entries related to the BioProject
    handle = Entrez.esearch(db="sra", term=f"{bioproject_id}[BioProject]")
    record = Entrez.read(handle)
    handle.close()
    
    # Get the list of IDs
    id_list = record["IdList"]
    
    # Fetch run info for each ID
    run_info_list = []
    for sra_id in id_list:
        handle = Entrez.efetch(db="sra", id=sra_id, rettype="runinfo", retmode="text")
        run_info = handle.read()
        handle.close()
        run_info_list.append(run_info)
        time.sleep(0.5)  # Be nice to NCBI servers
    
    return run_info_list

def parse_run_info(run_info_list):
    all_data = []
    headers = OrderedDict()

    for run_info in run_info_list:
        reader = csv.reader(run_info.splitlines())
        run_headers = next(reader)
        for header in run_headers:
            headers[header] = None
        all_data.extend(list(reader))

    return list(headers.keys()), all_data

def save_combined_run_info(all_headers, all_data, output_file):
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(all_headers)
        for row in all_data:
            # Pad the row with empty strings if it's shorter than all_headers
            padded_row = row + [''] * (len(all_headers) - len(row))
            writer.writerow(padded_row)

def main():
    parser = argparse.ArgumentParser(description="Download SRA run info for multiple BioProjects")
    parser.add_argument("input_file", help="Input CSV file with BioProject IDs")
    parser.add_argument("output_file", help="Output CSV file name")
    args = parser.parse_args()

    bioprojects = read_bioprojects(args.input_file)
    all_headers = []
    all_data = []

    for bioproject in bioprojects:
        print(f"Processing BioProject: {bioproject}")
        run_info_list = get_sra_run_info(bioproject)
        headers, data = parse_run_info(run_info_list)
        
        # Update all_headers with new unique headers
        all_headers.extend([h for h in headers if h not in all_headers])
        all_data.extend(data)

    save_combined_run_info(all_headers, all_data, args.output_file)
    print(f"Combined SRA run info for all BioProjects has been saved to {args.output_file}")

if __name__ == "__main__":
    main()