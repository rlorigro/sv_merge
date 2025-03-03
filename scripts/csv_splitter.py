import collections
import argparse
import csv
import os


def main(
    file_path,
    output_folder,
    bucket_name
):
    rows_by_sample_name = collections.defaultdict(list)

    # Open the input CSV file
    with open(file_path, 'r', newline='\n') as input_file:
        # Create a dictionary to store the rows for each sample name

        # Iterate over the rows in the CSV file
        for line in input_file:
            # Split the line into a list of values
            row = line.strip().split(',')

            # Split the name at the underscore to get the sample name
            sample_name = row[0].split('_')[0]

            # Add the row to the list for this sample name
            rows_by_sample_name[sample_name].append(row)

    # create output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # create output path with dir
    output_path = os.path.join(output_folder, 'upload.sh')

    # Open the shell script file
    with open(output_path, 'w') as script_file:
        # Iterate over the items in the dictionary
        for sample_name, rows in rows_by_sample_name.items():
            csv_path = os.path.join(output_folder, sample_name + '.csv')

            # Open a new CSV file with the sample name as the filename
            with open(csv_path, 'w') as output_file:
                for row in rows:
                    output_file.write(','.join(row) + '\n')

            # Write a line to the shell script to upload the CSV file to the GCS bucket
            script_file.write(f'gsutil cp {csv_path} gs://{bucket_name}/{sample_name}_haps_vs_chm13_bams.csv\n')


if __name__ == "__main__":
    # add an argparser which accepts the gs URI folder to upload to
    parser = argparse.ArgumentParser(description='Split a CSV file by sample name and generate a shell script to upload the files to a GCS bucket.')

    parser.add_argument('--input', type=str, help='The path to the input CSV file.', required=True)
    parser.add_argument('--output_folder', type=str, help='The path to the output folder.', required=True)
    parser.add_argument('--bucket_name', type=str, help='The name of the GCS bucket to upload the files to.', required=True)

    args = parser.parse_args()

    main(args.input, args.output_folder, args.bucket_name)

