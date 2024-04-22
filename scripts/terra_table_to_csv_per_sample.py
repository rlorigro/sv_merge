import argparse
import os


def main(
        tsv_path,
        output_folder,
        bucket_name,
        column_names: list,
        new_column_name: str
):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    script_path = os.path.join(output_folder, 'upload.sh')
    column_indices = list()

    gs_uris = list()

    with open(script_path, 'w', newline='\n') as script_file:
        # Parse the tsv file extract each of the columns indicated in the list of column_names, as well as their 0th column
        with open(tsv_path, 'r') as tsv_file:
            for l,line in enumerate(tsv_file):
                row = line.strip().split('\t')

                if l == 0:
                    for column in column_names:
                        if column not in row:
                            raise ValueError(f'Column {column} not found in header: {row}')

                        column_indices.append(row.index(column))

                    continue

                # Create a new CSV file for each row
                filename = f'{row[0]}_{new_column_name.strip("csv").strip("_")}.csv'

                with open(os.path.join(output_folder, filename), 'w', newline='\n') as output_file:
                    # write one line in the output CSV for each column in column_names
                    if len(column_indices) == 1:
                        output_file.write(row[0])
                        output_file.write(',')
                        output_file.write(row[column_indices[0]] + '\n')
                    else:
                        for i, column_index in enumerate(column_indices):
                            output_file.write(row[0] + '_' + str(i+1))
                            output_file.write(',')
                            output_file.write(row[column_index])
                            output_file.write('\n')

                # Write a line to the shell script to upload the CSV file to the GCS bucket
                gs_uri = f"gs://{bucket_name.strip('/')}/{filename}"
                script_file.write(f'gsutil cp {os.path.join(output_folder, filename)} {gs_uri}\n')

                gs_uris.append(gs_uri)

    # now write a modified version of the input TSV that only contains column[0] and a column indicating the gs uris for the CSVs
    with open(tsv_path, 'r') as tsv_file:
        with open(os.path.join(output_folder, 'modified.tsv'), 'w', newline='\n') as modified_tsv_file:
            for l,line in enumerate(tsv_file):
                row = line.strip().split('\t')

                if l == 0:
                    modified_tsv_file.write(row[0] + '\t' + new_column_name + '\n')
                    continue

                modified_tsv_file.write(row[0] + '\t' + gs_uris[l-1] + '\n')


if __name__ == "__main__":
    # add an argparser which accepts the gs URI folder to upload to
    parser = argparse.ArgumentParser(description='Split rows of a TSV file and generate a shell script to upload the files to a GCS bucket.')

    parser.add_argument('--tsv_path', type=str, help='The path to the input CSV file.', required=True)
    parser.add_argument('--output_folder', type=str, help='The path to the output folder.', required=True)
    parser.add_argument('--bucket_name', type=str, help='The name of the GCS bucket to upload the files to.', required=True)
    parser.add_argument('--column_names', type=str, help='names of columns to put into CSVs', required=False)
    parser.add_argument('--new_column_name', type=str, help='name of new column in Tera TSV', required=True)

    args = parser.parse_args()

    # parse column names as list
    args.column_names = args.column_names.split(',')

    main(args.tsv_path, args.output_folder, args.bucket_name, args.column_names, args.new_column_name)
