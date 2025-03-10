# coding=utf-8
# Author: Pham Huy Chau Long
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import pandas as pd

def parse_arguments():
    """
    Parse command-line arguments for the script.

    Returns:
        argparse.Namespace: Parsed command-line arguments.

    Examples:
        >>> python script.py -i data.parquet -o result.tsv
        ... args.input == 'data.parquet'
        ... args.output == 'result.tsv'
    """
    parser = argparse.ArgumentParser(
        description="Convert a parquet file to a tab-delimited (TSV) file."
    )
    parser.add_argument(
        '-i', '--input',
        default='report.parquet',
        type=str,
        help='Path to the input parquet file (default: report.parquet)'
    )
    parser.add_argument(
        '-o', '--output',
        default='report.tsv',
        type=str,
        help='Path to the output TSV file (default: report.tsv)'
    )
    return parser.parse_args()

def combine_intensities(df: pd.DataFrame) -> pd.DataFrame:
    """
    Combine the intensities from the individual fractions into a single column.

    Args:
        df (pd.DataFrame): The input DataFrame with the intensities in separate columns.
    
    Returns:
        pd.DataFrame: The DataFrame with an additional 'Intensities' column that combines 
                      the specified fraction columns.
        
    Examples:
        >>> df = pd.DataFrame({
        ...     'Fr.0.Quantity': [1, 2, 3],
        ...     'Fr.1.Quantity': [4, 5, 6],
        ...     'Fr.2.Quantity': [7, 8, 9],
        ...     # include all required columns...
        ... })
        >>> df = combine_intensities(df)
        >>> df['Intensities']
        0    1;4;7;...
        1    2;5;8;...
        2    3;6;9;...
    """
    required_columns = [
        'Fr.0.Quantity', 'Fr.1.Quantity', 'Fr.2.Quantity', 
        'Fr.3.Quantity', 'Fr.4.Quantity', 'Fr.5.Quantity', 
        'Fr.6.Quantity', 'Fr.7.Quantity', 'Fr.8.Quantity', 
        'Fr.9.Quantity', 'Fr.10.Quantity', 'Fr.11.Quantity'
    ]
    
    # Check if all required columns exist in the DataFrame
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError("Check the --export-quant option in DIA-NN. Missing required columns: " + ", ".join(missing_columns))
    
    # Combine the intensities into a single column
    df['Intensities'] = df[required_columns].astype(str).replace('^0\\.0*$','0', regex=True).agg(';'.join, axis=1)
    return df

def main():
    args = parse_arguments()

    # Provide feedback
    print(f"Input file  = {args.input}")
    print(f"Output file = {args.output}")

    # Read the parquet file
    df = pd.read_parquet(args.input)

    # Process the DataFrame
    df = combine_intensities(df)

    # Write to TSV file
    df.to_csv(args.output, sep='\t', index=False, quoting=False)

if __name__ == '__main__':
    main()