import pandas as pd
import argparse


def delta(col1, col2, prec):
    df = pd.DataFrame()
    df['diff'] = abs(col1 - col2)
    return df

def test_differences(parquet_file, tsv_file):

    parquet = pd.read_parquet(parquet_file)
    tsv = pd.read_csv(tsv_file, sep = '\t', header = 0)
    different = False

    # Blob to String
    for name in parquet.select_dtypes(include='object').columns.tolist():
        parquet[name] = parquet[name].str.decode('utf-8')

    parquet = parquet.sort_values(by = ['query','target'],ignore_index = True)
    tsv = tsv.sort_values(by = ['query','target'],ignore_index = True)

    # Check row numbers
    if (parquet.size != tsv.size):
        print("Different table size")
        print(parquet.size(), tsv.size())
        different = True


    # Check for NaNs difference
    if (parquet.isna().values.any() ^ tsv.isna().values.any()):
        print("Dissagrement on existance of NaN values")
        print(parquet.isna().any())
        print(tsv.isna().any())
        different = True

    # Compare values
    for p,t in zip(parquet.select_dtypes(exclude='object').columns.tolist(), tsv.select_dtypes(exclude='object').columns.tolist()):
        prec = 0.001
        dff = delta(parquet[p],tsv[t],prec)

        if sum(dff['diff'] >= prec) > 0:
            parq = parquet[dff['diff'] >= prec][p]
            tsvs = tsv[dff['diff'] >= prec][t]
            pout = pd.concat([parquet[dff['diff'] >= prec]["query"], parquet[dff['diff'] >= prec]["target"], parq,tsvs], axis = 1)
            pout.columns=["query","target",f"Parquet {p}", f"TSV {t}"]
            print(f"Difference bigger than {prec} in column {p}")
            print(pout.head())
            different = True

    # Compare Strings
    for p,t in zip(parquet.select_dtypes(include='object').columns.tolist(), tsv.select_dtypes(include='object').columns.tolist()):
        diff = parquet[p] != tsv[t]
        if (sum(diff) > 0):
            parq = parquet[diff][p]
            tsvs = tsv[diff][t]
            pout = pd.concat([parq,tsvs], axis = 1)
            pout.columns=[f"Parquet {p}", f"TSV {t}"]
            print(pout.head())
            different = True

    print(parquet)
    print(tsv)
    
    if different == False:
        print("No difference in table values found")
        return 1
    else:
        return 0


def main():
    parser = argparse.ArgumentParser(description="Checks for table differences")
    parser.add_argument("parquet", help="Parquet file", type=str)  # Positional
    parser.add_argument("tsv", help="TSV file", type=str)
    args = parser.parse_args()
    return test_differences(args.parquet, args.tsv)

main()