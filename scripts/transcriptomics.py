import pandas as pd
import cobra
import os.path

TRANSCRIPTOMICS_FN = "../data/SOP-TS9-F452 data.csv"

def read_transcriptomics(trans_fn = TRANSCRIPTOMICS_FN):
    df = pd.read_csv(trans_fn, sep = ";")
    print(df.columns)
    df.drop("id_0", inplace = True)
    # Remove rows that does not have a sco number
    df = df[~df["gene ID"].isnull()]
    print(df.head)


if __name__ == '__main__':
    read_transcriptomics()