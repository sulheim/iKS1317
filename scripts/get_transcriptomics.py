import pandas as pd
import cobra
import scipy

TRANSCRIPTOMICS_FN = "../data/SOP-TS9-F452 data.csv"
MODEL_FN = "../iKS1317.xml"
OFFLINE_DATA_FN = "../data/offline_data_F452.csv"

def read_transcriptomics(trans_fn = TRANSCRIPTOMICS_FN):
    df = pd.read_csv(trans_fn, sep = ";", decimal = ",")
    df.drop(["id_0", "Unnamed: 28"], axis = 1, inplace = True)

    # Remove rows that does not have a sco number
    df = df[~df["gene ID"].isnull()]
    df = df[df["gene ID"].str.contains("SC")]

    df.set_index("gene ID", inplace = True)
    df.columns = [float(x.split("_")[1]) for x in df.columns]
    # Split rows with more than one gene
    splitted_rows = []
    for idx, row in df[df.index.str.contains(",")].iterrows():
        for x in idx.split(","):
            splitted_rows.append([x] + list(row))
    new_df = pd.DataFrame(splitted_rows).set_index(0)
    new_df.columns = df.columns

    df = df[~df.index.str.contains(",")]
    df = pd.concat([new_df, df])


    return df


def get_model(model_fn = MODEL_FN):
    return cobra.io.read_sbml_model(model_fn)

def plot_transcriptomics_for_model_genes():
    df = read_transcriptomics()
    model = get_model()

    model_genes = [g.id for g in model.genes]
    model_genes.remove("s0001")


    print(len(df.index))
    print(df.columns)
    plt.ion()
    fig, ax = plt.subplots(figsize = (16, 10))
    cmap = plt.get_cmap("tab20")
    x = np.linspace(df.columns[0], df.columns[-1], 100)
    for i, gene_id in enumerate(model_genes):
        
        data = data_for_gene(df. gene_id)
        print(gene_id)
        # print(data)
        # try:
        #     spline = UnivariateSpline(data.index, data, k = 2, s = 0.5)
        # except:
        #     print(data)
        # ax.plot(x, spline(x), ls = "--", c = cmap(i%10), label = gene_id)
        # roll = data.rolling(3, min_periods = 1).mean()
        # ax.plot(roll, ls = "--", c = cmap(i%10), label = gene_id)

        data.index = pd.to_timedelta(data.index, "h")
        data = data.resample("10T").interpolate(method = "linear")
        data.plot(ax = ax, c = cmap(i%10), label = gene_id)
        b, a = scipy.signal.butter(N = 3, Wn = 0.05, btype = "lowpass") # It is the Wn parameter which defines the cut off 
        out = scipy.signal.filtfilt(b, a, data)
        ax.plot(data.index, out, ls = "--", c = cmap(i%10), label = gene_id)
        # plt.show()
        plt.draw()
        plt.pause(0.5)
        if i%5 == 0:
            ax.cla()

def data_for_gene(df, gene_id):
    if gene_id in ["SCP1233B", "SCP1233", "SCP153c", "SCP1301", "SCP155c", "SCP1299"]:
        key =  ".".join([gene_id[:4], gene_id[4:]])
    else:
        key = gene_id
    try:
        data = df.loc[key]
    except KeyError:
        print("No gene expression data for {0}".format(key))
        return None

    if isinstance(data, pd.core.frame.DataFrame):
        print("{0} rows for gene {1}".format(len(data.index), gene_id))
        data = data.agg("mean")
    return data

def resample_and_filter(data, resample_period, butter_wn, gene_id):
    data.index = pd.to_timedelta(data.index, "h")
    data = data.resample("10T").interpolate(method = "linear")
    b, a = scipy.signal.butter(N = 3, Wn = 0.05, btype = "lowpass") # It is the Wn parameter which defines the cut off 
    out = pd.Series(scipy.signal.filtfilt(b, a, data), name = gene_id)
    out.index = data.index
    return out


def filter_and_resample_and_derivative(resample_period = "10T", butter_wn = 0.05):
    df = read_transcriptomics()
    model = get_model()
    model_genes = [g.id for g in model.genes]
    model_genes.remove("s0001")

    resampled_and_smoothed = []
    for i, gene_id in enumerate(model_genes):
        data = data_for_gene(df, gene_id)
        if not data is None:
            smooth_data = resample_and_filter(data, resample_period, butter_wn, gene_id)
            resampled_and_smoothed.append(smooth_data)
    
    smooth_df = pd.concat(resampled_and_smoothed, axis = 1).T
    timedelta = ((smooth_df.columns[1]-smooth_df.columns[0]).seconds / 3600) # In hours

    derivate_df = smooth_df.diff(1, axis = 1) / timedelta
    return derivate_df, model


def convert_gene_data_to_reaction_data(model, derivate_df):
    reaction_data_list = []
    for r in model.reactions:
        grr = r.gene_reaction_rule
        if not len(grr):
            continue

        or_groups = grr.split(" or ")
        or_list = []
        for grr_part in or_groups:
            and_group = grr_part.split(" and ")
            and_group = [x.strip("() ") for x in and_group]
            and_data = derivative_df[derivative_df.index.isin(and_group)]
            if isinstance(and_data, pd.core.frame.DataFrame):
                and_data = and_data.agg("mean", axis = 0)
            or_list.append(and_data)

        reaction_data = pd.concat(or_list, axis = 1).agg("sum", axis = 1)
        reaction_data.name = r.id
        reaction_data_list.append(reaction_data)

    reaction_data_df = pd.concat(reaction_data_list, axis = 1).T
    return reaction_data_df