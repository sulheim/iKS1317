import pandas as pd
import cobra
import os.path
from matplotlib import pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
import numpy as np
import scipy.signal
import re
# import tim

TRANSCRIPTOMICS_FN = "../data/SOP-TS9-F452 data.csv"
MODEL_FN = "../iKS1317.xml"
OFFLINE_DATA_FN = "../data/offline_data_F452.csv"

# MAX_GLUCOSE_UPTAKE_RATE = -2.1
MAX_GLUCOSE_UPTAKE_RATE = -2.1
# CORRESPONDING_AMMONIUM_UPTAKE = -1.85
CORRESPONDING_AMMONIUM_UPTAKE = -2

MOLAR_MASS = {
    "PO4": 94.9714, #g/mol
}

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

# def func(x, A, b, c):
#     return c - A*np.exp(b*x)

def sigmoid(x, a, b, c, d):
    y = a / (1 + np.exp(c*(x-b))) + abs(d)
    return y

def fit_PO4(df):
    df = df[~df["PO4"].isnull()]
    fit, _ = curve_fit(sigmoid, df.index, df["PO4"], p0 = [480, 40, 1, 10])
    print(fit)
    return fit



def read_offline_data():
    df = pd.read_csv(OFFLINE_DATA_FN, sep = ";", decimal = ",")
    df.drop("Sample #", axis = 1, inplace = True)
    df.loc[0, "CDW"] = 0
    df.loc[30, "RED":"TBP"] = 0

    df.set_index("TAI", inplace = True)
    po4_fit = fit_PO4(df)

    df.index = pd.to_timedelta(df.index, "h")
    df_resampled = df.resample("10T").interpolate(method = "linear")

    df_resampled["Fitted PO4"] = sigmoid(df_resampled.index.total_seconds()/3600, *po4_fit)



    # # glc_fit = np.poly1d(np.polyfit(df.index.total_seconds(), df["Glucose"], 3))
    # glc_fit2, _ = curve_fit(func, df.index.total_seconds()/3600, np.array(df["Glucose"]), p0 = [1, 2, 40])
    # y = func(df.index.total_seconds()/3600, *glc_fit2)


    # # print(glc_fit(df.index.total_seconds()))
    # plt.plot(df_resampled.index.total_seconds()/3600, df_resampled["Fitted PO4"])
    # plt.plot(df_resampled.index.total_seconds()/3600, df_resampled["PO4"])
    # plt.show()


    df_resampled[df_resampled.isnull()] = 0

    # Make max rate columns
    df_resampled["Max glc"] = df_resampled["Glucose"].diff() * 6 # per hour
    df_resampled["Max glu"] = df_resampled["Glutamate"].diff() * 6 # per hour
    df_resampled["Max PO4"] = df_resampled["Fitted PO4"].diff() * 6 # per hour

    return df_resampled


def dFBA_trans(model, reaction_data_df, growth_data):
    pass

def dFBA(model, growth_data):
    # initial conditions
    # glc = 41 g/L
    # Glu = 57 g/L
    # Initial biomass = 0.03 g/L

    time_array = growth_data.index#.total_seconds()/3600 # Hours
    dt = (time_array[1]-time_array[0]).total_seconds()/3600
    
    N = len(time_array)
    biomass_arr = np.zeros(N)
    biomass_arr[0] = 0.03 #g/L

    PO4_arr = np.zeros(N)
    PO4_arr[0] = 491.7e-3 #g/L
    
    glc_arr = np.zeros(N)
    glc_arr[0] = 41
    
    glu_arr = np.zeros(N)
    glu_arr[0] = 57

    
    # S_glc_uptake_array = np.zeros(N)
    # S_glu_uptake_array = np.zeros(N)
    S_pi_uptake_array = np.zeros(N)
    growth_rate_array = np.zeros(N)

    # Units in mmol/g
    # No limitation on glucose and glutamate
    scale_sources_3(model, ["EX_glc__D_e", "EX_glu__L_e"])

    S_pi_available = -1000

    for i, timepoint in enumerate(time_array):
        if i > 0:
            biomass_arr[i] = biomass_arr[i-1] * np.exp(growth_rate * dt) # g/L
            PO4_arr[i] = PO4_arr[i-1] + S_pi_uptake*biomass_arr[i-1]*dt
        
            S_pi_available = min(-0.01, get_available_PO4(growth_data, timepoint) / biomass_arr[i])

        model.reactions.EX_pi_e.lower_bound = min(0, S_pi_available)
        solution = model.optimize()
        
        S_glc_uptake = solution.x_dict["EX_glc__D_e"]
        S_pi_uptake = solution.x_dict["EX_pi_e"]
        S_glut_uptake = solution.x_dict["EX_glu__L_e"]
        growth_rate = solution.objective_value

        growth_rate_array[i] = growth_rate
        S_pi_uptake_array[i] = S_pi_uptake
        

    fig, [ax1, ax2, ax3] = plt.subplots(1,3)
    ax1.plot(time_array.total_seconds()/3600, biomass_arr, c = "k", lw = 5, label = "Biomass")
    ax1.plot(time_array.total_seconds()/3600, growth_data["CDW"], lw = 4, label = "CDW")
    ax2.plot(time_array.total_seconds()/3600, S_pi_uptake_array)
    ax2.plot(time_array.total_seconds()/3600, growth_rate_array/max(growth_rate_array))
    ax2.plot(time_array.total_seconds()/3600, S_pi_uptake_array/max(S_pi_uptake_array))
    ax3.plot(time_array.total_seconds()/3600, growth_rate_array)
    # plt.legend()
    plt.show()


def get_available_PO4(growth_data, timepoint):
    max_PO4 = growth_data.loc[timepoint, "Max PO4"] #g/L/hr
    mmol_max_PO4 = 1e3*max_PO4/MOLAR_MASS["PO4"] #mmol/L/hr
    return mmol_max_PO4

    # except 
def scale_sources_3(model, sources):
    total_carbon_uptake = MAX_GLUCOSE_UPTAKE_RATE * 6 # Glucose has 6 carbons
    total_nitrogen_uptake = CORRESPONDING_AMMONIUM_UPTAKE * 1 # Ammonium has 1 nitrogen

    carbon_bound_exp = 0
    nitrogen_bound_exp = 0
    for source_id in sources:
        # print(source_id)
        source = model.reactions.get_by_id(source_id)
        n_carbon, n_nitrogen = get_number_of_carbons_and_nitrogens_in_metabolite(list(source.metabolites)[0])
        source.lower_bound = max(total_carbon_uptake / max(n_carbon, 1e-9), total_nitrogen_uptake / max(n_nitrogen, 1e-9))
        carbon_bound_exp += source.flux_expression * n_carbon
        nitrogen_bound_exp += source.flux_expression * n_nitrogen
    carbon_bound = model.problem.Constraint(carbon_bound_exp, lb = total_carbon_uptake)
    nitrogen_bound = model.problem.Constraint(nitrogen_bound_exp, lb = total_nitrogen_uptake)

    model.add_cons_vars(carbon_bound)
    model.add_cons_vars(nitrogen_bound)


def get_number_of_carbons_and_nitrogens_in_metabolite(metabolite):
    nitrogen_match = re.search(r"(N)(\d*)", metabolite.formula)
    if nitrogen_match is not None:
        if not nitrogen_match.group(2):
            n_nitrogen = 1
        else:
            n_nitrogen = int(nitrogen_match.group(2))
    else:
        n_nitrogen = 0

    # print(metabolite.formula)

    carbon_match = re.search(r"(C)(\d*)", metabolite.formula)
    if carbon_match is not None:
        if not carbon_match.group(2):
            n_carbon = 1
        else:
            n_carbon = int(carbon_match.group(2))
    else:
        n_carbon = 0
    return n_carbon, n_nitrogen

if __name__ == '__main__':
    # read_transcriptomics()
    # plot_transcriptomics_for_model_genes()

    if 0:
        derivative_df, model = filter_and_resample_and_derivative()
        derivative_df.to_csv("../data/derivative_transcr_genes.csv")

    if 0:
        derivative_df = pd.read_csv("../data/derivative_transcr_genes.csv", index_col = 0)
        # print(derivative_df)
        model = get_model()
        reaction_data_df = convert_gene_data_to_reaction_data(model, derivative_df)
        reaction_data_df.to_csv("../data/derivatve_transcr_rxns.csv")

    if 1:
        offline_df = read_offline_data()

        reaction_data_df = pd.read_csv("../data/derivatve_transcr_rxns.csv", index_col = 1)
        model = get_model()
        dFBA(model, offline_df)

    # print(derivative_df.head())