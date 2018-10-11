from matplotlib import pyplot as plt
import matplotlib
import pandas as pd
from pathlib import Path
import cobra
from transcriptomics import read_offline_data

MODEL_FN = Path("../iKS1317.xml")
OFFLINE_DATA_FN = "../data/offline_data_F452.csv"

plt.style.use('seaborn-white')
matplotlib.rcParams["font.size"] = 16
matplotlib.rcParams["font.weight"] = "bold"
matplotlib.rcParams["axes.labelweight"] = "bold"
matplotlib.rcParams["lines.linewidth"] = 3
matplotlib.rcParams['axes.titleweight'] = "bold"
matplotlib.rcParams['axes.titlesize'] = 18

def plot_dFBA_results(dFBA_fn, growth_data_fn):
    dFBA_df = pd.read_csv(dFBA_fn, sep = ";")
    growth_data_df = read_offline_data(growth_data_fn)


    fig, axes = plt.subplots(2,2, figsize = (20, 10))
    [ax1, ax2, ax3, ax4] = axes.flatten()
    dFBA_df.plot(x = "Hours", y = "Biomass", ax = ax1, c = "k")
    ax1.plot(dFBA_df["Hours"], growth_data_df["CDW"], c = "b", lw = 4, label = "Measured CDW")

    plt.show()

    # dFBA_df.plot(x = "Hours", y = "Glucose uptake", ax = ax2, c = "b")
    # dFBA_df.plot(x = "Hours", y = "Glutamate uptake", ax = ax2, c = "r")
    

    # dFBA_df.plot(x = "Hours", y = "PO4 uptake", ax = ax3, c = "r")
    # dFBA_df.plot(x = "Hours", y = "Max PO4 /gDW", ax = ax3, c = "b")
    # plt.show()

        # self.dFBA_df["Max PO4"].plot(ax = ax4, c = "b")
        # ax4.plot(self.time_array_hours, self.growth_data_df["Max PO4"] / self.model.metabolites.pi_e.formula_weight, c = "b")
        

        # ax1.plot(self.dFBA_df["Biomass"], c = "k", lw = 5, label = "Biomass")
        # ax1.plot(self.dFBA_df["Biomass"]["CDW"], lw = 4, label = "CDW")
        # ax2.plot(self.dFBA_df["Biomass"], c = "k", label = "PO4 uptake")
        # ax2.plot(self.dFBA_df["Biomass"], label = "PO4 uptake")

        # ax3.plot(self.dFBA_df["Biomass"])
        # ax4.plot(self.dFBA_df["Biomass"], c = "b", lw = 5)
        # ax4.plot(self.dFBA_df["Biomass"], c = "r", lw = 5)

        # ax5.plot(self.dFBA_df["Biomass"], c = "r", lw = 5)
        # ax5.plot(self.dFBA_df["Biomass"]["Fitted PO4"]*1e3/model.metabolites.pi_c.formula_weight, c = "k", lw = 5)
        # plt.legend()

def plot_offline_data():
    growth_data_df = read_offline_data(growth_data_fn)
    dFBA_df = pd.read_csv(dFBA_fn, sep = ";")

    fig, ax = plt.subplots(1, figsize = (8, 14))
    dFBA_df.plot(x = "Hours", y = "Biomass", ax = ax1, c = "k")
    ax1.plot(dFBA_df["Hours"], growth_data_df["CDW"], c = "b", label = "Measured CDW")
    ax1.plot(dFBA_df["Hours"], growth_data_df["PO4"], c = "r", label = "Measured phosphate")
    print(growth_data_df.columns)
    

def plot_exchange_results(fn, dFBA_fn, growth_data_fn = OFFLINE_DATA_FN):
    growth_data_df = read_offline_data(growth_data_fn)
    dFBA_df = pd.read_csv(dFBA_fn, sep = ";")
    df = pd.read_csv(fn, sep = ";")
    # df = df.drop("Unnamed: 0", axis = 1)
    # model = cobra.io.read_sbml_model(str(MODEL_FN))
    df.set_index("Hours", inplace = True)
    # df.columns = [r.id for r in model.exchanges]

    # Remove all zero columns
    df_00 = df.copy()
    df_00[abs(df_00)<1e-6] = 0
    df_00 = df_00.loc[:, (df_00 != 0).any(axis = 0)]
    print(len(df_00.columns))

    # df = df.loc[:, df.max().sort_values(ascending = False).index]

    # df.iloc[:, :15].plot()
    # plt.show()
    # exit()

    # fig1, ax = plt.subplots(1, figsize = (20, 10))
    # Figure 1
    fig, ax = plt.subplots(1, figsize = (16, 8))
    df_1 = df[["EX_o2_e", "EX_glu__L_e", "EX_glc__D_e", "EX_pi_e", "EX_co2_e", "EX_h2o_e"]].abs()
    df_1 = df_1/df_1.max()

    df_1.columns = ["Oxygen uptake", "Glutamate uptake", "Glucose uptake", "Phosphate uptake", "Water production", "CO2 production"]
    df_1.plot(ax = ax)
    print(len(dFBA_df["Growth rate"]), len(df["EX_o2_e"]))
    df_1["Growth rate"] = list(dFBA_df["Growth rate"])
    df_1.to_csv("../Plots/fig1.csv")
    ax.set_ylabel("Normalized rates")
    dFBA_df.plot(x = "Hours", y = "Growth rate", linestyle = "--", ax = ax, secondary_y = True, alpha = 0.8, c = "k")
    ax.set_ylabel("Growth rate [/h]")
    fig.savefig("../Plots/growth.png", bbox_inches = "tight", pad_inches = 0.5)
    fig.clf()

    fig2, [ax2, ax3] = plt.subplots(1,2, figsize = (20, 8))

    df["DM_germicidin_c"] = df[["DM_germicidinA_c", "DM_germicidinB_c", "DM_germicidinC_c", "DM_germicidinD_c"]].sum(axis = 1)
    df["EX_ACT"] = df[["DM_ACT_c", "EX_gACT_e"]].sum(axis = 1)
    df["EX_RED"] = df[["DM_strprub_c", "DM_RED_c"]].sum(axis = 1)

    df_2 = df[["DM_germicidin_c", "EX_ACT", "EX_RED", "EX_CDA_e", "DM_msh_c"]]
    df_2.columns = ["Germicidin", "Actinorhodin", "Red", "CDA", "Mycothiol"]
    df_2.plot(ax = ax2)
    # ax.plot(dFBA_df["Hours"], growth_data_df["RED"], label = "Measured RED")
    ax2.set_ylabel("Secretion rate [mmol/gDW/h]")
    ax2.set_title("Secondary metabolites")
    # plt.show()
    df_2.to_csv("../Plots/fig2.csv")
    # plt.savefig("../Plots/antibiotics.png", bbox_inches = "tight", pad_inches = 0.5)
    # fig2.clf()



    # fig3, ax3 = plt.subplots(1, figsize = (16, 10))

    df_3 = df[["EX_ala__L_e", "EX_pyr_e", "EX_ac_e", "EX_gly_e", "EX_thr__L_e", "EX_nh4_e"]]
    df_3.columns = ["Alanine", "Pyruvate", "Acetate", "Glycine", "Threonine", "Ammonium"]
    df_3.plot(ax = ax3)
    ax3.set_ylabel("Secretion rate [mmol/gDW/h]")
    df_3.to_csv("../Plots/fig3.csv")
    ax3.set_title("Primary metabolites")
    # plt.savefig("../Plots/aminoacids.png", bbox_inches = "tight", pad_inches = 0.5)
    fig2.savefig("../Plots/double_plot.png", bbox_inches = "tight", pad_inches = 0.5)
    # fig3.clf()
  
if __name__ == '__main__':
    file = "exchanges_df_results_20181009_2313.csv"
    fn_exchange = Path("../Results/") / file
    dFBA_fn = Path("../Results/") / "dFBA_df_results_20181009_2313.csv"
    if 1:
        plot_exchange_results(fn_exchange, dFBA_fn)
    if 0:
        OFFLINE_DATA_FN = "../data/offline_data_F452.csv"
        plot_dFBA_results(dFBA_fn, OFFLINE_DATA_FN)