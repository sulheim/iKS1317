from matplotlib import pyplot as plt
import pandas as pd
from pathlib import Path
import cobra

MODEL_FN = Path("../iKS1317.xml")


def plot_dFBA_results(dFBA_fn, growth_data_fn):
    dFBA_df = pd.read_csv(dFBA_fn, sep = ";")
    dFBA_df = pd.read_csv(growth_data_fn, sep = ";")


    fig, axes = plt.subplots(2,2, figsize = (20, 10))
    [ax1, ax2, ax3, ax4] = axes.flatten()
    dFBA_df.plot(x = "Hours", y = "Biomass", ax = ax1, c = "k")
    ax1.plot(time_array_hours, growth_data_df["CDW"], c = "b", lw = 4, label = "Measured CDW")

    dFBA_df.plot(x = "Hours", y = "Glucose uptake", ax = ax2, c = "b")
    dFBA_df.plot(x = "Hours", y = "Glutamate uptake", ax = ax2, c = "r")
    

    dFBA_df.plot(x = "Hours", y = "PO4 uptake", ax = ax3, c = "r")
    dFBA_df.plot(x = "Hours", y = "Max PO4 /gDW", ax = ax3, c = "b")

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
def plot_exchange_results(fn):
    df = pd.read_csv(fn, sep = ";")
    # df = df.drop("Unnamed: 0", axis = 1)
    # model = cobra.io.read_sbml_model(str(MODEL_FN))
    df.set_index("Hours", inplace = True)
    # df.columns = [r.id for r in model.exchanges]

    # Remove all zero columns
    # df[abs(df)<1e-6] = 0
    # df = df.loc[:, (df != 0).any(axis = 0)]

    df = df.loc[:, df.max().sort_values(ascending = False).index]

    df.iloc[:, 12:25].plot()
    plt.show()
    exit()

    # fig1, ax = plt.subplots(1, figsize = (20, 10))
    df_1 = df[["EX_o2_e", "EX_glu__L_e", "EX_glc__D_e", "EX_pi_e", "EX_co2_e", "EX_h2o_e", "EX_nh4_e"]].abs()
    df_1.plot()
    plt.show()

    df["DM_germicidin_c"] = df[["DM_germicidinA_c", "DM_germicidinB_c", "DM_germicidinC_c", "DM_germicidinD_c"]].sum(axis = 1)
    df["EX_ACT"] = df[["DM_ACT_c", "EX_gACT_e"]].sum(axis = 1)
    df["EX_RED"] = df[["DM_strprub_c", "DM_RED_c"]].sum(axis = 1)

    df_2 = df[["DM_germicidin_c", "EX_ACT", "EX_RED"]]
    df_2.plot()
    plt.show()

    df_3 = df[["EX_ala__L_e", "EX_pyr_e", "EX_nh4_e", "EX_ac_e", "EX_gly_e"]]


  
if __name__ == '__main__':
    file = "exchanges_df_results_20181009_2313.csv"
    fn_exchange = Path("../Results/") / file
    plot_exchange_results(fn_exchange)