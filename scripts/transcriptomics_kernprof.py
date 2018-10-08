import pandas as pd
import cobra
import os.path
from matplotlib import pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
import numpy as np
from itertools import chain
import re
import time
import sympy
add = sympy.Add._from_args
mul = sympy.Mul._from_args
# import tim

from get_transcriptomics import *



MAX_GLUCOSE_UPTAKE_RATE = -1.5
# MAX_GLUCOSE_UPTAKE_RATE = -0.8
# CORRESPONDING_AMMONIUM_UPTAKE = -1.85
CORRESPONDING_AMMONIUM_UPTAKE = 1.85*MAX_GLUCOSE_UPTAKE_RATE/2.1
PI_AVAILABLE_SCALE_FACTOR = 1

REACTION_CHANGE_CONSTANT = 1

SIGNIFICANCE_THRESHOLD = 0.5


def sigmoid(x, a, b, c, d):
    y = a / (1 + np.exp(c*(x-b))) + abs(d)
    return y

def fit_PO4(df):
    df = df[~df["PO4"].isnull()]
    fit, _ = curve_fit(sigmoid, df.index, df["PO4"], p0 = [480, 40, 3, 10])
    print(fit)
    return fit



def read_offline_data(offline_data_fn = OFFLINE_DATA_FN):
    df = pd.read_csv(OFFLINE_DATA_FN, sep = ";", decimal = ",")
    df.drop("Sample #", axis = 1, inplace = True)
    df.loc[0, "CDW"] = 0
    df.loc[30, "RED":"TBP"] = 0

    df.set_index("TAI", inplace = True)
    po4_fit = fit_PO4(df)

    df.index = pd.to_timedelta(df.index, "h")
    df_resampled = df.resample("10T").interpolate(method = "linear")


    df_resampled["Fitted PO4"] = sigmoid(df_resampled.index.total_seconds()/3600, *po4_fit)
    df_resampled["Fitted PO4"] = df_resampled["Fitted PO4"]*1e-3 #Convert from mg/L to g/L


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


class DFBA(object):
    @profile
    def __init__(self, model_fn, growth_data_fn, reaction_data_fn, solver = "gurobi"):
        self.get_model(model_fn)
        self.set_solver(solver)

        self.growth_data_df = read_offline_data(growth_data_fn)
        self.reaction_data_df = pd.read_csv(reaction_data_fn, index_col = 0)
        
        self.N = len(self.growth_data_df.index)
        self.time_array_hours = self.growth_data_df.index.total_seconds()/3600
        self.time_array = self.growth_data_df.index
        self.dt = self.time_array_hours[1]-self.time_array_hours[0]
        
        self._default_dFBA_settings()

        self.make_dFBA_df()
        self.init_model()
        self.exchange_storage = np.zeros((self.N, len(self.model.exchanges)))


    @profile
    def get_model(self, model_fn = None):
        if not model_fn:
            model_fn = self.model_fn
        else:
            self.model_fn = model_fn
            
        self.model = cobra.io.read_sbml_model(model_fn)
    
    def set_solver(self, solver):
        self.model.solver = solver

    def _default_dFBA_settings(self):
        self.use_measured_PO4_uptake = False
        self.minimum_PO4_uptake = 0.02

    def set_dFBA_settings(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
    @profile
    def make_dFBA_df(self):
        """
        Make a data frame to store values for:
        - Glucose uptake in mmol/gDW/hr/L
        - Glutamate uptake in mmol/gDW/hr/L
        - PO4 uptake mmol/gDW/hr/L
        - Predicted glucose medium concentration in mmol/L
        - Predicted glutamate medium concentration in mmol/L
        - Predicted PO4 medium concentration in mmol/L
        - Biomass in gDW/L
        - Growth rate in /hr
        - Max PO4 uptake in mmol/gDW/hr/L
        - Max PO4 uptake in mmol/hr/L


        """
        dFBA_columns = ["Hours",
                        "Glucose uptake",
                        "Glutamate uptake",
                        "PO4 uptake",
                        "Glucose in medium",
                        "Glutamate in medium",
                        "PO4 in medium",
                        "Biomass",
                        "Growth rate",
                        "Max PO4 /gDW",
                        "Max PO4",
                        ]

        arr = np.zeros((self.N, 11))
        arr[:, 0] = self.time_array_hours
        self.dFBA_df = pd.DataFrame(arr)
        self.dFBA_df.columns = dFBA_columns
        # self.dFBA_df.set_index("Hours", inplace = True)
        # print(self.dFBA_df)

    def set_initial_amounts(self, glucose = None, glutamate = None, PO4 = None, biomass = 0.01):
        """
        Converts values from g/L to mmol/L and sets as initial amount
        """
        if glucose:
            self.dFBA_df["Glucose in medium"][0] = 1e3* glucose / self.model.metabolites.glc__D_e.formula_weight
        if glutamate:
            self.dFBA_df["Glutamate in medium"][0] = 1e3* glutamate / self.model.metabolites.glu__L_e.formula_weight
        if PO4:
            self.dFBA_df["PO4 in medium"][0] = 1e3* PO4 / self.model.metabolites.pi_e.formula_weight

        self.dFBA_df["Biomass"][0] = biomass

    @profile
    def init_model(self, glucose_uptake = MAX_GLUCOSE_UPTAKE_RATE):
        self.model.reactions.EX_nh4_e.lower_bound = 0
        scale_carbon_sources(self.model, ["EX_glc__D_e", "EX_glu__L_e"], glucose_uptake)
        blocked_reactions = cobra.flux_analysis.find_blocked_reactions(self.model)
        self.model.remove_reactions(blocked_reactions, remove_orphans = True)
        self.non_exchange_reaction_ids = [r.id for r in self.model.reactions if not r in self.model.exchanges]

    def calc_available_PO4(self, i):
        if self.use_measured_PO4_uptake:
            raise NotImplementedError
        else:
            self.dFBA_df["PO4 in medium"][i]
            self.dFBA_df["Biomass"][i]
            self.dFBA_df["Max PO4"][i] = self.dFBA_df["PO4 in medium"][i] / self.dt
            self.dFBA_df["Max PO4 /gDW"][i] = max(self.minimum_PO4_uptake, self.dFBA_df["Max PO4"][i]/self.dFBA_df["Biomass"][i])


    @profile
    def update_biomass_and_medium(self, i):
        if i > 0:
            self.dFBA_df["Glucose in medium"][i] =  self.dFBA_df["Glucose in medium"][i-1] - self.dFBA_df["Glucose uptake"][i-1] * self.dt * self.dFBA_df["Biomass"][i-1]
            self.dFBA_df["Glutamate in medium"][i] =  self.dFBA_df["Glutamate in medium"][i-1] - self.dFBA_df["Glutamate uptake"][i-1] * self.dt * self.dFBA_df["Biomass"][i-1]
            self.dFBA_df["PO4 in medium"][i] =  self.dFBA_df["PO4 in medium"][i-1] - self.dFBA_df["PO4 uptake"][i-1] * self.dt * self.dFBA_df["Biomass"][i-1]
            self.dFBA_df["Biomass"][i] = self.dFBA_df["Biomass"][i-1]*np.exp(self.dFBA_df["Growth rate"][i-1]*self.dt)

    def pFBA(self, solution = None):
        if solution is None:
            solution = self.model.optimize()
        with self.model as model:
            model.reactions.BIOMASS_SCO.lb = solution.x_dict["BIOMASS_SCO"] * 0.95
            model.reactions.BIOMASS_SCO.objective_coefficient = 0
            reactions_flux_exp = [(r.forward_variable, r.reverse_variable) for r in model.reactions]
            temp = chain(*reactions_flux_exp)
            pFBA_exp = add([mul((x, sympy.singleton.S.One)) for x in temp])
            objective = model.problem.Objective(pFBA_exp, direction = "min", sloppy = True)
            model.objective = objective
            solution = model.optimize()
        return solution
        # reaction_vars = chain(*((r.forward_variable, r.reverse_variable) for r in model.reactions))

    @profile
    def transFBA(self, i, solution = None):
        timepoint = str(self.time_array[i])
        if (str(timepoint) in list(self.reaction_data_df.columns)) and (solution is not None):
            solution = self._transFBA2(timepoint, solution)
        else:
            solution = self.model.optimize()
        return solution
    @profile
    def collect_exchanges(self, i, solution):
        for j, r in enumerate(self.model.exchanges):
            self.exchange_storage[i,j] = solution.x_dict[r.id]


    
    def _transFBA(self, timepoint, solution, fraction_of_optimum = 0.95):
        with self.model as model:
            current_reaction_data = self.reaction_data_df[timepoint]
            current_reaction_data = current_reaction_data[abs(current_reaction_data) > SIGNIFICANCE_THRESHOLD]
            print(timepoint, len(current_reaction_data))
            objective_list = []
            for r_id, r_change in current_reaction_data.iteritems():
                try:
                    r_current_flux = solution.x_dict[r_id] 
                except KeyError:
                    continue
                r_target_flux = r_current_flux + r_change*REACTION_CHANGE_CONSTANT
                pos_var = model.problem.Variable(r_id+"_pos_diff", lb = 0)
                neg_var = model.problem.Variable(r_id+"_neg_diff", lb = 0)
                reaction = model.reactions.get_by_id(r_id)
                cons = model.problem.Constraint(reaction.flux_expression + pos_var - neg_var - r_target_flux, lb = 0, ub = 0)
                model.add_cons_vars([pos_var, neg_var, cons])
                objective_list += [pos_var, neg_var]
            if len(objective_list):
                cobra.util.solver.fix_objective_as_constraint(model, fraction = fraction_of_optimum)    
                objective_exp = add(objective_list)
                # print(objective_exp)
                obj = model.problem.Objective(objective_exp, direction = "min")
                model.objective = obj
                solution = model.optimize(objective_sense = None)
                print(model.summary())
            else:
                print("No reaction data above threshold for time: {0}".format(timepoint))
                solution = model.optimize()

        return solution                

    @profile          
    def _transFBA2(self, timepoint, solution, fraction_of_optimum = 0.95):
        current_reaction_data = self.reaction_data_df[abs(self.reaction_data_df[timepoint]) > SIGNIFICANCE_THRESHOLD][timepoint]
        target_flux = self.get_target_flux(current_reaction_data, solution)
        objective_list = []
        print(timepoint, len(current_reaction_data), solution.x_dict["BIOMASS_SCO"])
        with self.model as model:
            for i, r_id in enumerate(self.non_exchange_reaction_ids):
                pos_var = model.problem.Variable(r_id+"_pos_diff", lb = 0)
                neg_var = model.problem.Variable(r_id+"_neg_diff", lb = 0)
                reaction = model.reactions.get_by_id(r_id)
                cons = model.problem.Constraint(reaction.flux_expression + pos_var - neg_var - target_flux[i], lb = 0, ub = 0)
                model.add_cons_vars([pos_var, neg_var, cons])
                objective_list += [pos_var, neg_var]

            cobra.util.solver.fix_objective_as_constraint(model, fraction = fraction_of_optimum)    
            objective_exp = add(objective_list)
            # print(objective_exp)
            model.objective = model.problem.Objective(objective_exp, direction = "min")
            solution = model.optimize(objective_sense = None)
            print(model.summary())
        return solution
    @profile
    def get_target_flux(self, current_reaction_data, solution):
        r_target_flux = np.zeros(len(self.non_exchange_reaction_ids))
        for i, r_id in enumerate(self.non_exchange_reaction_ids):
            try:
                flux_change = current_reaction_data[r_id]*REACTION_CHANGE_CONSTANT
            except KeyError:
                flux_change = 0
            r_target_flux[i] = solution.x_dict[r_id] + flux_change
        return r_target_flux

        
    def set_FBA_bounds(self, i):
        self.model.reactions.EX_pi_e.lower_bound = -self.dFBA_df["Max PO4 /gDW"][i]
    @profile
    def store_solution(self, solution, i):
        self.dFBA_df["Glucose uptake"][i] = - solution.x_dict["EX_glc__D_e"]
        self.dFBA_df["Glutamate uptake"][i] = - solution.x_dict["EX_glu__L_e"]
        self.dFBA_df["PO4 uptake"][i] = - solution.x_dict["EX_pi_e"]
        self.dFBA_df["Growth rate"][i] = solution.x_dict["BIOMASS_SCO"]

    def write_results_to_file(self, filename = None):
        if not filename:
            filename = "_results_{0}.csv".format(time.strftime("%Y%m%d_%H%M%"))
        self.dFBA_df.to_csv("dFBA_df"+filename)
        
        exchanges_df = pd.DataFrame(self.exchange_storage)
        exchanges_df.to_csv("exchanges_df"+filename)


    def run_dFBA(self):
        for i in range(self.N):
            self.update_biomass_and_medium(i)
            self.calc_available_PO4(i)
            self.set_FBA_bounds(i)
            solution = self.model.optimize()
            self.store_solution(solution, i)

    def run_pdFBA(self):
        solution = None
        for i in range(self.N):
            print(i)
            self.update_biomass_and_medium(i)
            self.calc_available_PO4(i)
            self.set_FBA_bounds(i)
            solution = self.pFBA(solution)
            self.collect_exchanges(i, solution)
    @profile
    def run_transdFBA(self):
        solution = None
        for i in range(self.N):
            # print(i, self.time_array[i])
            self.update_biomass_and_medium(i)
            self.calc_available_PO4(i)
            self.set_FBA_bounds(i)
            solution = self.transFBA(i, solution)
            if solution.status == "infeasible":
                print("Infeasible at timepoint {0}".format(self.time_array[i]))
                solution = self.model.optimize()
            self.collect_exchanges(i, solution)
            

    def plot_dFBA_results(self):
        fig, axes = plt.subplots(2,2, figsize = (20, 10))
        [ax1, ax2, ax3, ax4] = axes.flatten()
        self.dFBA_df.plot(x = "Hours", y = "Biomass", ax = ax1, c = "k")
        ax1.plot(self.time_array_hours, self.growth_data_df["CDW"], c = "b", lw = 4, label = "Measured CDW")

        self.dFBA_df.plot(x = "Hours", y = "Glucose uptake", ax = ax2, c = "b")
        self.dFBA_df.plot(x = "Hours", y = "Glutamate uptake", ax = ax2, c = "r")
        

        self.dFBA_df.plot(x = "Hours", y = "PO4 uptake", ax = ax3, c = "r")
        self.dFBA_df.plot(x = "Hours", y = "Max PO4 /gDW", ax = ax3, c = "b")

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
        plt.show()

    def plot_exchanges(self):
        fig, ax = plt.subplots(1, figsize = (20, 10))
        df = pd.DataFrame(self.exchange_storage)
        df.index = self.dFBA_df["Hours"]
        df.columns = [r.id for r in self.model.exchanges]

        # Remove all zero columns
        df = df.loc[:, (df != 0).any(axis = 0)]

        # Normalize each row
        df = df / df.abs().max()
        df.plot(ax = ax)
        plt.show()




def get_available_PO4_2(growth_data, timepoint, model):
    max_PO4 = growth_data.loc[timepoint, "Fitted PO4"] #g/L
    mmol_max_PO4 = 1e3*max_PO4/model.metabolites.pi_c.formula_weight #mmol/L
    return mmol_max_PO4

def get_available_PO4(growth_data, timepoint, model):
    max_PO4 = growth_data.loc[timepoint, "Max PO4"] #g/L/hr
    mmol_max_PO4 = 1e3*max_PO4/model.metabolites.pi_c.formula_weight #mmol/L/hr
    return mmol_max_PO4

    # except 
def scale_carbon_sources(model, sources, max_glucose_uptake_rate):
    total_carbon_uptake = max_glucose_uptake_rate * 6 # Glucose has 6 carbons
    carbon_bound_exp = 0
    for source_id in sources:
        # print(source_id)
        source = model.reactions.get_by_id(source_id)
        n_carbon, n_nitrogen = get_number_of_carbons_and_nitrogens_in_metabolite(list(source.metabolites)[0])
        source.lower_bound = total_carbon_uptake / n_carbon
        carbon_bound_exp += source.flux_expression * n_carbon
        
    carbon_bound = model.problem.Constraint(carbon_bound_exp, lb = total_carbon_uptake)

    model.add_cons_vars(carbon_bound)
    # model.add_cons_vars(nitrogen_bound)

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

    if 0:
        offline_df = read_offline_data()

        reaction_data_df = pd.read_csv("../data/derivatve_transcr_rxns.csv", index_col = 0)
        model = get_model()
        dFBA_trans(model, reaction_data_df, offline_df)

    if 1:
        reaction_data_fn = "../data/derivatve_transcr_rxns.csv"
        dFBA = DFBA(MODEL_FN, OFFLINE_DATA_FN, reaction_data_fn)
        dFBA.set_initial_amounts(glucose = 41, glutamate = 57, PO4 = 491.7e-3, biomass = 0.03) #g/L
        print(dFBA.reaction_data_df.abs().std().mean())
        print(dFBA.reaction_data_df.abs().mean().mean())


        dFBA.run_transdFBA()
        # dFBA.run_pdFBA()
        # dFBA.plot_dFBA_results()
        # dFBA.plot_exchanges()
        dFBA.write_results_to_file()



    # print(derivative_df.head())