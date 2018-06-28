import cobra

"""

"""  
AMINO_ACID_DICT = {"glutamate": ("glu__L_c", "glutrna_c", "trnaglu_c"), 
                   "alanine": ("ala__L_c", "alatrna_c", "trnaala_c"), 
                   "arginine": ("arg__L_c", "argtrna_c", "trnaarg_c"),
                   "asparagine": ("asn__L_c", "asntrna_c", "trnaasn_c"),
                   "aspartate": ("asp__L_c", "asptrna_c", "trnaasp_c"),
                   "cysteine": ("cys__L_c", "cystrna_c", "trnacys_c"),
                   "Methionine": ("met__L_c", "mettrna_c", "trnamet_c"),
                    
                   "Glutamine": ("gln__L_c", "glntrna_c", "trnagln_c"),
                   "Glycine": ("gly_c", "glytrna_c", "trnagly_c"),
                   "Histidine": ("his__L_c", "histrna_c", "trnahis_c"),
                   "Isoleucine": ("ile__L_c", "iletrna_c", "trnaile_c"),
                   "Leucine": ("leu__L_c", "leutrna_c", "trnaleu_c"),
                   "Lysine": ("lys__L_c", "lystrna_c", "trnalys_c"),
                   "phenylalanine": ("phe__L_c", "phetrna_c", "trnaphe_c"),
                   "Proline": ("pro__L_c", "protrna_c", "trnapro_c"),
                   "Serine": ("ser__L_c", "sertrna_c", "trnaser_c"),
                   "Threonine": ("thr__L_c", "thrtrna_c", "trnathr_c"),
                   "Tryptophan": ("trp__L_c", "trptrna_c", "trnatrp_c"),
                   "Tyrosine": ("tyr__L_c", "tyrtrna_c", "trnatyr_c"),
                   "Valine": ("val__L_c", "valtrna_c", "trnaval_c"),
                    
                    # ""fmettrna_c", # N-formylmethionine
}

def get_trna_premodules(model):
    for m in model.metabolites:
        if m.id[:4] == "trna":
            yield m

def get_trna_complexes(model):
    for m in model.metabolites:
        if m.id[-6:] == "trna_c":
            yield m

def make_tRNA_biomass(model):
    new_biomass = model.reactions.BIOMASS_SCO.copy()
    new_biomass.id = "BIOMASS_SCO_tRNA"
    new_biomass.name = "Biomass function with tRNA included"
    new_biomass.objective_coefficient = 0


    trna_complexes = get_trna_complexes(model)
    trna_complexes_dict = {}
    for m in trna_complexes:
        trna_complexes_dict[m] = -1e-6
    new_biomass.add_metabolites(trna_complexes_dict)

    trna_premodules_dict = {}
    trna_premodules = get_trna_premodules(model)
    for m in trna_premodules:
        if m.id == "trnamet_c":
            trna_premodules_dict[m] = 2 * 1e-6
        else:
            trna_premodules_dict[m] = 1e-6
    new_biomass.add_metabolites(trna_premodules_dict)

    model.add_reaction(new_biomass)

def replace_aminoacids_in_biomass(model):
    new_biomass = model.reactions.BIOMASS_SCO.copy()
    new_biomass.id = "BIOMASS_SCO_tRNA"
    new_biomass.name = "S. coelicolor biomass function with AA replaced by tRNA-AA"

    summed_coefficient = 0
    for key, tup in AMINO_ACID_DICT.items():
        base_met = model.metabolites.get_by_id(tup[0])
        trna_complex_met = model.metabolites.get_by_id(tup[1])
        trna_premodule_met = model.metabolites.get_by_id(tup[2])

        # Get Biomass coefficient for that metabolite
        biomass_coeff = new_biomass.get_coefficient(base_met.id)
        summed_coefficient += biomass_coeff

        # Add trnA metabolites and remove base met
        new_biomass.add_metabolites({base_met: -biomass_coeff, trna_complex_met: biomass_coeff, trna_premodule_met: -biomass_coeff})

    # Remove excessive ATP, ADP, proton and PI and h2o nowfrom biomass
    # H2O??
    m_atp = model.metabolites.atp_c
    m_adp = model.metabolites.adp_c
    m_pi = model.metabolites.pi_c
    m_h = model.metabolites.h_c
    m_h2o = model.metabolites.h2o_c

    new_biomass.add_metabolites({m_atp: -summed_coefficient, m_h2o: -summed_coefficient, m_adp: summed_coefficient, m_pi: summed_coefficient, m_h: summed_coefficient})
    print(new_biomass.build_reaction_string())
    model.add_reaction(new_biomass)
    model.reactions.BIOMASS_SCO_tRNA.objective_coefficient = 0
    return model 



if __name__ == '__main__':
    model = cobra.io.read_sbml_model("C:/Users/snorres/git/gem_sco/iKS1317.xml")
    if 0:
        make_tRNA_biomass(model)
        model.optimize()
        model.summary()
        cobra.io.write_sbml_model(model, "C:/Users/snorres/git/gem_sco/iKS1317.xml")
    if 1:
        model = replace_aminoacids_in_biomass(model)
        model.optimize()
        model.summary()
        model.reactions.BIOMASS_SCO_tRNA.objective_coefficient = 1
        model.reactions.BIOMASS_SCO.objective_coefficient = 0
        model.summary()
        model.optimize()
        cobra.io.write_sbml_model(model, "C:/Users/snorres/git/gem_sco/iKS1317.xml")

