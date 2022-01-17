from FullAnalysis_functions import *
from parameter_and_initial_values import *

from scipy.optimize import *
    
def CreateReactionFunction(information):
    full_labels = ["D", "Cd", "Nd", "m", "Cm", "Nm", "Ym", "YC", "P"]
    mRNA_labels = ["D", "Cd", "Nd", "m"]
    t_max = 100000
    n = 100000
    
    from sympy.utilities.lambdify import lambdify, implemented_function
    
    transcription_rates_syms = symbols("km_p, km_m, km, km_c")
    translation_rates_syms = symbols("kd_p, kd_m, kd, kd_c")
    miBin_rate_syms = symbols("kbind")
    ExoPasBin_rates_syms = symbols("kX_m")
    ExoDeg_rates_syms = symbols("kdeg")
    miPasDil_rates_syms = symbols("lam_mi")
    ProPasDil_rates_syms = symbols("lam_P")
    mi3Rec_rate_syms = symbols("kX_3")
    mi5Rec_rate_syms = symbols("kX_5")


    full_rates_syms = list(transcription_rates_syms+translation_rates_syms)+[ExoPasBin_rates_syms]+[ExoDeg_rates_syms]+[ProPasDil_rates_syms]
    mRNA_rates_syms = list(transcription_rates_syms)+[miPasDil_rates_syms]
    mi3Deg_rates_syms = [miBin_rate_syms] + [mi3Rec_rate_syms] + [ExoDeg_rates_syms] 
    mi5Deg_rates_syms = [miBin_rate_syms] + [mi5Rec_rate_syms] + [ExoDeg_rates_syms] 
    all_rates_syms = list(transcription_rates_syms+translation_rates_syms)+[miBin_rate_syms]+[ExoPasBin_rates_syms]+[mi3Rec_rate_syms]+[ExoDeg_rates_syms]+[ProPasDil_rates_syms]+[miPasDil_rates_syms]

    species = {}
    species_syms = {}
    species_syms_sens = {}
    ODEs = {}
    rates_syms = {}
    resources_syms = list(symbols("Pol, Rib, Exo"))
    total_species_syms = []

    species["Pol"], species["Rib"], species["Exo"] = 0, 0, 0
    species_syms["Pol"], species_syms["Rib"], species_syms["Exo"] = resources_syms[0], resources_syms[1], resources_syms[2]
    species_syms_sens["Pol"], species_syms_sens["Rib"], species_syms_sens["Exo"] = symbols("Pol"), symbols("Rib"), symbols("Exo")
    ODEs["Pol"], ODEs["Rib"], ODEs["Exo"] = 0, 0, 0

    # Iterate through each species 
    for i in range(len(information)):

        #The information for the current species 
        info = information[i]

        #If current species has a full production pathway, DNA->RNA->Protein
        if info[1] == "full":
            #List for all sub-species
            temp_species_syms = []
            #Iterate through all sub-species
            for j in range(len(full_labels)):
                species[full_labels[j]+"_"+info[0]]=0
                species_syms[full_labels[j]+"_"+info[0]] = symbols(full_labels[j]+"_"+info[0])
                species_syms_sens["S"+full_labels[j]+"_"+info[0]] = symbols("S"+full_labels[j]+"_"+info[0])
                ODEs[full_labels[j]+"_"+info[0]]=0
                temp_species_syms = temp_species_syms + [symbols(full_labels[j]+"_"+info[0])]
            IN = temp_species_syms + resources_syms + full_rates_syms
            out = production_full_pathway(*IN) #Re name function
            for j in range(len(full_labels)):
                ODEs[full_labels[j]+"_"+info[0]] = ODEs[full_labels[j]+"_"+info[0]] + out[j]
            ODEs["Pol"] = ODEs["Pol"] + out[-3]
            ODEs["Rib"] = ODEs["Rib"] + out[-2]
            ODEs["Exo"] = ODEs["Exo"] + out[-1]
            for j in range(len(full_rates_syms)):
                if full_rates_syms[j] not in rates_syms:
                    rates_syms[str(full_rates_syms[j])] = full_rates_syms[j]

        #If current species has a mRNA production pathway, DNA->RNA
        if info[1] == "mRNA":
            #List for all sub-species
            temp_species_syms = []
            #Iterate through all sub-species
            for j in range(len(mRNA_labels)):
                species[mRNA_labels[j]+"_"+info[0]]=0
                species_syms[mRNA_labels[j]+"_"+info[0]] = symbols(mRNA_labels[j]+"_"+info[0])
                species_syms_sens["S"+mRNA_labels[j]+"_"+info[0]] = symbols("S"+mRNA_labels[j]+"_"+info[0])
                ODEs[mRNA_labels[j]+"_"+info[0]]=0
                temp_species_syms = temp_species_syms + [symbols(mRNA_labels[j]+"_"+info[0])]
            IN = temp_species_syms + resources_syms + mRNA_rates_syms
            out = production_mRNA_pathway(*IN) #Re name function
            for j in range(len(mRNA_labels)):
                ODEs[mRNA_labels[j]+"_"+info[0]] = ODEs[mRNA_labels[j]+"_"+info[0]] + out[j]
            ODEs["Pol"] = ODEs["Pol"] + out[-1]
            for j in range(len(mRNA_rates_syms)):
                if mRNA_rates_syms[j] not in rates_syms:
                    rates_syms[str(mRNA_rates_syms[j])] = mRNA_rates_syms[j]

    for i in range(len(information)):

        #The information for the current species 
        info = information[i]

        #If current species is regulated by 3'miRNA
        if info[2] == "1TS3":
            mi3_labels = [["m_"+info[0], "m_3", "Exo", "Rib", "Lm"+info[2]+"_"+info[0], "Ym"+info[2]+"_"+info[0], "Nm_"+info[0], "Cm_"+info[0]], ["Cm_"+info[0], "m_3", "Exo", "Rib", "LC"+info[2]+"_"+info[0], "YC"+info[2]+"_"+info[0], "Nm_"+info[0], "m_"+info[0]]]
            #Iterate through all sub-species
            for k in range(len(mi3_labels)):
                temp_species_syms = []
                for j in range(len(mi3_labels[k])):
                    if mi3_labels[k][j] in species:
                        temp_species_syms = temp_species_syms + [species_syms[mi3_labels[k][j]]]
                    else:
                        species[mi3_labels[k][j]]=0
                        species_syms[mi3_labels[k][j]] = symbols(mi3_labels[k][j])
                        species_syms_sens["S"+mi3_labels[k][j]] = symbols("S"+mi3_labels[k][j])
                        ODEs[mi3_labels[k][j]]=0
                        temp_species_syms = temp_species_syms + [symbols(mi3_labels[k][j])]
                IN = temp_species_syms + mi3Deg_rates_syms
                out =  mRNA_3miRNAdegradation_fullpathway(*tuple(list(IN)+[k]+[1]))
                for j in range(len(mi3_labels[k])):
                    ODEs[mi3_labels[k][j]] = ODEs[mi3_labels[k][j]] + out[j]
                for j in range(len(mi3Deg_rates_syms)):
                    if mi3Deg_rates_syms[j] not in rates_syms:
                        rates_syms[str(mi3Deg_rates_syms[j])] = mi3Deg_rates_syms[j]


        if info[2] == "1TS5":
            mi5_labels = ["m_"+info[0], "m_5", "Exo", "Rib", "Lm"+info[2]+"_"+info[0], "Ym"+info[2]+"_"+info[0], "Nm_"+info[0], "K_"+info[0], "Cm_"+info[0]]
            #Iterate through all sub-species
            temp_species_syms = []
            for j in range(len(mi5_labels)):
                if mi5_labels[j] in species:
                    temp_species_syms = temp_species_syms + [species_syms[mi5_labels[j]]]
                else:
                    species[mi5_labels[j]]=0
                    species_syms[mi5_labels[j]] = symbols(mi5_labels[j])
                    species_syms_sens["S"+mi5_labels[j]] = symbols("S"+mi5_labels[j])
                    ODEs[mi5_labels[j]]=0
                    temp_species_syms = temp_species_syms + [symbols(mi5_labels[j])]
            IN = temp_species_syms + mi5Deg_rates_syms
            out =  mRNA_5miRNAdegradation_fullpathway(*tuple(list(IN)+[0]+[1]))
            for j in range(len(mi5_labels)):
                ODEs[mi5_labels[j]] = ODEs[mi5_labels[j]] + out[j]
            for j in range(len(mi5Deg_rates_syms)):
                if mi5Deg_rates_syms[j] not in rates_syms:
                    rates_syms[str(mi5Deg_rates_syms[j])] = mi5Deg_rates_syms[j]
                    
                    
        if info[2] == "3TS3":
            mi3_labels = [["m_"+info[0], "m_3", "Exo", "Rib", "Lm"+info[2]+"_"+info[0], "Ym"+info[2]+"_"+info[0], "Nm_"+info[0], "Cm_"+info[0]], ["Cm_"+info[0], "m_3", "Exo", "Rib", "LC"+info[2]+"_"+info[0], "YC"+info[2]+"_"+info[0], "Nm_"+info[0], "m_"+info[0]]]
            #Iterate through all sub-species
            for k in range(len(mi3_labels)):
                temp_species_syms = []
                for j in range(len(mi3_labels[k])):
                    if mi3_labels[k][j] in species:
                        temp_species_syms = temp_species_syms + [species_syms[mi3_labels[k][j]]]
                    else:
                        species[mi3_labels[k][j]]=0
                        species_syms[mi3_labels[k][j]] = symbols(mi3_labels[k][j])
                        species_syms_sens["S"+mi3_labels[k][j]] = symbols("S"+mi3_labels[k][j])
                        ODEs[mi3_labels[k][j]]=0
                        temp_species_syms = temp_species_syms + [symbols(mi3_labels[k][j])]
                IN = temp_species_syms + mi3Deg_rates_syms
                out =  mRNA_3miRNAdegradation_fullpathway(*tuple(list(IN)+[k]+[3]))
                for j in range(len(mi3_labels[k])):
                    ODEs[mi3_labels[k][j]] = ODEs[mi3_labels[k][j]] + out[j]
                for j in range(len(mi3Deg_rates_syms)):
                    if mi3Deg_rates_syms[j] not in rates_syms:
                        rates_syms[str(mi3Deg_rates_syms[j])] = mi3Deg_rates_syms[j]


        if info[2] == "3TS5":
            mi5_labels = ["m_"+info[0], "m_5", "Exo", "Rib", "Lm"+info[2]+"_"+info[0], "Ym"+info[2]+"_"+info[0], "Nm_"+info[0], "K_"+info[0], "Cm_"+info[0]]
            #Iterate through all sub-species
            print("mi5_labels",mi5_labels)
            temp_species_syms = []
            for j in range(len(mi5_labels)):
                if mi5_labels[j] in species:
                    temp_species_syms = temp_species_syms + [species_syms[mi5_labels[j]]]
                else:
                    species[mi5_labels[j]]=0
                    species_syms[mi5_labels[j]] = symbols(mi5_labels[j])
                    species_syms_sens["S"+mi5_labels[j]] = symbols("S"+mi5_labels[j])
                    ODEs[mi5_labels[j]]=0
                    temp_species_syms = temp_species_syms + [symbols(mi5_labels[j])]
            IN = temp_species_syms + mi5Deg_rates_syms
            print(IN)
            print(tuple(list(IN)+[0]+[3]))
            out =  mRNA_5miRNAdegradation_fullpathway(*tuple(list(IN)+[0]+[3]))
            for j in range(len(mi5_labels)):
                ODEs[mi5_labels[j]] = ODEs[mi5_labels[j]] + out[j]
            for j in range(len(mi5Deg_rates_syms)):
                if mi5Deg_rates_syms[j] not in rates_syms:
                    rates_syms[str(mi5Deg_rates_syms[j])] = mi5Deg_rates_syms[j]



    ODEs_list = list(ODEs.values())
    ODEs_array = np.array(ODEs_list)
    reaction_function = lambdify(list(species_syms.values())  + list(rates_syms.values()), ODEs_array)     
    def model(y, t, reaction_function, args):
        return reaction_function(*y, *args)
    Initial_values = [Initial_values_dict[i] for i in species_syms.keys()]
    Rates_values = [Rates_values_dict[i] for i in rates_syms.keys()]
    
    Rates_indexes = []
    for i in range(len(list(rates_syms.keys()))):
        Rates_indexes = np.concatenate((Rates_indexes,[list(Rates_values_dict.keys()).index(list(rates_syms.keys())[i])]))
    
    return {"Model": model, "Reaction_functions": reaction_function, "Initial_values": Initial_values,"Rates_values": Rates_values, "Rates_indexes":Rates_indexes, "Protien_indexes": [ list(species_syms.keys()).index("P_G"), list(species_syms.keys()).index("P_R")], "Rates_string": list(rates_syms.keys())}