from FullAnalysis_functions import *
from parameter_and_initial_values import *

from scipy.optimize import *

def FullAnalysisRun(information, information_label):
    import numpy as np
    import scipy.integrate
    from sympy import symbols, Matrix
    import math
    from sympy.utilities.lambdify import lambdify, implemented_function
    import matplotlib.pyplot as plt
    from parameter_and_initial_values import time_points

    full_labels = ["D", "Cd", "Nd", "m", "Cm", "Nm", "Ym", "YC", "P"]
    mRNA_labels = ["D", "Cd", "Nd", "m"]
    t_max = 100000
    n = 100000
    
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
    mi5Deg_rates_syms = [miBin_rate_syms] + [mi3Rec_rate_syms] + [ExoDeg_rates_syms] 
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


    t = np.linspace(0, t_max, n)    
    simulation_results = scipy.integrate.odeint(model, Initial_values, t, args=(reaction_function, Rates_values))
    
    print("---------------------------------------------------------")
    print("|    Species \t|    Initial values \t|  Final values\t|")
    print("---------------------------------------------------------")
    for i in range(len(species_syms)):
        print("|",list(species_syms.keys())[i],"\t|" if len(str(list(species_syms.keys())[i]))>=5 else "\t\t|", round(simulation_results[0][i],3),"\t\t|" if len(str(round(simulation_results[0][i],3)))>=5 else "\t\t\t|" , round(simulation_results[-1][i],3),"\t|" if len(str(round(simulation_results[-1][i],3)))>=5 else "\t\t|")
    print("---------------------------------------------------------")
    print("|     Rates\t|  Rate values\t|")
    print("---------------------------------")
    for i in range(len(rates_syms)):
        print("|", list(rates_syms.keys())[i], "\t|" if len(str(list(rates_syms.keys())[i]))>=5 else "\t\t|" , Rates_values[i],"\t|" if len(str(Rates_values[i]))>=5 else "\t\t|")
    print("---------------------------------")
    
    
    for i in range(len(ODEs_array)):
        print("d/dt(",list(species_syms.keys())[i], ") = ", ODEs_array[i])
    
    num_rates = len(list(rates_syms.keys()))
    num_species = len(list(species_syms.keys()))
    Heatmap = np.zeros([num_rates,num_species])
    def model_sensitivity(y, t, sensitivity_functions, parameter):
        return np.squeeze(sensitivity_functions(*y))
    print("Sentivities calculated for: ")
    for i in range(num_rates):
        f = Matrix(ODEs_list)
        y = list(species_syms.values())
        parameter = list(rates_syms.values())[i]
        J = f.jacobian(y)
        F = Matrix(f.diff(parameter))
        dS_dt = J*Matrix(list(species_syms_sens.values())) + F
        sensitivity_functions = lambdify(list(species_syms.values())+list(rates_syms.keys()), dS_dt)
        sensitivity_functions = sensitivity_functions(*simulation_results[:][-1], *Rates_values)
        sensitivity_functions = lambdify(list(species_syms_sens.values()), sensitivity_functions)
        sensitivity_results = scipy.integrate.odeint(model_sensitivity, np.zeros(len(simulation_results[:][-1])), t, args=(sensitivity_functions, parameter))
        Heatmap[i,:] = sensitivity_results[-1,:]*Rates_values[i]/(simulation_results[-1,:] + (1*10**-50))
        print( i+1, " of ", num_rates, "rates")
    ylabels = list(rates_syms.keys())
    xlabels = list(species_syms.keys())

    fig, ax = plt.subplots(figsize=(35, 10))
    ax.set_title("Heatmap 1: Sensitivities", fontsize=30)
    max_pos_val = abs(np.max(Heatmap))
    max_neg_val = abs(np.min(Heatmap))
    if max_pos_val<max_neg_val:
        max_val = max_neg_val
    else:
        max_val = max_pos_val
    img = ax.imshow(Heatmap, cmap='RdBu', vmax=max_val, vmin=-max_val)
    ax.set_yticks([i for i in range(num_rates)])
    ax.set_yticklabels(ylabels, fontsize=20)
    ax.set_xticks([i for i in range(num_species)])
    plt.xticks(rotation=90)
    ax.set_xticklabels(xlabels, fontsize=20)
    fig.colorbar(img)

    for i in range(num_species):
        max_pos_val = abs(max(Heatmap[:,i]))
        max_neg_val = abs(min(Heatmap[:,i]))
    #     print("max_pos_val", max_pos_val,"max_neg_val", max_neg_val)
        if max_pos_val<max_neg_val:
            max_val = max_neg_val
        else:
            max_val = max_pos_val
    #     print("heatmap:",Heatmap[:,i])
    #     print("max_val:",max_val)
        Heatmap[:,i] = np.array([Heatmap[:,i]])/(max_val+ (1*(10**-50)))
    #     print(Heatmap[:,i])
    #     print(max_val)

    fig, ax = plt.subplots(figsize=(35, 10))
    ax.set_title("Heatmap 2: Normalised Sensitivities", fontsize=30)
    max_pos_val = abs(np.max(Heatmap))
    max_neg_val = abs(np.min(Heatmap))
    if max_pos_val<max_neg_val:
        max_val = max_neg_val
    img = ax.imshow(Heatmap, cmap='RdBu', vmax=max_val, vmin=-max_val)
    ax.set_yticks([i for i in range(num_rates)])
    ax.set_yticklabels(ylabels, fontsize=20)
    ax.set_xticks([i for i in range(num_species)])
    plt.xticks(rotation=90)
    ax.set_xticklabels(xlabels, fontsize=20)
    fig.colorbar(img)

    
    experimental_data_indecies = []
    experimental_estimation_data = []
    for i in range(len(experimental_data[information_label])):
        experimental_data_index = list(species_syms.keys()).index(list(experimental_data[information_label].keys())[i])
        experimental_data_indecies = experimental_data_indecies + [experimental_data_index]
        experimental_estimation_data = np.concatenate((experimental_estimation_data, [list(experimental_data[information_label].values())[i]][0]))

    rates_string = ""
    for i in range(len(list(rates_syms.keys()))):
        rates_string = rates_string + list(rates_syms.keys())[i]+","
        
    exec('''def model_predict(t,'''+rates_string+'''):\n
            \t simulation_results_parameter_estimation = scipy.integrate.odeint(model, Initial_values, t, args=(reaction_function, Rates_values))\n        
            \t output_data = []\n
            \t for i in range(len(experimental_data_indecies)):\n
            \t\t output_data = np.concatenate((output_data, simulation_results[time_points,experimental_data_indecies[i]]))\n
            \t return output_data \n''', locals(), globals())

    bounds = np.zeros([2,len(list(rates_syms.keys()))])
    for i in range(len(list(rates_syms.keys()))):
        bounds[:,i] = Rates_bounds_dict[list(rates_syms.keys())[i]]

    popt, pcov = curve_fit(model_predict, t, experimental_estimation_data, bounds=bounds)
    print("-----------------------------------------------------------------")
    print("|      Rates\t|     Initial Values\t|   Estimated Values\t|")
    print("-----------------------------------------------------------------")
    for i in range(len(popt)):
        initial_rate_value = round(Rates_values[i], 3 - int(math.floor(math.log10(abs(Rates_values[i])))) - 1)
        final_rate_value = round(popt[i], 3 - int(math.floor(math.log10(abs(popt[i])))) - 1)
        print("|",list(rates_syms.keys())[i],"\t\t|" if len(list(rates_syms.keys())[i])<=4 else "\t|" , initial_rate_value,"\t\t\t|" if len(str(initial_rate_value))<=4 else "\t\t|" , final_rate_value,"\t\t\t|" if len(str(final_rate_value))<=4 else "\t\t|")
    print("-----------------------------------------------------------------")
    
    return