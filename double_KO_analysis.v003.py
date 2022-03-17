import os, pickle, cobra, datetime

import multiprocessing, multiprocessing.pool
from multiprocessing import Process, Queue

###
### FUNCTIONS
###

def essential_genes_retriever(simulation_folder):

    #print(simulation_folder, end=' ')

    pairs = {}

    # 1. read the jar file
    pickle_files = os.listdir(simulation_dir + simulation_folder + '/results/')

    if len(pickle_files) == 0:
        pass

    elif len(pickle_files) == 1:
        condition_name = pickle_files[0]
        pairs[condition_name] = []
        jar = simulation_dir + simulation_folder + '/results/' + condition_name
        f = open(jar,'rb')
        [sampleID, result, double_ko_results] = pickle.load(f)
        f.close()

        # 2. add gene essentiality to each pair
        filtered = double_ko_results[double_ko_results['growth'] < original_growth_value/2]
        for pair in filtered['ids']:
            cobra_pair_tuple = tuple(pair)
            if len(cobra_pair_tuple) == 2:
                pairs[condition_name].append(cobra_pair_tuple)

    else:
        raise ValueError('\t found a diffent number of expected files')

    return pairs

def printt(message):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))

    return None

###
### MAIN
###

#
# 0. user-defined variables
#
number_of_threads = 18

simulation_dir = '/home/adrian/projects/vinland/data/hpc_results/constrained/double_KO/deployments/'
model_file = '/home/adrian/projects/vinland/data/model/Recon3DModel_301.mat'
heatmap_info_file = '/home/adrian/projects/vinland/results/heatmap.doubleKO.pickle'

#
# 1. read and simulate the model
#
#printt('load model')
#model = cobra.io.load_matlab_model(model_file)

#printt('simulate model')
#optimization_results = model.optimize()
#original_growth_value = optimization_results.objective_value
original_growth_value = 755.0032155506631
printt('biomass {}'.format(original_growth_value))

#
# 2. read BSC folders
#
simulation_folders = next(os.walk(simulation_dir))[1]
simulation_folders.sort()

simulation_folders = simulation_folders[:20]


printt('entering a parallel world of {} threads'.format(number_of_threads))
hydra = multiprocessing.pool.Pool(number_of_threads)
hydra_output = hydra.map(essential_genes_retriever, simulation_folders)
hydra.close()
printt('completed {} results'.format(len(hydra_output)))

#
# 3. find conditional essentiality
#
conditional_essentiality = {}
for result in hydra_output:
    if result != {}:
        condition_name = list(result.keys())[0]

        pairs = result[condition_name]
        for pair in pairs:
            if pair not in conditional_essentiality:
                conditional_essentiality[pair] = [condition_name]
            else:
                conditional_essentiality[pair].append(condition_name)

#
# 4. store dictionary of pairs containing conditions where they are essential
#
printt('store essential gene pairs and their conditions')
for pair in list(conditional_essentiality.keys())[:20]:
    conditional_essentiality_prob = len(conditional_essentiality[pair])/len(simulation_folders)
    print(pair, conditional_essentiality[pair][:4], len(conditional_essentiality[pair]), conditional_essentiality_prob)

f = open(heatmap_info_file,'wb')
pickle.dump(conditional_essentiality, f)
f.close()
