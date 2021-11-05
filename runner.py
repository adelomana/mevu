import cobra
import pandas, numpy, sys, multiprocessing, datetime, pickle

###
### FUNCTIONS
###

def printt(message):

    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S \t {}".format(message)))

    return None

def boundary_setter(reaction, reaction_weights, unique_reaction_set):

    # f.1. define presence of logic rules
    grr = reaction.gene_reaction_rule
    code = 0
    if 'or' in grr:
        code = code + 1
    if 'and' in grr:
        code = code + 2

    # f.2. procceed to change

    # if there are no rules
    if code == 0:
        w = reaction_weights.values[0]
        nlb = reaction.bounds[0]*w
        nub = reaction.bounds[1]*w
        new_reaction_bounds = (nlb, nub)

    # there are only OR rules
    elif code == 1: # the mean of all reaction genes. if no expression available, 1
        working_w = list(reaction_weights.values)
        missing_elements = len(unique_reaction_set) - len(working_w)
        for element in range(missing_elements):
            working_w.append(1)
        w = numpy.mean(working_w)
        nlb = reaction.bounds[0]*w
        nub = reaction.bounds[1]*w
        new_reaction_bounds = (nlb, nub)

    # there are only AND rules
    elif code == 2:
        w = numpy.min(reaction_weights.values)
        nlb = reaction.bounds[0]*w
        nub = reaction.bounds[1]*w
        new_reaction_bounds = (nlb, nub)

    # there are both AND and OR rules
    elif code == 3:
        all_cases = []
        or_cases = grr.split(' or ')
        isoforms = []
        for case in or_cases:
            if 'and' in case:
                and_cases = case.split(' and ')
                and_cases_cleanish = [element.replace('(', '') for element in and_cases]
                and_cases_cleaned = [element.replace(')', '') for element in and_cases_cleanish]
                and_cases_bounds = []
                for element in and_cases_cleaned:
                    if element in reaction_weights:
                        and_cases_bounds.append(reaction_weights[element])

                if len(and_cases_bounds) != 0:
                    z = numpy.min(and_cases_bounds)
                    all_cases.append(z)
            else:
                isoform_root = case.split('.')[0]
                if case in reaction_weights:
                    all_cases.append(reaction_weights[case])
                else:
                    if isoform_root not in isoforms:
                        all_cases.append(1)
                # adding at the level of isoform_root =
                isoforms.append(isoform_root)

        w = numpy.mean(all_cases)
        nlb = reaction.bounds[0]*w
        nub = reaction.bounds[1]*w
        new_reaction_bounds = (nlb, nub)
    else:
        raise ValueError('logic rule not anticipated')

    return new_reaction_bounds

def constrain_evaluator(task):

    sampleID = task[0]
    model = task[1]
    flux_weights = task[2]

    #! iterate reactions
    with model as model:
        for reaction in model.reactions:

            # f.1. define a new boundary for each reaction
            local_reaction_genes = [gene.id for gene in reaction.genes]
            unique_reaction_set = list(set([element.split('.')[0] for element in local_reaction_genes]))

            # f.1.1. check if the reaction contains genes and if so, if they are ammenable to constraining
            if len(local_reaction_genes) >= 1:
                ammenable_genes = [element for element in local_reaction_genes if element in flux_weights.index]
                if len(ammenable_genes) >= 1:

                    # f.1.2. actual constraining
                    reaction_weights = flux_weights[ammenable_genes]
                    new_reaction_bounds = boundary_setter(reaction, reaction_weights, unique_reaction_set)
                    reaction.bounds = new_reaction_bounds

        #! f.2. perform analysis
        model.optimize()
        result = model.objective.value
        double_ko_results = cobra.flux_analysis.double_gene_deletion(model, model.genes, processes=threads)

    #! f.3. store results
    jar_file = 'results/{}.pickle'.format(sampleID)
    f = open(jar_file,'wb')
    pickle.dump([sampleID, result, double_ko_results], f)

    return result

###
### MAIN
###

if __name__ == '__main__':

    printt('WELCOME')

    ###
    ### 0. user defined variables
    ###

    flux_weights_file = '/home/hpce17/hpce17000/projects/hpc/data/expression/formatted_flux_change_results_file.csv'
    model_file = '/home/hpce17/hpce17000/projects/hpc/data/Recon3DModel_301.mat'

    threads = 48

    cobra_config = cobra.Configuration()
    cobra_config.processes = threads
    print(cobra.Configuration())

    ###
    ### 1. read information
    ###

    printt('read information')

    ### 1.1. read boundary weights
    printt('read boundary expression weights')
    weights = pandas.read_csv(flux_weights_file, index_col=0)
    print(weights.shape)
    print(weights.head())

    ### 1.2. read model
    printt('read model')
    model = cobra.io.load_matlab_model(model_file)

    ###
    ### 2. constraining
    ###

    printt('constraining')

    sampleID = weights.index[30]
    printt('sample ID {}'.format(sampleID))
    task = [sampleID, model, weights.loc[sampleID]]
    result = constrain_evaluator(task)
    printt('completed job for one sample with growth value {}'.format(result))

    ###
    ### 4. farewell
    ###
    printt('goodbye')
