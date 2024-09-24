import warnings
warnings.filterwarnings("ignore")
import networkExpansionPy.folds as nf
import pandas as pd
from pathlib import PurePath
import concurrent.futures
from multiprocessing import Manager
import random


def fold_gated_expansion(metabolism, rn2rules, seed, algorithm, write, write_tmp, custom_write_path, ordered_outcome,
                         ignore_reaction_versions, run_counter, lock):
    for _ in range(1000):  # Run the script up to 1000 times
        RANDOMSEED = random.randint(0, 99999)
        random.seed(RANDOMSEED)
        # set file name
        str_to_append_to_fname = "shuffle2_%s" % RANDOMSEED

        # Inititalize fold metabolism with shuffle 1 (key:value pair)
        # rn2rules_shuffle = shuffle_dict_pairs(rn2rules)  # Shuffle the dictionary pairs
        # foldrules = nf.FoldRules.from_rn2rules(rn2rules_shuffle)
        # fm = nf.FoldMetabolism(metabolism, foldrules, seed)

        # Inititalize fold metabolism with shuffle 1 (key:value pair)
        rn2rules_shuffle2 = {}
        while len(rn2rules_shuffle2) != len(rn2rules):
            x2rulesNum = {}
            for rn, rules in rn2rules.items():
                for rule in rules:
                    for x in rule:
                        if x not in x2rulesNum.keys():
                            x2rulesNum[x] = 1
                        else:
                            x2rulesNum[x] += 1
            shuffleList = []
            for x in x2rulesNum:
                for i in range(int(x2rulesNum[x])):
                    shuffleList.append(x)

            rn2rules_shuffle2 = shuffle2(rn2rules, shuffleList)

        foldrules = nf.FoldRules.from_rn2rules(rn2rules_shuffle2)
        fm = nf.FoldMetabolism(metabolism, foldrules, seed)


        result = fm.rule_order(
            algorithm=algorithm,
            write=write,
            write_tmp=write_tmp,
            path=custom_write_path,
            str_to_append_to_fname=str_to_append_to_fname,
            ordered_outcome=ordered_outcome,
            ignore_reaction_versions=ignore_reaction_versions
        )
        with lock:  # Ensure atomic increment
            run_counter.value += 1
            print(f'Run count: {run_counter.value}')
        # Add additional logic if needed
    return result

def shuffle_dict_pairs(d):
    # Extract keys and values
    keys = list(d.keys())
    values = list(d.values())

    # Shuffle the keys
    random.shuffle(keys)

    # Create a new dictionary with shuffled key-value pairs
    shuffled_dict = dict(zip(keys, values))

    return shuffled_dict


def newRule(num_folds, new_rules, shuffleList):
    new_rule = random.choices(shuffleList, k=num_folds)

    # prevent duplicate folds
    new_rule = list(set(new_rule))
    while len(new_rule) < num_folds:
        new_sample = random.choices(shuffleList, k=num_folds - len(new_rule))
        new_rule.extend(new_sample)
        new_rule = list(set(new_rule))

    if frozenset(new_rule) in new_rules:  # it's a duplicate rule => run it again
        new_rule = newRule(num_folds, new_rules, shuffleList)
    return new_rule


def shuffle2(rn2rules, shuffleList):
    rn2rules_shuffle2 = {}
    for rn, rules in rn2rules.items():
        new_rules = set()
        while len(new_rules) < len(rules):
            for rule in rules:
                num_folds = len(rule)

                # check for impossible cases => redo the process
                if len(set(shuffleList)) < num_folds:
                    print('Ran out of uniques folds! Restarting shuffle process.')
                    return []

                new_rule = newRule(num_folds, new_rules, shuffleList)

                for x in new_rule:
                    shuffleList.remove(x)
                new_rules.add(frozenset(new_rule))

        rn2rules_shuffle2[rn] = new_rules
    return rn2rules_shuffle2


if __name__ == '__main__':
    preNONE = pd.read_pickle(
        '../fold_pre-expansion-results.26Sep2023/2023-09-26_19-25-59_rules_to_rn_randomSeed_42_preExpansion_NONE.pkl.gz')
    preATP = pd.read_pickle(
        '../fold_pre-expansion-results.26Sep2023/2023-09-26_20-37-47_rules_to_rn_randomSeed_42_preExpansion_C00002_iter20.pkl.gz')
    print(len(set(preATP.cpds_subiter[0])))

    asset_path = nf.asset_path

    # for vanilla
    METABOLISM_PATH = PurePath(asset_path, "metabolic_networks",
                               "metabolism.v8.02Jun2024_peroxide_fix.pkl")  # path to metabolism object pickle
    RN2RULES_PATH = PurePath(asset_path, "rn2fold", "rn2rules.20230224.pkl")  # path to rn2rules object pickle
    SEED_CPDS_PATH = PurePath(asset_path, "compounds", "seeds.Goldford2022.csv")  # path to seed compounds csv

    # for FOLD-GATED
    ALGORITHM = "no_look_ahead_rules"
    WRITE = True  # write result to disk
    WRITE_TMP = False  # write after each iteration
    CUSTOM_WRITE_PATH = None  # if writing result, custom path to write to
    ORDERED_OUTCOME = False  # ignore random seed and always choose folds based on sort order
    IGNORE_REACTION_VERSIONS = True  # when maximizing for reactions, don't count versioned reactions

    ## Metabolism
    metabolism = pd.read_pickle(METABOLISM_PATH)

    ## FoldRules
    rn2rules = pd.read_pickle(RN2RULES_PATH)

    ## Modify seeds with AA and GATP_rns
    aa_cids = set(["C00037",
                   "C00041",
                   "C00065",
                   "C00188",
                   "C00183",
                   "C00407",
                   "C00123",
                   "C00148",
                   "C00049",
                   "C00025"])

    GATP_rns = {'R00200_gATP_v1',
                'R00200_gATP_v2',
                'R00430_gGTP_v1',
                'R00430_gGTP_v2',
                'R01523_gATP_v1',
                'R04144_gATP_v1',
                'R04208_gATP',
                'R04463_gATP',
                'R04591_gATP_v1',
                'R06836_gATP',
                'R06974_gATP',
                'R06975_gATP_v1'}

    # Seed (modify with pre_expansion seed)
    seed = nf.Params(
        rns = set(metabolism.network["rn"]) - set(rn2rules) | GATP_rns,
        cpds = set(preATP.cpds_subiter[0]) | aa_cids,
        folds = set(['spontaneous'])
    )

    ## Seed (no pre-expansion)
    # seed = nf.Params(
    #     rns=set(metabolism.network["rn"]) - set(rn2rules) | GATP_rns,  # start with non-fold-gated reactions
    #     cpds=set((pd.read_csv(SEED_CPDS_PATH)["ID"])) | aa_cids,
    #     folds=set(['spontaneous'])
    # )


    # run FOLD-GATED expansion up to 1000 runs, using multiprocessing
    # Counter for the number of runs
    manager = Manager()
    run_counter = manager.Value('i', 0)  # Shared counter
    lock = manager.Lock()

    with concurrent.futures.ProcessPoolExecutor(max_workers=3) as executor:
        futures = [executor.submit(fold_gated_expansion, metabolism, rn2rules, seed, ALGORITHM, WRITE, WRITE_TMP, CUSTOM_WRITE_PATH,
                                   ORDERED_OUTCOME, IGNORE_REACTION_VERSIONS, run_counter, lock)
                   for _ in range(3)]

        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
                # Result processing if needed
            except Exception as exc:
                print(f'Generated an exception: {exc}')