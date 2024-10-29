import os
import warnings
warnings.filterwarnings("ignore")
import networkExpansionPy.folds as nf
import pandas as pd
from pathlib import PurePath
import concurrent.futures
from multiprocessing import Manager
import random

NUM_RUNS = 1000  # number of simulations to run before stopping
NUM_CORES = 8  # number of cores to use for multiprocessing
PRE_EXP = "C00004"  # pre-expansion type
# ["NONE", "C00002", "C00004", "C00010", "C00016", "C00019", "Z00009", "Z00035", "Z00047", "ALL"]


def fold_gated_expansion(fm, algorithm, write, write_tmp, custom_write_path, ordered_outcome,
                         ignore_reaction_versions, run_counter, lock, max_runs):
    while True:
        with lock:  # Ensure atomic check and increment
            if run_counter.value >= max_runs:
                break
            run_counter.value += 1
            current_run = run_counter.value
            print(f'Run count: {current_run}')

        # set random seed for reproducibility
        RANDOMSEED = random.randint(0, 99999)
        random.seed(RANDOMSEED)

        # set file name
        # str_to_append_to_fname = "random_fold_ordering_%s" % RANDOMSEED
        str_to_append_to_fname = f"no_lookahead_preExpansion_{PRE_EXP}_{RANDOMSEED}"

        result = fm.rule_order(
            algorithm=algorithm,
            write=write,
            write_tmp=write_tmp,
            path=custom_write_path,
            str_to_append_to_fname=str_to_append_to_fname,
            ordered_outcome=ordered_outcome,
            ignore_reaction_versions=ignore_reaction_versions
        )

    return result


if __name__ == '__main__':
    asset_path = nf.asset_path

    # for vanilla
    METABOLISM_PATH = PurePath(asset_path, "metabolic_networks",
                               "metabolism.v8.01May2023.pkl")  # path to metabolism object pickle
    RN2RULES_PATH = PurePath(asset_path, "rn2fold", "rn2rules.20230224.pkl")  # path to rn2rules object pickle
    SEED_CPDS_PATH = PurePath(asset_path, "compounds", "seeds.Goldford2022.csv")  # path to seed compounds csv

    # for FOLD-GATED
    ALGORITHM = "no_look_ahead_rules"
    # ALGORITHM = "random_fold_order"
    WRITE = True  # write result to disk
    WRITE_TMP = False  # write after each iteration
    CUSTOM_WRITE_PATH = None  # if writing result, custom path to write to
    ORDERED_OUTCOME = False  # ignore random seed and always choose folds based on sort order
    IGNORE_REACTION_VERSIONS = True  # when maximizing for reactions, don't count versioned reactions

    ## Metabolism
    metabolism = pd.read_pickle(METABOLISM_PATH)
    # remove reactions that produce H2O2 before O2
    H2O2_rns = ['R00017_v1', 'R03532_v1', 'R09507_v1', 'R09740_v1', 'R09741_v1', 'R11522', 'R12455', 'R12454']
    condition = ((metabolism.network['rn'].isin(H2O2_rns)) & (metabolism.network['direction'] == 'reverse'))
    metabolism.network = metabolism.network[~condition]
    assert 'R00017_v1' not in list(metabolism.network['rn'])
    print(len(list(metabolism.network['rn'])))

    ## FoldRules
    rn2rules = pd.read_pickle(RN2RULES_PATH)
    foldrules = nf.FoldRules.from_rn2rules(rn2rules)

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

    # Seed
    pre = pd.read_pickle(f"../data/pre-expansion_seed_cpds/pre-expansion_seed_cpds_{PRE_EXP}.pkl")
    print(PRE_EXP, len(pre))
    seed = nf.Params(
        rns = set(metabolism.network["rn"]) - set(rn2rules) | GATP_rns,
        cpds = set(pre) | aa_cids,
        folds = set(['spontaneous'])
    )

    ## Inititalize fold metabolism
    fm = nf.FoldMetabolism(metabolism, foldrules, seed)

    # run FOLD-GATED expansion up to 1000 runs, using multiprocessing
    manager = Manager()
    run_counter = manager.Value('i', 0)  # Shared counter
    lock = manager.Lock()
    # change 'max_workers' and 'range()' to the number of cores to use
    # change the last parameter for the number of runs to perform
    with concurrent.futures.ProcessPoolExecutor(max_workers=NUM_CORES) as executor:
        futures = [executor.submit(fold_gated_expansion, fm, ALGORITHM, WRITE, WRITE_TMP, CUSTOM_WRITE_PATH,
                                   ORDERED_OUTCOME, IGNORE_REACTION_VERSIONS, run_counter, lock, NUM_RUNS)
                   for _ in range(NUM_CORES)]

        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
                # Result processing if needed
            except Exception as exc:
                print(f'Generated an exception: {exc}')