import pandas as pd
import math


def _set_empty_idx_dictionary(df):
    """
    Helper function that takes in the concatenated output csv from `bootstrap_LD.py`
    and returns a dictionary of "comparison_name;individual_ID" keys and empty list
    values. This dictionary will be later loaded with indeces to filter the dataframe by.

    :param df: dataframe from concatenated output csv from `bootstrap_LD.py`
    :type df: pandas.DataFrame

    :returns: dictionary that stores all possible comparision;individual pairs as keys with [] values
    :rtype: dict
    """

    #setting empty dictiory
    empty_dict = {}

    #finding unique comparision names
    comp_names = set(df["comparison_name"].to_list())

    #finding unique individual names
    all_indvs_str = ":".join(df["Samples"].to_list())
    indv_names = set(all_indvs_str.split(":"))

    #finding all possible pairs of comparisions and indidivuals
    for comp in comp_names:
        for indv in indv_names:
            empty_dict[f"{comp};{indv}"] = []

    return empty_dict


def _find_indeces(df, idx_dictionary):
    """
    Helper function runs through the dataframe and fills an empty idx dictionary with row
    indeces that contain the matching comparison and individual names

    :param df: dataframe from concatenated output csv from `bootstrap_LD.py`
    :param idx_dictionary: empty dictionary from `_set_empty_idx_dictionary`

    :type df: pandas.DataFrame
    :type idx_dictionary: dict

    :returns: filled in dictionary with values of lists of row indeces for each matching key
    :rtype: dict
    """

    #iterate through dataframe
    for row in df.itertuples():
        row_idx = row.Index
        row_comp_name = row.comparison_name
        row_samples_str = row.Samples

        #iterate through samples
        row_samples = row_samples_str.split(":")
        for sample in row_samples:

            #update idx dictionary
            idx_dictionary[f"{row_comp_name};{sample}"].append(row_idx)

    return idx_dictionary





def build_per_individual_LDstat_table(df):
    """
    Function builds a table of per individual LDstat values by filtering the larger dataframe
    by comparison & individual names and their corresponding row indeces and averaging the
    remaining Gmean_Dprime values. Thus, each individual's LDstat value is the average of the 
    Dprime value of the boostraps it appears in

    :param df: dataframe from concatenated output csv from `bootstrap_LD.py`

    :type df: pandas.DataFrame

    :returns: dataframe of per individual LDstat values
    :rtype: pandas.DataFrame
    """

    #setting up dictionary of row indeces for filtering
    empty_dictionary = _set_empty_idx_dictionary(df)
    idx_dictionary = _find_indeces(df, empty_dictionary)

    #setting up mapping dictionary from comparison name to comparison type
    comp_type_map = dict(zip(df["comparison_name"], df["comparison_type"]))

    #empty lists for new dataframe
    comparison_names = []
    sample_names = []
    ldstats_g = []
    ldstats_m = []
    ldstats_mnz = []

    #iterating through dictionary of indeces and computing per individual LDstat
    for k, v in idx_dictionary.items():

        #filtering df
        filtered_df = df.iloc[v]

        #computing LDstat
        ld_stat_g = filtered_df["dprime_gmean"].mean()
        ld_stat_m = filtered_df["dprime_mean"].mean()
        ld_stat_mnz = filtered_df["dprime_mean_nonzero"].mean()

        #appending data
        comparison_names.append(k.split(";")[0])
        sample_names.append(k.split(";")[1])
        ldstats_g.append(ld_stat_g)
        ldstats_m.append(ld_stat_m)
        ldstats_mnz.append(ld_stat_mnz)

    #constructing df
    outdf = pd.DataFrame({"comparison_name":comparison_names, "Sample":sample_names, "LDstat_gmean":ldstats_g, "LDstat_mean":ldstats_m, "LDstat_mean_nonzero":ldstats_mnz})

    #adding comparison type column
    outdf["comparison_type"] = outdf["comparison_name"].map(comp_type_map)

    #reorder columns
    outdf = outdf[["comparison_name", "comparison_type", "Sample", "LDstat_gmean", "LDstat_mean", "LDstat_mean_nonzero"]]

    return outdf





def mainfunc(bootstrapLD_csv, outLDstat_csv):
    """
    Main Function that takes in a bootstrapLD output csv and saves a per individual LDstat csv

    :param bootstrapLD_csv: concatenated output csv from `bootstrap_LD.py`
    :param outLDstat_csv: name for output csv

    :type bootstrapLD_csv: str
    :type outLDstat_csv: str
    """

    bootstrapLD_df = pd.read_csv(bootstrapLD_csv)
    LDstat_df = build_per_individual_LDstat_table(bootstrapLD_df)
    LDstat_df.to_csv(outLDstat_csv, index=False)







mainfunc("bootstrapLD_all_comparisions.csv", "individual_bootstapped_LDstats.csv")
