import sys
import pyodbc
import pandas as pd
#import matplotlib.pyplot as plt

#####################################################################################################################
#   Usage: python3 calo_production_event_number.py <dataset tag> [start run number] [end run number]                #
#                                                                                                                   #
#   dataset tag is required (ex. ana430_2024p007)                                                                   #
#   start/end run numbers are optional if no start/end run numbers are provided all runs in the dataset are used    #
#   for use with calorimeter datasets only - does not currently query tracking DSTs
#####################################################################################################################

runs = [47289,47293,47306,47310,47313,47314,47315,47316,47333,47334,47352,47357,47360,47363,47376,47377,47378,47379,47380,47381,47382,47391,47393,47394,47395,47396,47440,47454,47457,
47458,47459,47461,47463,47464,47474,47478,47480,47481,47482,47484,47485,47489,47491,47492,47494,47495,47497,47501,47502,47503,47506,47507,47509,47513,47514,47516,47517,47519,47521,
47522,47524,47525,47540,47547,47549,47552,47557,47568,47582,47601,47628,47634,47635,47636,47637,47638,47650,47657,47658,47659,47661,47662,47664,47666,47667,47714,47715,47716,47719,
47720,47722,47723,47724,47725,47726,47727,47730,47732,47733,47783,47791,47810,47824,47831,47834,47836,47852,47856,47869,47873,47889,47890,47892,47893,47940,47944,47946,47972,47975,
47976,47977,47978,47979,47982,47999,48001,48002,48008,48011,48073,48079,48082,48084,48085,48086,48089,48097,48099,48166,48181,48185,48230,48233,48236,48237,48239,48240,48244,48245,
48246,48253,48255,48256,48257,48258,48259,48260,48261,48262,48263,48265,48290,48291,48293,48294,48295,48307,48312,48313,48318,48320,48323,48325,48326,48327,48336,48337,48341,48342,
48343,48346,48347,48348,48349,48352,48356,48358,48359,48366,48367,48368,48369,48398,48409,48410,48412,48416,48417,48418,48421,48423,48455,48456,48461,48462,48469,48471,48510,48511,
48512,48513,48514,48516,48518,48519,48520,48521,48522,48529,48530,48531,48532,48533,48635,48636,48645,48656,48657,48658,48660,48701,48720,48721,48722,48725,48726,48727,48731,48732,
48734,48742,48743,48745,48746,48801,48803,48805,48806,48807,48810,48813,48824,48826,48828,48829,48832,48836,48838,48839,48863,48864,48865,48867,48883,48885,48891,48892,48893,48894,
48895,48897,48899,48900,48901,48902,48903,48935,48938,48940,48943,48946,48949,48951,48971,48982,48983,48984,48985,48986,48987,48988,48990,48991,49023,49026,49027,49029,49030,49031,
49035,49038,49042,49045,49048,49053,49054,49060,49061,49062,49063,49065,49066,49067,49069,49070,49071,49072,49073,49098,49099,49100,49111,49113,49116,49118,49121,49125,49127,49128,
49133,49136,49138,49140,49141,49143,49148,49157,49159,49165,49217,49218,49219,49224,49226,49227,49228,49229,49230,49233,49240,49241,49244,49245,49247,49248,49249,49250,49251,49254,
49262,49263,49264,49265,49266,49267,49268,49270,49307,49308,49309,49310,49311,49312,49313,49314,49316,49317,49324,49329,49330,49331,49337,49339,49340,49341,49345,49346,49347,49348,
49349,49351,49352,49356,49358,49359,49360,49361,49362,49365,49366,49367,49368,49369,49374,49376,49377,49378,49379,49380,49381,49382,49383,49384,49385,49390,49432,49433,49434,49435,
49437,49438,49439,49440,49445,49446,49447,49448,49449,49451,49452,49453,49454,49455,49456,49457,49464,49466,49468,49652,49653,49655,49656,49658,49659,49660,49661,49662,49664,49736,
49737,49742,49743,49749,49751,49752,49756,49758,49760,49761,49762,49919,49922,49923,49925,49926,49928,49930,49934,49959,49962,49964,49966,49967,49968,49969,49970,50033,50034,50036,
50045,50046,50047,50068,50070,50073,50076,50253,50262,50283,50330,50343,50356,50376,50463,50464,50465,50466,50467,50468,50470,50471,50483,50491,50492,50496,50497,50500,50501,50520,
50522,50526,50527,50528,50529,50531,50535,50536,50546,50548,50549,50550,50551,50554,50561,50563,50564,50569,50599,50600,50601,50602,50603,50605,50606,50607,50613,50615,50650,50655,
50658,50660,50661,50662,50663,50664,50666,50668,50670,50671,50673,50674,50676,50677,50853,50857,50874,50875,50885,50886,50887,50889,50897,50962,50965,50969,50981,50987,50992,50996,
51012,51015,51093,51095,51152,51154,51160,51161,51162,51163,51166,51167,51171,51172,51176,51177,51179,51183,51188,51200,51201,51202,51208,51209,51210,51211,51212,51213,51214,51215,
51220,51222,51225,51229,51230,51234,51236,51237]

def get_dst_query(tag, start_run, end_run):
    query = f"SELECT runnumber, dsttype, SUM(events) FROM datasets WHERE tag = '{tag}' AND dsttype LIKE '%DST%' "
    if start_run.strip() and end_run.strip():
        query += f"AND runnumber >= {start_run} AND runnumber <= {end_run} ";
    query += "GROUP BY runnumber, dsttype ORDER BY runnumber, dsttype;"
    return query

def get_rundb_query(start_run, end_run):
    query = f"""
    SELECT run.runnumber, run.eventsinrun 
    FROM run 
    JOIN event_numbers ON run.runnumber = event_numbers.runnumber
    WHERE run.runtype = 'physics' AND run.runnumber >= {start_run} AND run.runnumber <= {end_run} 
    AND event_numbers.hostname IN ('seb00','seb01','seb02','seb03','seb04','seb05','seb06','seb07','seb08','seb09','seb10','seb11','seb12','seb13','seb14','seb15','seb16','seb17') 
    GROUP BY run.runnumber, run.eventsinrun 
    HAVING count(DISTINCT event_numbers.hostname) = 18
    ORDER BY run.runnumber;
    """
    return query

def get_data(cursor, query):
    cursor.execute(query)
    data = [row for row in cursor.fetchall()]
    return data

def main():

    # require arguments 
    if len(sys.argv) < 2:
        print("Usage: python3 calo_production_event_number.py <dataset tag> [start run number] [end run number]")
        sys.exit(1)

    # take in arguments 
    tag = sys.argv[1]
    start_run = sys.argv[2] if len(sys.argv) > 2 else ""
    end_run = sys.argv[3] if len(sys.argv) > 3 else ""


    print("start run and end run", start_run, end_run)

    # connect to file catalog database and extract DST information
    file_catalog_conn = pyodbc.connect("DSN=FileCatalog;UID=phnxrc;READONLY=True")
    file_catalog_cursor = file_catalog_conn.cursor()
    file_catalog_query = get_dst_query(tag, start_run, end_run)
    file_catalog_data = get_data(file_catalog_cursor, file_catalog_query)
    if len(file_catalog_data) > 0: 
        print("Found entries in the File Catalog")
    file_catalog_conn.close()

    if len(file_catalog_data) > 0:
        start_run = file_catalog_data[0][0]
        end_run = file_catalog_data[-1][0]
    else:
        print("No entries found, exiting...")
        sys.exit(1)

    print("start run and end run", start_run, end_run)

    # connect to daq database and extract run information
    rundb_conn = pyodbc.connect("DSN=daq;READONLY=True")
    rundb_cursor = rundb_conn.cursor()
    rundb_query = get_rundb_query(start_run, end_run)
    rundb_data = get_data(rundb_cursor, rundb_query)
    if len(rundb_data) > 0: 
        print("Found entries in the RunDB")
    rundb_conn.close()

    # transform into dataframes and merge data
    file_catalog_data = [tuple(row) for row in file_catalog_data] 
    df1 = pd.DataFrame(file_catalog_data, columns=['runnumber', 'dsttype', 'sum_events'])
    pivot_df = df1.pivot_table(index='runnumber', columns='dsttype', values='sum_events', aggfunc='sum')
    pivot_df = pivot_df.rename(columns={
        'DST_CALOFITTING': 'dst_calofitting',
        'DST_TRIGGERED_EVENT': 'dst_triggered_event',
        'DST_CALO': 'dst_calo',
        'DST_JET': 'dst_jet',
        'DST_JETCALO': 'dst_jetcalo'
    })
    pivot_df = pivot_df.reset_index()

    print(pivot_df.head())
    filter_pivot_df = pivot_df[pivot_df['runnumber'].isin(runs)]
    filter_pivot_df.to_csv(f'{tag}_event_production_number.csv', index=False)

    rundb_data = [tuple(row) for row in rundb_data] 
    df2 = pd.DataFrame(rundb_data, columns=['runnumber','sum_events'])
    df = pd.merge(pivot_df, df2, on='runnumber', how='outer')
    df = df.rename(columns={'sum_events': 'rundb'})

    has_dst_triggered_event = 'dst_triggered_event' in df.columns
    has_dst_calo = 'dst_calo' in df.columns
    has_dst_calofitting = 'dst_calofitting' in df.columns
    
    if has_dst_triggered_event:
        fdf = df[(df['rundb'] > 5) & (df['dst_triggered_event'].isna()) & (df['rundb'].notna())]
        filtered_df = df[(df['rundb'] > 10000) & (df['dst_triggered_event'].notna()) & (df['rundb'].notna())]
    else:
        fdf = df[(df['rundb'] > 5) & (df['rundb'].notna())]
        filtered_df = df[(df['rundb'] > 10000) & (df['rundb'].notna())]
    allrundf = df[(df['rundb'] > 10000) & (df['rundb'].notna())]
    NumNotProduced = len(fdf)
    EventsNotProducted = "{:.3e}".format(fdf['rundb'].sum())
    print()
    print(f"Number of calo physics runs in RunDB not passed to production: {NumNotProduced}")
    print(EventsNotProducted,"events not passed to production")
    print()



    if has_dst_triggered_event:
        print("Fraction of dst_triggered_event/rundb events: {0:.3f}".format(filtered_df['dst_triggered_event'].sum()/filtered_df['rundb'].sum()))
        print("Fraction of dst_triggered_event/rundb events including not produced runs: {0:.3f}".format(allrundf['dst_triggered_event'].sum()/allrundf['rundb'].sum()))
        if has_dst_calofitting:
            print("Fraction of dst_calo_fitting/dst_triggered_event events: {0:.3f}".format(filtered_df['dst_calofitting'].sum()/filtered_df['dst_triggered_event'].sum()))
        if has_dst_calo:
            print("Fraction of dst_calo/rundb events: {0:.3f}".format(filtered_df['dst_calo'].sum()/filtered_df['rundb'].sum()))
            print("Fraction of dst_calo/dst_triggered_event events: {0:.3f}".format(filtered_df['dst_calo'].sum()/filtered_df['dst_triggered_event'].sum()))
    column_sums = filtered_df[['rundb'] + [col for col in ['dst_triggered_event', 'dst_calofitting', 'dst_calo'] if col in df.columns]].sum()
    allcolumn_sums = allrundf[['rundb'] + [col for col in ['dst_triggered_event', 'dst_calofitting', 'dst_calo'] if col in df.columns]].sum()
    long_runs = filtered_df[filtered_df['rundb'] > 10000000]   


    '''
    # first plot - number of events dropped from rundb event number through production to dst_triggered_event and eventually to dst_calo 
    plt.figure(figsize=(8, 6)) 
    ax = column_sums.plot(kind='bar', color='skyblue')
    for i, value in enumerate(column_sums):
        ax.text(i, value + 0.01 * value, f'{int(value)}', ha='center', va='bottom')
    plt.title(f'{tag} Calorimeter Events for Runs {start_run}-{end_run}')
    plt.ylabel('Events')
    plt.xticks(rotation=0, ha='center')

    # second plot - number of events dropped from rundb event number through production to dst_triggered_event and eventually to dst_calo 
    # including runs that have not started production
    plt.figure(figsize=(8, 6)) 
    ax = allcolumn_sums.plot(kind='bar', color='skyblue')
    for i, value in enumerate(allcolumn_sums):
        ax.text(i, value + 0.01 * value, f'{int(value)}', ha='center', va='bottom')
    plt.title(f'{tag} Calorimeter Events for All Calo Phyiscs Runs {start_run}-{end_run}')
    plt.ylabel('Events')
    plt.xticks(rotation=0, ha='center')

    # third plot - number of events in runs that are not being produced 
    plt.figure(figsize=(10, 6))  # Create a new figure
    plt.hist(fdf['rundb'], bins=30, edgecolor='black')
    plt.title(f'Run DB event number for runs with no production events for {tag} Runs {start_run}-{end_run}')
    plt.xlabel('Event number')
    plt.ylabel('Frequency')
    plt.grid(True)

    # fourth plot - production event number for long runs (tells us where we cap out in the production)
    plt.figure(figsize=(10, 6))
    plt.hist(long_runs['dst_triggered_event'], bins=30, edgecolor='black')
    plt.title(f'Production event number for runs with > 10M events in production tag {tag} Runs {start_run}-{end_run}')
    plt.xlabel('Prod. event number')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()
    '''
    
if __name__ == "__main__":
    main()