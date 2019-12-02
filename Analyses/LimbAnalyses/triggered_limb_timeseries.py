# This script loads in csv files from matlab and converts them into a pandas dataFrame and computes associated timeseries

import os
import numpy as np
import pandas as pd

#TODO: Figure out how to get this to run faster and update code accordingly
def make_timeseries(file_names):

    # Define the save_path
    save_path = '/Users/bdeangelis/Desktop/Datasets/OptoMethodDatasets/csv/trials-withBody.csv'

    # Define the path to the desired file
    parent_path = '/Users/bdeangelis/Desktop/Datasets/OptoMethodDatasets/csv'

    # Define the window length
    win_len = 15 # NOTE: Matches (DeAngelis et al. 2019 eLife)

    # Define the hit variables for later reference
    hit_list = ['L1_hit', 'L2_hit', 'L3_hit', 'R1_hit', 'R2_hit', 'R3_hit']

    # Define a trial counter
    counter = 1

    # Loop through each of the files
    for idx, val in enumerate(file_names):

        # Read the csv into a pandas dataFrame
        full_path = os.path.join(parent_path, val)
        data = pd.read_csv(full_path)

        # # TEST: Truncate the dataset
        # data = data[0:1000]

        # Add a column that is true if any limb was hit
        data['any_hit'] = (data.L1_hit| data.L2_hit | data.L3_hit | data.R1_hit| data.R2_hit | data.R3_hit)

        # Check if this is the first iteration
        if idx == 0:

            # Create an empty dataFrame for storing the output time-series
            time_cols = list(data.columns)
            time_cols.extend(['trial_id', 'condition', 'hit_type'])
            timeseries = pd.DataFrame(columns=time_cols)

        # Groupby id and get all the trajectories
        grouped = data.groupby('uniqueFlyTrajID')

        for id, group in grouped:

            # Check that the id has enough data to have a valid trajectory
            if group.shape[0] >= (win_len*2)+1:

                # Append the new data to the dataFrame
                temp = [group[int(x-win_len):int(x+win_len+1)].reset_index(drop=True) for x in np.argwhere(group.any_hit == 1)]

                # Check the length of each section that is being appended
                for v in temp:
                    if v.shape[0] == ((win_len*2) + 1):

                        # Add a variable that is the trial_id
                        v['trial_id'] = counter

                        # Add a variable for the condition
                        v['condition'] = val

                        # Add a variable that defines the hit_type for the current trial
                        # TODO: Define the conditions for each hit_type
                        if sum(v.loc[win_len,hit_list]) > 1:
                            v['hit_type'] = 'Multi'
                        elif v['L1_hit'][win_len] == 1:
                            v['hit_type'] = 'L1'
                        elif v['L2_hit'][win_len] == 1:
                            v['hit_type'] = 'L2'
                        elif v['L3_hit'][win_len] == 1:
                            v['hit_type'] = 'L3'
                        elif v['R1_hit'][win_len] == 1:
                            v['hit_type'] = 'R1'
                        elif v['R2_hit'][win_len] == 1:
                            v['hit_type'] = 'R2'
                        elif v['R3_hit'][win_len] == 1:
                            v['hit_type'] = 'R3'
                        else:
                            print('Error in code, need to check why there is no hit type')

                        # Append the current trial to the new dataset
                        timeseries = timeseries.append(v)

                        # Increment the counter
                        counter = counter + 1

        # Print a status update
        print(f'Dataset {idx+1} of {len(file_names)} complete.')

    # Convert the output into a useable format for downstream analyses
    val_list = ['L1_xPlot_mm', 'L2_xPlot_mm', 'L3_xPlot_mm',
                'R1_xPlot_mm', 'R2_xPlot_mm', 'R3_xPlot_mm',
                'L1_yPlot_mm', 'L2_yPlot_mm', 'L3_yPlot_mm',
                'R1_yPlot_mm', 'R2_yPlot_mm','R3_yPlot_mm',
                'angVel_radPerSec', 'forwardSpeed_mmPerSec',
                'translationalSpeed_mmPerSec']
    # val_list = ['L1_xPlot_mm', 'L2_xPlot_mm', 'L3_xPlot_mm',
    #             'R1_xPlot_mm', 'R2_xPlot_mm', 'R3_xPlot_mm',
    #             'L1_yPlot_mm', 'L2_yPlot_mm', 'L3_yPlot_mm',
    #             'R1_yPlot_mm', 'R2_yPlot_mm','R3_yPlot_mm']
    timeseries['time'] = timeseries.index
    trials = timeseries.pivot_table(index=['trial_id','condition', 'hit_type'], columns='time', values=val_list)

    # Write the trials dataFrame to disk
    trials.to_csv(save_path, index=True)

if __name__ == '__main__':

    # Define the file_name
    file_names = ['Bristle-Chrimson.csv',
                  'Control.csv',
                  'Gr5a-Chrimson.csv']

    # Call the function above
    make_timeseries(file_names)