import os

# Base directory path containing the .list files
base_directory = '/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mbd_eff_list_files'

# Get all files in the directory with the .list extension
file_list = [f for f in os.listdir(base_directory) if f.endswith('.list')]

# Separate files into two categories based on their names
calo_files = sorted([f for f in file_list if 'dst_calofitting' in f])
track_files = sorted([f for f in file_list if 'dst_trkr_tracks' in f])
track_cluster_files = sorted([f for f in file_list if 'dst_trkr_cluster' in f])
track_seed_files = sorted([f for f in file_list if 'dst_trkr_seed' in f])

# Number of files per output file
files_per_runlist = 10000

for i in range(0, len(track_files), files_per_runlist):
    chuck = track_files[i:i + files_per_runlist]

    calo_runlist_filename = f'dst_calofitting_runlist_{i // files_per_runlist}.txt'
    track_runlist_filename = f'dst_trkr_tracks_runlist_{i // files_per_runlist}.txt'
    track_cluster_filename = f'dst_trkr_cluster_runlist_{i // files_per_runlist}.txt'
    track_seed_filename = f'dst_trkr_seed_runlist_{i // files_per_runlist}.txt'

    with open(track_runlist_filename, "w") as track_runlist_file, \
         open(calo_runlist_filename, "w") as calo_runlist_file, \
         open(track_cluster_filename, "w") as track_cluster_runlist_file, \
         open(track_seed_filename, "w") as track_seed_runlist_file:

        for file in chuck:
            track_runlist_file.write(os.path.join(base_directory, file) + '\n')

            parts = file.replace('.list', '').split('-')
            run_number = parts[-2]
            segment_number = parts[-1]

            track_cluster_file = f'dst_trkr_cluster-{run_number}-{segment_number}.list'
            track_cluster_runlist_file.write(os.path.join(base_directory, track_cluster_file) + '\n')
            track_seed_file = f'dst_trkr_seed-{run_number}-{segment_number}.list'
            track_seed_runlist_file.write(os.path.join(base_directory, track_seed_file) + '\n')

            # Convert to integer, divide by 10, and reformat with leading zeros
            segment_int = int(segment_number)
            calo_segment = f"{segment_int // 10:05d}"

            # Construct the calo filename
            calo_file = f"dst_calofitting-{run_number}-{calo_segment}.list"
            calo_runlist_file.write(os.path.join(base_directory, calo_file) + '\n')

'''

# Function to create runlist files
def create_runlist_files(file_list, category):
    for i in range(0, len(file_list), files_per_runlist):
        chunk = file_list[i:i + files_per_runlist]
        
        # Create the filename for the runlist file
        runlist_filename = f'{category}_runlist_{i // files_per_runlist}.txt'
        
        # Write the chunk of file paths to the runlist file
        with open(runlist_filename, 'w') as runlist_file:
            for file in chunk:
                runlist_file.write(os.path.join(base_directory, file) + '\n')

# Create runlist files for each category
create_runlist_files(calo_files, 'dst_calofitting')
create_runlist_files(track_files, 'dst_trkr_tracks')

print("Runlist files created successfully!")

'''
