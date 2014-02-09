import h5py
from glob import glob



def loadTrack(rep_num, h5_path_regx):
    h5_files = sorted(glob(h5_path_regx))

    energy_track = [3, 1, 4, 1, 5, 9, 2, 6, 9, 10]  # used to compare to remove the duplicates

    track = []

    for h5_file in h5_files:
        h = h5py.File(h5_file, 'r')
        h_dset = h['dset']

        rep = h_dset[rep_num][0]  
        # e.g. query into 3rd step
        # step_num = 2
        # step = rep[step_num]
        # idx_rep, idx_tmp, idx_lig, idx_prt = step[0]
        # energies = step[1][0]
        # mcc = step[1][1]
        # move_vector = step[2]

        for i in rep:
            if (energy_track[-1] == i[1][0]).all():
                pass
            else:
                energy_track.append(i[1][0])
                track.append((i[3], i[1][0][9], i[1][1][0], i[2]))  # step, total_energy, mcc, move_vector
                
    return track
    
if __name__ == "__main__":
    rep_num = 0
    h5_path_regx = "out*/*.h5"
    track = loadTrack(rep_num, h5_path_regx)

    print "energy,mcc"
    for i in track:
        print "%f,%f" % (i[1], i[2])



