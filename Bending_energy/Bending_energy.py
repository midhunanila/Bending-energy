import sys
import numpy as np
import os
import MDAnalysis as mda


################################################################################
# Read command-line arguments
################################################################################

if len(sys.argv) < 4:
    print("""
Usage:
python bending_energy.py GRO_FILE XTC_FILE1 XTC_FILE2 N_GRIDS [ATOM_SELECTION]

Arguments
---------
GRO_FILE        : GROMACS topology (.gro)
XTC_FILE       : First trajectory (.xtc)
N_GRIDS         : Number of grid divisions
ATOM_SELECTION  : MDAnalysis atom selection (optional)

Example
-------
python bending_energy.py system.gro traj1.xtc traj2.xtc 8 "resname POPC CHOL DBSM POPS DPG3"
""")
    sys.exit(1)

gro_file = sys.argv[1]

xtc_files = [
    sys.argv[2]
]

N_of_grids = int(sys.argv[3])

if len(sys.argv) > 4:
    lipid_selection = sys.argv[4]
else:
    lipid_selection = "resname POPC CHOL DBSM POPS DPG3"

def safe_xtc_list(xtc_files):
    """Filter out missing or empty XTC files."""
    valid = []
    for f in xtc_files:
        if os.path.exists(f) and os.path.getsize(f) > 0:
            valid.append(f)
        else:
            print(f"Skipping invalid file: {f}")
    return valid
#define a function which calculate the bending energy of membrane
def calculate_energy_binding(gro_file,
                             xtc_files,
                             n_grids,
                             lipid_selection):

  valid_xtcs = safe_xtc_list(xtc_files)

  if not valid_xtcs:
        print(f"No usable XTC files found for {fij}")
        return [], [], [], [], []

    # Load GROMACS trajectory
  universe = mda.Universe(gro_file, *valid_xtcs)
  #select the lipid tail region
  lipid_tail = universe.select_atoms(lipid_selection)

  print('We have {} beads in the lipid_tail'.format(len(lipid_tail)))

  #list for bending energy of each frame
  Bending_energy_for_each_frame_sga = []
  Bending_energy_for_each_frame_sga2 = []
  Bending_energy_for_each_frame_gamma = []
  Bending_energy_for_each_frame_gamma2 = []
  energy_sga = []
  energy_sga2 = []
  energy_gamma = []
  energy_gamma2 = []
  times = []



  #loop for all the frames
  for ts in universe.trajectory[:5]:
      a = universe.dimensions[0] / N_of_grids
      b = universe.dimensions[1] / N_of_grids
      coord_table = lipid_tail.atoms.positions

      # Create a zero matrix with order of number of grids
      M = [[[] for _ in range(N_of_grids)] for _ in range(N_of_grids)]

      # Adjust coordinates for PBC and place them in the grid
      for coord in coord_table:
          i_x, i_y = int(coord[0] / a), int(coord[1] / b)
          if 0 <= i_x < N_of_grids and 0 <= i_y < N_of_grids:
              M[i_x][i_y].append(coord)

      # Find mean height along z axis
      M2 = np.zeros((N_of_grids, N_of_grids))

      for i in range(N_of_grids):
          for j in range(N_of_grids):
              sublist = M[i][j]
              if sublist:
                  M2[i, j] = np.mean(np.array(sublist)[:, 2])


      gradient = np.zeros((N_of_grids, N_of_grids, 2))
      gradient2 = np.zeros((N_of_grids, N_of_grids, 2))
      gradient_square = np.zeros((N_of_grids, N_of_grids))
      gradient_square2 = np.zeros((N_of_grids, N_of_grids))
      factor = np.zeros((N_of_grids, N_of_grids, 2))
      factor2 = np.zeros((N_of_grids, N_of_grids, 2))
      gamma = np.zeros((N_of_grids, N_of_grids))

      # Compute gradient
      for i in range(N_of_grids):
          for j in range(N_of_grids):
              gradient[i, j] = np.array([(M2[(i+1) %N_of_grids, j] - M2[(i) %N_of_grids, j]) / (a), (M2[i, (j+1) %N_of_grids] - M2[i, (j) %N_of_grids]) / (b)])
              gradient2[i, j] = np.array([(M2[(i+1) %N_of_grids, j] - M2[(i-1) %N_of_grids, j]) / (2*a), (M2[i, (j+1) %N_of_grids] - M2[i, (j-1) %N_of_grids]) / (2*b)])
      # Compute gradient square
      for i in range(N_of_grids):
          for j in range(N_of_grids):
              gradient_square[i, j] = np.dot(gradient[i, j], gradient[i, j])
              gradient_square2[i, j] = np.dot(gradient2[i, j], gradient2[i, j])


      # Compute factor
      for i in range(N_of_grids):
          for j in range(N_of_grids):
              factor[i, j]= gradient[i, j] / np.sqrt(1 + gradient_square[i,j])
              factor2[i, j]= gradient2[i, j] / np.sqrt(1 + gradient_square2[i,j])


      curvature_gamma = 0
      curvature_sga = 0
      curvature_sga2 = 0
      curvature_gamma2 = 0
      for i in range(N_of_grids):
        for j in range(N_of_grids):
          #find square of laplacian to find the bending energy
          curvature_gamma = curvature_gamma + (np.sqrt(1 + gradient_square[i,j]) * ((factor[(i) %N_of_grids, j, 0] - factor[(i-1) %N_of_grids, j, 0])/(a) + (factor[i, (j) %N_of_grids, 1] - factor[i, (j-1) %N_of_grids, 1])/(b))**2)
          curvature_gamma2 = curvature_gamma2 + (np.sqrt(1 + gradient_square2[i,j]) * ((factor2[(i+1) %N_of_grids, j, 0] - factor2[(i-1) %N_of_grids, j, 0])/(2*a) + (factor2[i, (j+1) %N_of_grids, 1] - factor2[i, (j-1) %N_of_grids, 1])/(2*b))**2)
          curvature_sga = curvature_sga + np.square((M2[(i+1)%N_of_grids,j] + M2[(i-1)%N_of_grids,j] - 2*M2[i,j])/((a)**2) + (M2[i,(j+1)%N_of_grids] + M2[i,(j-1)%N_of_grids] - 2*M2[i,j])/((b)**2))
          curvature_sga2 = curvature_sga2 + np.square((M2[(i+2)%N_of_grids,j] + M2[(i-2)%N_of_grids,j] - 2*M2[i,j])/((2*a)**2) + (M2[i,(j+2)%N_of_grids] + M2[i,(j-2)%N_of_grids] - 2*M2[i,j])/((2*b)**2))
      Bending_energy_sga = curvature_sga*a*b
      Bending_energy_sga2 = curvature_sga2*a*b
      Bending_energy_gamma = curvature_gamma*a*b
      Bending_energy_gamma2 = curvature_gamma2*a*b


      energy_sga.append(Bending_energy_sga)
      energy_sga2.append(Bending_energy_sga2)
      energy_gamma.append(Bending_energy_gamma)
      energy_gamma2.append(Bending_energy_gamma2)

      times.append(ts.frame)



  return energy_sga, energy_sga2, energy_gamma, energy_gamma2, times
  
if __name__ == "__main__":


    energy_sga, energy_sga2, energy_gamma, energy_gamma2, time = \
        calculate_energy_binding(
            gro_file,
            xtc_files,
            N_of_grids,
            lipid_selection
        )


    print("\nAnalysis completed successfully.")
    print(f"Frames analyzed : {len(time)}")


    print(f"Average SGA Energy      : {np.mean(energy_sga):.6f}")
    print(f"Average SGA2 Energy     : {np.mean(energy_sga2):.6f}")
    print(f"Average Gamma Energy    : {np.mean(energy_gamma):.6f}")
    print(f"Average Gamma2 Energy   : {np.mean(energy_gamma2):.6f}")



    # Save energy of every frame

    with open("energy_binding_results.txt", "w") as f:

        f.write("Frame\tSGA\tSGA2\tGamma\tGamma2\n")


        for i in range(len(time)):

            f.write(
                f"{time[i]}\t"
                f"{energy_sga[i]:.6f}\t"
                f"{energy_sga2[i]:.6f}\t"
                f"{energy_gamma[i]:.6f}\t"
                f"{energy_gamma2[i]:.6f}\n"
            )


    print("\nFrame-by-frame energies saved to energy_binding_results.txt")
