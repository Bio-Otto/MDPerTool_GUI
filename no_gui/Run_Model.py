from MD_1 import *
from MD_2 import *
from Classic_MD import *
from Velocity_Changer import *
from get_positions_from_trajectory_file import *
from Energy_decomposition_from_pdb_trajectory import *
from response_time_creator import *
from simtk import unit
from simtk.openmm import *
import time

start_time = time.time()

## CLASSIC MD PARAMETERS
pdb_file = '2j0w.pdb'  # Inıtıal pdb structure
protein_ff = 'charmm36'  # [amber03, amber10, amber96, amber99sb, amber99sbildn] === [charmm36]
water_ff = 'tip5p'  # [tip3p, tip5p, spce, tip4pew] === [spce, tip3p-pme-b, tip3p-pme-f, tip5p]
classic_md_time_step = 2.0  # femtosecond
nonbonded_cutoff = 12.0  # angstrom
water_padding = 10  # angstrom
Device_Index = False
Device_Index_Number = 1
classic_md_total_step = 3000
temperature = 310
platform_to_use = 'OpenCL'
properties = None
precision = 'single'
friction_cofficient = 1.0  # picosecond^-1
minimize = True
minimize_max_step = 500
CPU_Thread = 2
equilibrate = True
equilibration_step = 500
report_interval = 100
write_system_xml = False
system_file_name = 'system.xml'
state_file_name = 'state.xml'
last_pdb = 'last_structure.pdb'
write_to_dcd_trajectory = True
dcd_trajectory_write_period = 100
write_to_xtc_trajectory = False
xtc_trajectory_write_period = 100

## PERTURBATION DYNAMICS PARAMETERS
pertured_residue = ['SER345']
speed_factor = 4
perturbation_total_Steps = 2000
perturbation_time_step = 1.0
perturbation_report_interval = 10
dissipated_trajectory_name = 'energy_perturbation_trajectory'
undissipated_trajectory_name = 'without_energy_perturbation_trajectory'


"""
################################################  CLASSIC MD PROCEDURE  ################################################
Classic_MD_Engine(pdb_path=pdb_file, protein_ff=protein_ff, water_ff=water_ff, time_step=classic_md_time_step,
                  nonbondedCutoff=nonbonded_cutoff, water_padding=water_padding, Device_Index=Device_Index,
                  Device_Index_Number=Device_Index_Number, total_Steps=classic_md_total_step, temp=temperature,
                  platform_name=platform_to_use, properties=properties, precision=precision,
                  friction_cofficient=friction_cofficient, minimize=minimize, minimize_steps=minimize_max_step,
                  CPU_Threads=CPU_Thread, equilibrate=equilibrate, equilibration_step=equilibration_step,
                  report_interval=report_interval, write_system_xml=write_system_xml, system_file_name=system_file_name,
                  state_file_name=state_file_name, last_pdb_filename=last_pdb, write_to_dcd=write_to_dcd_trajectory,
                  dcd_write_period=dcd_trajectory_write_period, write_to_xtc=write_to_xtc_trajectory,
                  xtc_write_period=xtc_trajectory_write_period)

##############################################  DISSIPATION MD PROCEDURE  ##############################################
modify_atoms = convert_res_to_atoms(last_pdb, pertured_residue, 'CA')
print(modify_atoms)
name_of_changed_state_xml = change_velocity(state_file_name, 4, modify_atoms)

Dissipation_MD_Engine(pdb_path=last_pdb, state_file=name_of_changed_state_xml, protein_ff=protein_ff, water_ff=water_ff,
                      time_step=perturbation_time_step, nonbondedCutoff=nonbonded_cutoff, Device_Index=Device_Index,
                      Device_Index_Number=Device_Index_Number, dissipation_total_Steps=perturbation_total_Steps,
                      platform_name=platform_to_use, properties=properties, precision=precision,
                      CPU_Threads=CPU_Thread, report_interval=perturbation_report_interval,
                      write_to_dcd=write_to_dcd_trajectory, dcd_write_period=1, write_to_xtc=write_to_xtc_trajectory,
                      xtc_write_period=1, dissipated_traj_name=dissipated_trajectory_name)

Reference_MD_Engine(pdb_path=last_pdb, state_file=state_file_name, protein_ff=protein_ff, water_ff=water_ff,
                    time_step=perturbation_time_step, nonbondedCutoff=nonbonded_cutoff, Device_Index=Device_Index,
                    Device_Index_Number=Device_Index_Number, reference_total_Steps=perturbation_total_Steps,
                    platform_name=platform_to_use, properties=properties, precision=precision, CPU_Threads=CPU_Thread,
                    report_interval=perturbation_report_interval, write_to_dcd=write_to_dcd_trajectory,
                    dcd_write_period=1, write_to_xtc=write_to_xtc_trajectory, xtc_write_period=1,
                    undissipated_traj_name=undissipated_trajectory_name)


if write_to_dcd_trajectory:
    reference_traj_file_for_pos = undissipated_trajectory_name + '.dcd'
    dissipation_traj_file_for_pos = dissipated_trajectory_name + '.dcd'

if write_to_xtc_trajectory:
    reference_traj_file_for_pos = undissipated_trajectory_name + '.xtc'
    dissipation_traj_file_for_pos = dissipated_trajectory_name + '.xtc'

position_list, unwrap_pdb = get_openmm_pos_from_traj(last_pdb, reference_traj_file_for_pos,
                                                     dissipation_traj_file_for_pos, write_dcd=False)

main(unwrap_pdb, position_list, 0, 250)

getResidueResponseTimes('reference_energy_file.csv', 'modified_energy_file.csv')

print("\n--- %s seconds ---" % (time.time() - start_time))

"""

if __name__ == '__main__':
    import argparse
    """
    ## CLASSIC MD PARAMETERS
    
    classic_md_time_step = 2.0  # femtosecond
    nonbonded_cutoff = 12.0  # angstrom
    water_padding = 10  # angstrom
    Device_Index = False
    Device_Index_Number = 1
    classic_md_total_step = 3000
    temperature = 310
    
    """

    parser = argparse.ArgumentParser(description='The Program applying Energy Dissipation Concept using powerfull '
                                                 'OpenMM Molecular Dynamic Toolkit, which also supports the Cuda '
                                                 'platform. Each residual energy calculation required for the concept '
                                                 'can be calculated using OpenMM''s flexible and useful infrastructure.' 
                                                 'In addition, you can use the package only for energy decomposition. '
                                                 'For this, it will be sufficient to specify a XTC or a DCD file '
                                                 'in the script.')

    parser.add_argument('-p', '--topology', type=str, nargs='+', help='Need *.pdb file for loading trajectory file',
                        required=True)

    parser.add_argument('-pff', '--protein_forcefield', choices=['amber03', 'amber10', 'amber96', 'amber99sb',
                                                                 'amber99sbildn', 'charmm36'], nargs='?', type=str,
                        help='Protein Forcefield (The program defaultly will use "amber96" forcefield)',
                        required=False)

    parser.add_argument('-wff', '--water_forcefield', choices=['tip3p', 'tip5p', 'spce', 'tip4pew'], nargs='?',
                        type=str, help='Water Forcefield (The program defaultly will use "tip3p" forcefield)',
                        required=False)





    parsed = parser.parse_args()
    print('Result:', vars(parsed))



    print("\n--- %s seconds ---" % (time.time() - start_time))
