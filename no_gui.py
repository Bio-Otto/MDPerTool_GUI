from distutils.util import strtobool
from no_gui.MD_1 import *
from no_gui.MD_2 import *
from no_gui.Classic_MD import *
from no_gui.Velocity_Changer import *
from no_gui.get_positions_from_trajectory_file import *
from no_gui.Energy_decomposition_from_pdb_trajectory import *
from no_gui.response_time_creator import *
from simtk import unit
from simtk.openmm import *
import time
from no_gui.write_outputs import *
from logging.config import dictConfig

if __name__ == '__main__':
    import argparse
    import multiprocessing as mp

    ### VARIABLES
    state_file_name = 'state.xml'
    last_pdb = 'last.pdb'
    dissipated_trajectory_name = 'energy_perturbation_trajectory'
    undissipated_trajectory_name = 'without_energy_perturbation_trajectory'
    reference_traj_file_for_pos = str
    dissipation_traj_file_for_pos = str
    OUTPUT_DIRECTORY = str
    created_file_for_work = str
    OUTPUT_FOLDER_NAME = str

    ### PARSER
    parser = argparse.ArgumentParser(description='The Program applying Energy Dissipation Concept using powerfull '
                                                 'OpenMM Molecular Dynamic Toolkit, which also supports the Cuda '
                                                 'platform. Each residual energy calculation required for the concept '
                                                 'can be calculated using OpenMM''s flexible and useful infrastructure.'
                                                 'In addition, you can use the package only for energy decomposition. '
                                                 'For this, it will be sufficient to specify a XTC or a DCD file '
                                                 'in the script.')

    parser.add_argument('-p', '--topology', type=str, help='Need *.pdb file for loading trajectory file',
                        required=True)

    parser.add_argument('-pff', '--protein_ff', choices=['amber03', 'amber10', 'amber96', 'amber99sb', 'amber99sbildn',
                                                         'charmm36'], default='amber96', nargs='?', type=str,
                        help='Protein Forcefield (The program defaultly will use "amber96" forcefield)', required=False)

    parser.add_argument('-wff', '--water_ff', choices=['tip3p', 'tip5p', 'spce', 'tip4pew'], default='tip3p', nargs='?',
                        type=str, help='Water Forcefield (The program defaultly will use "tip3p" forcefield)',
                        required=False)

    parser.add_argument('-lts', '--long_md_time_step', default=2.0, nargs='?', type=float,
                        help='The aim is to obtain the equilibrium state of the protein, whose population is the '
                             'highest in the initial ensemble (The program defaultly will use "2 femtosecond")',
                        required=False)

    parser.add_argument('-ts', '--long_md_total_step', default=300000, nargs='?', type=int,
                        help='It is the total number of steps the simulation wants to run. (The program defaultly will '
                             'use "300000")', required=False)

    parser.add_argument('-nbc', '--nonbonded_cutoff', default=12.0, nargs='?', type=float, required=False,
                        help='cut-off was applied to the non-covalent interactions. (The program defaultly will use '
                             '"12 Å ") NOT: Switching distance must satisfy 0 <= r_switch < r_cutoff')

    parser.add_argument('-swd', '--switch_distance', default=10.0, type=float, required=False,
                        help='The program defaultly will use 10 Å NOT: Switching distance must satisfy '
                             '0 <= r_switch < r_cutoff')

    parser.add_argument('-wp', '--water_padding', default=15, nargs='?', type=int, required=False,
                        help='The program determining largest dimension of protein, and a cubic box of size(largest '
                             'dimension)+2*padding is used. (The program defaultly will use "15 Å")')

    # parser.add_argument('-dnx-use', '--use_device_index', choices=[True, False], default=False, nargs='?', type=bool,
    #                     help='This option can only be used with OpenCL or CUDA platform. You can also specify the gpu '
    #                          'number you want on systems with more than one GPU. NOTE: OpenCL must use only one gpu. '
    #                          '(eg: <- gpu_id 0> or <- gpu-id 0,1>)', required=False)

    parser.add_argument('-dnx-use', '--use_device_index', type=lambda x: bool(strtobool(x)), choices=[True, False],
                        default=False, nargs='?', required=False, help='This option can only be used with OpenCL or '
                                                                       'CUDA platform. You can also specify the gpu '
                                                                       'number you want on systems with more than one '
                                                                       'GPU. NOTE: OpenCL must use only one gpu. '
                                                                       '(eg: <- gpu_id 0> or <- gpu-id 0,1>)')

    parser.add_argument('-dnx', '--device_index', nargs='?', type=int, required=False,
                        help='This option can only be used with OpenCL or CUDA platform. NOTE: OpenCL must use only one'
                             ' gpu. (eg. for "OpenCL": <- gpu_id 0> and example for "CUDA": <- gpu-id 0,1>)')

    parser.add_argument('-temp', '--temperature', nargs='?', default=310, type=float, required=False,
                        help='The temperature unit is kelvin. (The program defaultly will use "310 Kelvin")')

    parser.add_argument('-plt', '--platform', nargs='?', type=str, default='OpenCL', required=False,
                        help='The program defaultly will use "OpenCL" platform')

    parser.add_argument('-precision', '--plt_precision', nargs='?', type=str, default='single', required=False,
                        help='The program defaultly will use platform precision as "single"')

    parser.add_argument('-nt', '--cpu_thread', nargs='?', type=int, default=mp.cpu_count() / 2, required=False,
                        help='If you chose "CPU" for simulation platform, the program automatically will use half of '
                             'all threads. For this issue you can specify threads number by indicating <-nt>')

    parser.add_argument('-friction', '--friction_coff', nargs='?', type=float, default=1.0, required=False,
                        help='The program defaultly will use "1.0 /picosecond"')

    # parser.add_argument('-minim', '--minimize', choices=[True, False], default=True, nargs='?', type=bool,
    #                     required=False, help='The program defaultly will minimize system for 500 steps automatically. '
    #                                          'But you can by pass the minimize with -minim False')

    parser.add_argument('-minim', '--minimize', type=lambda x: bool(strtobool(x)), choices=[True, False], default=True,
                        nargs='?', required=False, help='The program defaultly will minimize system for 500 steps '
                                                        'automatically. But you can by pass the minimize with -minim '
                                                        'False')

    parser.add_argument('-minim_step', '--minimize_step', default=500, nargs='?', type=int, required=False,
                        help='The program defaultly will minimize system for 500 steps if mimimize option is not '
                             '"False"')

    # parser.add_argument('-equ', '--equilibrate', choices=[True, False], default=True, nargs='?', type=bool,
    #                     required=False, help='The program defaultly will equilibrate system for 500 steps. But you can'
    #                                          ' by pass the equilibrate with -equ False')

    parser.add_argument('-equ', '--equilibrate', type=lambda x: bool(strtobool(x)), choices=[True, False], default=True,
                        nargs='?', required=False, help='The program defaultly will equilibrate system for 500 steps. '
                                                        'But you can by pass the equilibrate with -equ False')

    parser.add_argument('-equ-step', '--equilibrate_step', default=500, nargs='?', type=int, required=False,
                        help='The program defaultly will equilibrate system for 500 steps if equilibrate option is not'
                             ' "False"')

    parser.add_argument('-ri', '--report_interval', default=100, nargs='?', type=int, required=False,
                        help='The program defaultly will report situations every 100 steps.')

    # parser.add_argument('-wdcd', '--write_dcd', choices=[True, False], default=True, type=bool,
    #                     required=True, help='The program defaultly will use dcd reporting. But you can exchange it with'
    #                                         ' XTC file format')

    parser.add_argument('-wdcd', '--write_dcd', type=lambda x: bool(strtobool(x)), choices=[True, False], default=True,
                        required=False, help='The program defaultly will use dcd reporting. But you can exchange it '
                                             'with XTC file format')

    parser.add_argument('-dcd-per', '--dcd_period', default=100, nargs='?', type=int, required=False,
                        help='The program defaultly will report trajectories every 100 steps.')

    # parser.add_argument('-wxtc', '--write_xtc', choices=[True, False], default=False, nargs='?', type=bool,
    #                     required=False, help='You can use XTC file format for output', )

    parser.add_argument('-wxtc', '--write_xtc', type=lambda x: bool(strtobool(x)), choices=[True, False], default=False,
                        required=False, help='You can use XTC file format for output')

    parser.add_argument('-xtc_per', '--xtc_period', default=100, nargs='?', type=int, required=False,
                        help='If you use XTC write condition The program defaultly will report trajectories every 100'
                             ' steps. But you can change it by using <-xtc-per> argument')

    parser.add_argument('-o', '--output', nargs='?', default=None, type=str, required=False,
                        help='The Output Directory that you provide as string will be location of the all output files.'
                             ' (The program defaultly will use script running location for output)')

    #################################    ENERGY PERTURBATION SIMULATION ARGUMENTS    #################################

    parser.add_argument('-pert_res', '--perturbed_residues', nargs='+', required=True,
                        help='You must list the residue or residues you want to perturbed.')

    parser.add_argument('-speed_factor', '--velocity_speed_factor', type=int, required=True,
                        help='Indicate how many times you want to increase the velocity of the residue atoms you want '
                             'to perturbed.', )

    parser.add_argument('-pts', '--perturbation_total_step', default=2000, type=int, required=False,
                        help='It is the total number of steps the perturbation simulation wants to run. '
                             '(The program defaultly will use "2000")')

    parser.add_argument('-per_ts', '--perturbation_time_step', default=1.0, nargs='?', type=float, required=False,
                        help='(The program defaultly will use "1 femtosecond")')

    parser.add_argument('-per_ri', '--perturbation_report_interval', default=10, nargs='?', type=int, required=False,
                        help='The program defaultly will report perturbation simulation reports every 10 steps.')

    parsed = parser.parse_args()
    start_time = time.time()

    if parsed.output is None:
        OUTPUT_DIRECTORY = os.getcwd()
        created_file_for_work, OUTPUT_FOLDER_NAME = write_folder(OUTPUT_DIRECTORY)

    ### LOG FILE
    logging_config = lof_file_settings(file_name=OUTPUT_FOLDER_NAME)
    dictConfig(logging_config)
    api_logger = logging.getLogger('api_logger')
    batch_process_logger = logging.getLogger('batch_process_logger')

    ### LOG FILE CREATION
    api_logger.info("Topology File:\t%s" % parsed.topology)
    api_logger.info("Protein Forcefield:\t%s" % parsed.protein_ff)
    api_logger.info("Water Forcefield:\t%s" % parsed.water_ff)
    api_logger.info("First MD Time Step (fs):\t%s" % parsed.long_md_time_step)
    api_logger.info("Nonbonded Cutoff (Angstrom):\t%s" % parsed.nonbonded_cutoff)
    api_logger.info("Switching Distance (Angstrom):\t%s" % parsed.switch_distance)
    api_logger.info("Water Box Padding (Angstrom):\t%s" % parsed.water_padding)
    api_logger.info("GPU Utilization Condition (For OpenCL or CUDA platforms):\t%s" % parsed.use_device_index)
    api_logger.info("GPU Device Index:\t%s" % parsed.device_index)
    api_logger.info("First MD Total Step:\t%s" % parsed.long_md_total_step)
    api_logger.info("Temperature:\t%s" % parsed.temperature)
    api_logger.info("Platform:\t%s" % parsed.platform)
    api_logger.info("Precision:\t%s" % parsed.plt_precision)
    api_logger.info("Friction Cofficient:\t%s" % parsed.friction_coff)
    api_logger.info("Minimization Condition:\t%s" % parsed.minimize)
    api_logger.info("Maximum Minimization Step:\t%s" % parsed.minimize_step)
    api_logger.info("Number of CPU Threads:\t%s" % parsed.cpu_thread)
    api_logger.info("Equilibrate Condition:\t%s" % parsed.equilibrate)
    api_logger.info("Equilibrate Step:\t%s" % parsed.equilibrate_step)
    api_logger.info("Report Interval:\t%s" % parsed.report_interval)
    api_logger.info("Write to DCD Condition:\t%s" % parsed.write_dcd)
    api_logger.info("DCD Writing Period:\t%s" % parsed.dcd_period)
    api_logger.info("Write to XTC Condition:\t%s" % parsed.write_xtc)
    api_logger.info("XTC Writing Period:\t%s" % parsed.xtc_period)
    api_logger.info("Output Destination:\t%s" % created_file_for_work)
    api_logger.info("Perturbed Residue(s):\t%s" % parsed.perturbed_residues)
    api_logger.info("Speed Factor:\t%s" % parsed.velocity_speed_factor)
    api_logger.info("Perturbation Total Step:\t%s" % parsed.perturbation_total_step)
    api_logger.info("Perturbation Time Step:\t%s" % parsed.perturbation_time_step)
    api_logger.info("Dissipation Simulation Report Interval:\t%s" % parsed.perturbation_report_interval)

    Classic_MD_Engine(pdb_path=parsed.topology, protein_ff=parsed.protein_ff, water_ff=parsed.water_ff,
                      time_step=parsed.long_md_time_step, nonbondedCutoff=parsed.nonbonded_cutoff,
                      switching_distance=parsed.switch_distance, water_padding=parsed.water_padding,
                      Device_Index=parsed.use_device_index, Device_Index_Number=parsed.device_index,
                      total_Steps=parsed.long_md_total_step, temp=parsed.temperature, platform_name=parsed.platform,
                      precision=parsed.plt_precision, friction_cofficient=parsed.friction_coff,
                      minimize=parsed.minimize, minimize_steps=parsed.minimize_step, CPU_Threads=parsed.cpu_thread,
                      equilibrate=parsed.equilibrate, equilibration_step=parsed.equilibrate_step,
                      report_interval=parsed.report_interval, write_to_dcd=parsed.write_dcd,
                      dcd_write_period=parsed.dcd_period, write_to_xtc=parsed.write_xtc,
                      xtc_write_period=parsed.xtc_period, output_directory=created_file_for_work)

    last_pdb_file_path = os.path.join(created_file_for_work, last_pdb)
    modify_atoms = convert_res_to_atoms(last_pdb_file_path, parsed.perturbed_residues, 'CA')
    print(modify_atoms)
    state_file_path = os.path.join(created_file_for_work, state_file_name)

    name_of_changed_state_xml = change_velocity(state_file_path, parsed.velocity_speed_factor, modify_atoms)

    if not parsed.write_dcd and not parsed.write_xtc:
        parsed.write_dcd = True

    Dissipation_MD_Engine(pdb_path=last_pdb_file_path, state_file=name_of_changed_state_xml,
                          protein_ff=parsed.protein_ff, water_ff=parsed.water_ff,
                          time_step=parsed.perturbation_time_step, nonbondedCutoff=parsed.nonbonded_cutoff,
                          switching_distance=parsed.switch_distance, Device_Index=parsed.use_device_index,
                          Device_Index_Number=parsed.device_index,
                          dissipation_total_Steps=parsed.perturbation_total_step,
                          platform_name=parsed.platform, precision=parsed.plt_precision,
                          CPU_Threads=parsed.cpu_thread, report_interval=parsed.perturbation_report_interval,
                          write_to_dcd=parsed.write_dcd, dcd_write_period=1, write_to_xtc=parsed.write_xtc,
                          xtc_write_period=1, dissipated_traj_name=dissipated_trajectory_name,
                          output_directory=created_file_for_work)

    Reference_MD_Engine(pdb_path=last_pdb_file_path, state_file=state_file_path, protein_ff=parsed.protein_ff,
                        water_ff=parsed.water_ff, time_step=parsed.perturbation_time_step,
                        nonbondedCutoff=parsed.nonbonded_cutoff, switching_distance=parsed.switch_distance,
                        Device_Index=parsed.use_device_index, Device_Index_Number=parsed.device_index,
                        reference_total_Steps=parsed.perturbation_total_step, platform_name=parsed.platform,
                        precision=parsed.plt_precision, CPU_Threads=parsed.cpu_thread,
                        report_interval=parsed.perturbation_report_interval, write_to_dcd=parsed.write_dcd,
                        dcd_write_period=1, write_to_xtc=parsed.write_xtc, xtc_write_period=1,
                        undissipated_traj_name=undissipated_trajectory_name, output_directory=created_file_for_work)

    if parsed.write_dcd:
        reference_traj_file_for_pos = os.path.join(created_file_for_work, undissipated_trajectory_name + '.dcd')
        dissipation_traj_file_for_pos = os.path.join(created_file_for_work, dissipated_trajectory_name + '.dcd')

    if parsed.write_xtc:
        reference_traj_file_for_pos = os.path.join(created_file_for_work, undissipated_trajectory_name + '.xtc')
        dissipation_traj_file_for_pos = os.path.join(created_file_for_work, dissipated_trajectory_name + '.xtc')

    position_list, unwrap_pdb = get_openmm_pos_from_traj(last_pdb_file_path, reference_traj_file_for_pos,
                                                         dissipation_traj_file_for_pos, write_dcd=False)

    main(unwrap_pdb, position_list, 0, 250)

    getResidueResponseTimes(os.path.join(created_file_for_work, 'reference_energy_file.csv'),
                            os.path.join(created_file_for_work, 'modified_energy_file.csv'))

    print("\n--- %s seconds ---" % (time.time() - start_time))

    # python no_gui.py -p no_gui/2j0w.pdb -pff amber10 -wff tip3p -lts 3.0 -ts 3000 -nbc 10.0 -wp 10 -dnx-use True -dnx 0 -temp 300 -wdcd True

    # python no_gui.py -p no_gui/2j0w.pdb -pff amber10 -wff tip3p -lts 3.0 -ts 3000 -nbc 10.0 -wp 10 -dnx-use True -dnx 0 -temp 300 -wdcd True
