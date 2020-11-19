# Energy Dissipation Concept with MDPERTOOL v0.1

A Software Tool for Investigation of Allosteric Communication within Protein Structures via Energy Dissipation in Molecular Dynamics Simulations

For terminal usage just use **_python no_gui.py -p pdb_file -wdcd True -pert_res 'SER345' -speed_factor 4_**

* **_-p_**    --> Give the absolute path of your pdb file. 
* **_-wdcd_**    --> The program defaultly will use dcd reporting. But you can exchange it with XTC file format using **_-wdcd False -wxtc True_** argument
* **_-pert_res_**   --> You must list the residue or residues you want to perturbed.
* **_-speed_factor_**   --> Indicate how many times you want to increase the velocity of the residue atoms you want to perturbed.


## Then run it.

### Also you can check full functional parameters with typing <**_python no_gui.py -h_**>


The Program applying Energy Dissipation Concept using powerfull OpenMM Molecular Dynamic Toolkit, which also supports the Cuda platform. Each residual energy calculation required for the concept can be calculated using OpenMMs flexible and useful infrastructure.In addition, you can use the package only for energy decomposition. For this, it will be sufficient to specify a XTC or a DCD file in the script.

optional arguments:
  -h, --help            show this help message and exit
  -p TOPOLOGY, --topology TOPOLOGY
                        Need .pdb file for loading trajectory file
  -pff [{amber03,amber10,amber96,amber99sb,amber99sbildn,charmm36}]
                        Protein Forcefield (The program defaultly will use
                        "amber96" forcefield)
  -wff [{tip3p,tip5p,spce,tip4pew}]
                        Water Forcefield (The program defaultly will use
                        "tip3p" forcefield)
  -lts [LONG_MD_TIME_STEP], --long_md_time_step [LONG_MD_TIME_STEP]
                        The aim is to obtain the equilibrium state of the
                        protein, whose population is the highest in the
                        initial ensemble (The program defaultly will use "2
                        femtosecond")
  -ts [LONG_MD_TOTAL_STEP], --long_md_total_step [LONG_MD_TOTAL_STEP]
                        It is the total number of steps the simulation wants
                        to run. (The program defaultly will use "300000")
  -nbc [NONBONDED_CUTOFF], --nonbonded_cutoff [NONBONDED_CUTOFF]
                        cut-off was applied to the non-covalent interactions.
                        (The program defaultly will use "12 Å ") NOT:
                        Switching distance must satisfy 0 <= r_switch <
                        r_cutoff
  -wp [WATER_PADDING], --water_padding [WATER_PADDING]
                        The program determining largest dimension of protein,
                        and a cubic box of size(largest dimension)+2.padding
                        is used. (The program defaultly will use "15 Å")
  -dnx-use [{True,False}]
                        This option can only be used with OpenCL or CUDA
                        platform. You can also specify the gpu number you want
                        on systems with more than one GPU. NOTE: OpenCL must
                        use only one gpu. (eg: <- gpu_id 0> or <- gpu-id 0,1>)
  -dnx [DEVICE_INDEX], --device_index [DEVICE_INDEX]
                        This option can only be used with OpenCL or CUDA
                        platform. NOTE: OpenCL must use only one gpu. (eg. for
                        "OpenCL": <- gpu_id 0> and example for "CUDA": <- gpu-
                        id 0,1>)
  -temp [TEMPERATURE], --temperature [TEMPERATURE]
                        The temperature unit is kelvin. (The program defaultly
                        will use "310 Kelvin")
  -plt [PLATFORM], --platform [PLATFORM]
                        The program defaultly will use "OpenCL" platform
  -precision [PLT_PRECISION], --plt_precision [PLT_PRECISION]
                        The program defaultly will use platform precision as
                        "single"
  -nt [CPU_THREAD], --cpu_thread [CPU_THREAD]
                        If you chose "CPU" for simulation platform, the
                        program automatically will use half of all threads.
                        For this issue you can specify threads number by
                        indicating <-nt>
  -friction [FRICTION_COFF], --friction_coff [FRICTION_COFF]
                        The program defaultly will use "1.0 /picosecond"
  -minim [{True,False}], --minimize [{True,False}]
                        The program defaultly will minimize system for 500
                        steps automatically. But you can by pass the minimize
                        with -minim False
  -minim_step [MINIMIZE_STEP], --minimize_step [MINIMIZE_STEP]
                        The program defaultly will minimize system for 500
                        steps if mimimize option is not "False"
  -equ [{True,False}], --equilibrate [{True,False}]
                        The program defaultly will equilibrate system for 500
                        steps. But you can by pass the equilibrate with -equ
                        False
  -equ-step [EQUILIBRATE_STEP], --equilibrate_step [EQUILIBRATE_STEP]
                        The program defaultly will equilibrate system for 500
                        steps if equilibrate option is not "False"
  -ri [REPORT_INTERVAL], --report_interval [REPORT_INTERVAL]
                        The program defaultly will report situations every 100
                        steps.
  -wdcd [{True,False}], --write_dcd [{True,False}]
                        The program defaultly will use dcd reporting. But you
                        can exchange it with XTC file format
  -dcd-per [DCD_PERIOD], --dcd_period [DCD_PERIOD]
                        The program defaultly will report trajectories every
                        100 steps.
  -wxtc [{True,False}], --write_xtc [{True,False}]
                        You can use XTC file format for output
  -xtc_per [XTC_PERIOD], --xtc_period [XTC_PERIOD]
                        If you use XTC write condition The program defaultly
                        will report trajectories every 100 steps. But you can
                        change it by using <-xtc-per> argument
  -pert_res PERTURBED_RESIDUES, --perturbed_residues PERTURBED_RESIDUES
                        You must list the residue or residues you want to
                        perturbed.
  -speed_factor VELOCITY_SPEED_FACTOR, --velocity_speed_factor VELOCITY_SPEED_FACTOR
                        Indicate how many times you want to increase the
                        velocity of the residue atoms you want to perturbed.
  -pts PERTURBATION_TOTAL_STEP, --perturbation_total_step PERTURBATION_TOTAL_STEP
                        It is the total number of steps the perturbation
                        simulation wants to run. (The program defaultly will
                        use "2000")
  -per_ts [PERTURBATION_TIME_STEP], --perturbation_time_step [PERTURBATION_TIME_STEP]
                        (The program defaultly will use "1 femtosecond")
  -per_ri [PERTURBATION_REPORT_INTERVAL], --perturbation_report_interval [PERTURBATION_REPORT_INTERVAL]
                        The program defaultly will report perturbation
                        simulation reports every 10 steps.
