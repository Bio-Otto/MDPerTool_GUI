import os.path

import yaml
import io


def config_template():
    template = {
        'Inputs': {
            'topology': '',
            'output directory': '',
            'run duration': 0.002,
            'run duration unit': 'nanosecond',
            'threading number': 2,
            'all cpu is active': False,
            'r factor': '3,4',
            'selected residue': list,
            'perturbation run duration': 1000,
            'perturbation run duration unit': 'step',
            'perturbation threading number': 2,
            'perturbation all cpu is active': False,
        },
        'Simulation': {
            'equilibrium platform': 'CUDA',
            'equilibrium precision': 'single',
            'equilibrium device id': 1,
            'equilibrium use device id': False,
            'perturbation platform': 'CUDA',
            'perturbation precision': 'double',
            'perturbation device id': 1,
            'perturbation use device id': False,
            'protein forcefield': 'amber03',
            'water forcefield': 'tip3p',
            'water geometry padding': '10 * angstrom',
            'equilibrium integrator': 'Langevin',
            'equilibrium time step': 2.0,
            'equilibrium time step unit': 'femtosecond',
            'equilibrium additional integrator parameters': True,
            'friction': '91.0/picosecond',
            'temperature': '310.0 * kelvin',
            'nonbonded method': 'PME',
            'constraints': None,
            'rigid water is active': True,
            'nonbonded cutoff': '1.2 * nanometer',
            'switching distance': '1 * nanometer',
            'use switching distance': True,
            'number of simulation steps': True,
            'minimize': True,
            'minimization max iterations': 500,
            'Equilibrate': True,
            'Equilibration steps': 500,
            'DCD report': '1.2 * nanometer',
            'XTC report': '1 * nanometer',
            'state data report': True,
            'DCD writing frequency': 100,
            'DCD output name': 'output.dcd',
            'XTC writing frequency': 100,
            'XTC output name': 'output.xtc',
            'state data frequency': 100,

        }
    }
    return template


def read_output_configuration_file(file_path):
    try:
        with open(file_path, 'r') as stream:
            data_loaded = yaml.safe_load(stream)
            return data_loaded
    except yaml.YAMLError as exc:
        print(exc)


def write_output_configuration_file(file_path, template_yml):
    if file_path.split('.')[-1] != 'yaml':
        file_path = file_path+'.yaml'

    with io.open(file_path, 'w', encoding='utf8') as outfile:
        yaml.dump(template_yml, outfile, default_flow_style=False, allow_unicode=True, sort_keys=True)
