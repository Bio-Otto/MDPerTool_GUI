#!/usr/bin/env python
from __future__ import division, print_function

import sys

# OpenMM Imports
import openmm as mm
import openmm.app as app
from collections import namedtuple
# ParmEd Imports
from parmed.charmm import CharmmCrdFile, CharmmParameterSet, CharmmRstFile
from parmed.openmm.reporters import StateDataReporter
from parmed import unit as u
from parmed.charmm.psf import CharmmPsfFile
import parmed as pmd
from parmed.namd import NamdBinVel, NamdBinCoor


def from_vel(path):
    vel = NamdBinVel.read(path)
    velocities = u.Quantity(vel.velocities[0], unit=u.angstroms/u.picosecond)
    return velocities


def from_coor(path):
    namdbincoor = NamdBinCoor.read(path)
    positions = u.Quantity(namdbincoor.coordinates[0], unit=u.angstroms)
    return positions

def from_rst(path):
    namdrstfile = CharmmRstFile()
    positions = namdrstfile.positions # Angstrom
    velocities = namdrstfile.velocities # Angstrom/picosecond
    coordinates = namdrstfile.coordinates
    resnames = namdrstfile.resname
    box = namdrstfile.box


def from_xsc(path):
    """ Returns u.Quantity with box vectors from XSC file """

    def parse(path):
        """
        Open and parses an XSC file into its fields
        Parameters
        ----------
        path : str
            Path to XSC file
        Returns
        -------
        namedxsc : namedtuple
            A namedtuple with XSC fields as names
        """
        with open(path) as f:
            lines = f.readlines()
        NamedXsc = namedtuple('NamedXsc', lines[1].split()[1:])
        return NamedXsc(*[float(x) for x in lines[2].split()])

    xsc = parse(path)
    return u.Quantity([[xsc.a_x, xsc.a_y, xsc.a_z],
                       [xsc.b_x, xsc.b_y, xsc.b_z],
                       [xsc.c_x, xsc.c_y, xsc.c_z]], unit=u.angstroms)


pdb_file = 'pilb_ionized.pdb'
psf_file = 'pilb_ionized.psf'
binVel = 'pilb_out01.vel'
binCoord = 'pilb_out01.restart.coor'
xsc_file = 'pilb_out01.restart.xsc'


velocities = from_vel(binVel)
coords = from_coor(binCoord) 
box_shape = from_xsc(xsc_file)

print(box_shape)

# Load the CHARMM files
print('Loading CHARMM files...')
params = CharmmParameterSet('par_all36_prot.prm', 'toppar_water_ions.str')
params.condense(do_dihedrals=True)

# pdb_file = pmd.load_file(pdb_file)
# print(pdb_file.box_vectors)

namd_out = CharmmPsfFile(psf_file)

namd_out.load_parameters(params, copy_parameters=False)

"""
# DEFINE MANUALLY BOX SHAPE FROM LIMITS OF STRUCRUTE

min_crds = [coords[0][0], coords[0][1], coords[0][2]]
max_crds = [coords[0][0], coords[0][1], coords[0][2]]

for coord in coords:
    min_crds[0] = min(min_crds[0], coord[0])
    min_crds[1] = min(min_crds[1], coord[1])
    min_crds[2] = min(min_crds[2], coord[2])
    max_crds[0] = max(max_crds[0], coord[0])
    max_crds[1] = max(max_crds[1], coord[1])
    max_crds[2] = max(max_crds[2], coord[2])

"""
print('Loading Box Vectors to System')
namd_out.box_vectors= box_shape

# Create the OpenMM system
print('Creating OpenMM System')
system = namd_out.createSystem(params, nonbondedMethod=app.PME,
                                nonbondedCutoff=12.0*u.angstroms,
                                constraints=app.HBonds,
                                switchDistance=10.0*u.angstroms,
)

# Create the integrator to do Langevin dynamics
integrator = mm.LangevinIntegrator(
                        310*u.kelvin,       # Temperature of heat bath
                        5.0/u.picoseconds,  # Friction coefficient
                        2.0*u.femtoseconds, # Time step
)

# Define the platform to use; CUDA, OpenCL, CPU, or Reference. Or do not specify
# the platform to use the default (fastest) platform
platform = mm.Platform.getPlatformByName('CUDA')
prop = dict(CudaPrecision='mixed') # Use mixed single/double precision

# Create the Simulation object
sim = app.Simulation(ala2_solv.topology, system, integrator, platform, prop)


# Set the particle positions
print("Loading atoms positions")
sim.context.setPositions(coords)
sim.context.setVelocities(velocities)

# Minimize the energy
print('Minimizing energy')
sim.minimizeEnergy(maxIterations=500)

# Set up the reporters to report energies and coordinates every 100 steps
sim.reporters.append(
        StateDataReporter(sys.stdout, 100, step=True, potentialEnergy=True,
                          kineticEnergy=True, temperature=True,
                          volume=True, density=True)
)
sim.reporters.append(app.DCDReporter('namd_out.dcd', 100))

# Run dynamics
print('Running dynamics')
sim.step(10000)