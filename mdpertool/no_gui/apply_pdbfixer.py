from pdbfixer import PDBFixer
from openmm.app.pdbfile import PDBFile
from openmm import app
import os


def _build_fixed_pdb_filename(source_path_or_id, pH):
    source_name = os.path.basename(str(source_path_or_id))
    stem, _ = os.path.splitext(source_name)
    return f"{stem}_fixed_pH_{pH}.pdb"


def _log(log_obj, message):
    if log_obj is not None:
        log_obj.info(message)


def _prepare_fixed_structure(source_pdb, pH, log_obj=None):
    _log(log_obj, "Creating PDBFixer...")
    fixer = PDBFixer(source_pdb)
    _log(log_obj, "Finding missing residues...")
    fixer.findMissingResidues()
    _log(log_obj, "Finding nonstandard residues...")
    fixer.findNonstandardResidues()
    _log(log_obj, "Replacing nonstandard residues...")
    fixer.replaceNonstandardResidues()
    _log(log_obj, "Removing heterogens...")
    fixer.removeHeterogens(keepWater=False)
    _log(log_obj, "Finding missing atoms...")
    fixer.findMissingAtoms()
    _log(log_obj, "Adding missing atoms...")
    fixer.addMissingAtoms()
    _log(log_obj, "Adding missing hydrogens...")
    fixer.addMissingHydrogens(pH)
    _log(log_obj, "Writing PDB file...")
    return fixer


def _write_fixed_structure(fixer, target_path):
    with open(target_path, "w") as output_handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, output_handle, keepIds=True)


def apply_mutation(pdb_path, output=None, mut_region=None, chain_id=None):
    """
    Make a mutant protein easy.

        Parameters
        ----------

        pdb_path: Give your pdb whole path to this parameter

        mut_region : list of strings
            Each string must include the resName (original), index,
            and resName (target).  For example, ALA-133-GLY will mutate
            alanine 133 to glycine.

        chain_id : str
            Chain ID to apply mutation.

        Example
        ----------

        mutate('C:/Users/HIbrahim/Desktop/MolDynAnalyze/test/last.pdb', mut_region=['ASP-306-ARG'], chain_id='A')

    """

    try:
        pdb_name = os.path.basename(pdb_path).split('.')[0]
        output_directory = output or os.path.dirname(pdb_path) or os.getcwd()

        mut_file_name = pdb_name + '_chain' + chain_id + '_' + str(mut_region[0]) + '.pdb'
        mut_file_path = os.path.join(output_directory, mut_file_name)

        fixer = PDBFixer(pdb_path)
        fixer.applyMutations(mut_region, chain_id)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        with open(mut_file_path, 'w') as w_file:
            app.PDBFile.writeFile(fixer.topology, fixer.positions, w_file, keepIds=True)

        return mut_file_path

    except ValueError as error:
        print("CANNOT APPLIED MUTATION !", error)
        # return pdb_path
        print('Please Check Your Input Parameters !!')


def fix_pdb(pdb_id, pH=7.4, output=None, log_obj=None, mutate=None, mut_region=None, chain_id=None):
    if output is not None:
        path = output
    else:
        path = os.getcwd()

    if len(pdb_id) == 4:
        return None

    source_pdb = pdb_id
    if mutate:
        _log(log_obj, "Mutating {} in Chain {} ...".format(mut_region, chain_id))
        source_pdb = apply_mutation(
            pdb_path=pdb_id,
            output=output,
            mut_region=mut_region,
            chain_id=chain_id,
        )

    fixer = _prepare_fixed_structure(source_pdb=source_pdb, pH=pH, log_obj=log_obj)
    fixed_pdb_filename = _build_fixed_pdb_filename(source_pdb, pH)
    fixed_pdb_path = os.path.join(path, fixed_pdb_filename)
    _write_fixed_structure(fixer=fixer, target_path=fixed_pdb_path)
    return fixed_pdb_path
