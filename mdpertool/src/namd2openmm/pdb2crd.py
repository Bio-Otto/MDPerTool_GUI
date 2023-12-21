import argparse
import re
import sys

_s = re.compile('\s+')

def parse_psf(psffile):
    psf = open(psffile)
    psf.readline(); psf.readline()
    ntitle = int(psf.readline().strip().split()[0])

    title = []
    for i in range(ntitle):
        title.append(psf.readline().rstrip())

    psf.readline()
    natoms = int(psf.readline().strip().split()[0])

    atoms = []
    for i in range(natoms):
        line = psf.readline()
        num, segid, resid, resname, name = _s.split(line)[1:6]
        atoms.append({'num': int(num), 'segid': segid, 'resid': int(resid), 'resname': resname, 'name': name})
    return title, atoms

def parse_pdb(pdbfile):
    pdb = open(pdbfile)
    atoms = []
    for line in pdb.readlines():
        if not line.startswith('ATOM'): continue
        x,y,z = map(float, [line[k:k+8] for k in (30, 38, 46)])
        b = float(line[60:66])
        atoms.append((x,y,z,b))
    return atoms

def write_crd(title, psf, pdb, crdfile):
    crd.write("\n".join(title))
    crd.write("\n*\n")
    crd.write("%10d  EXT\n" % len(psf))

    prev_resid = -1
    resid = 0
    for i,atom in enumerate(psf):
        if prev_resid != atom['resid']: resid += 1
        crd.write("  %8d  %8d  %-8s  %-8s  %18.10f  %18.10f  %18.10f  %-8s  %-8d  %18.10f\n" % (atom['num'], resid, atom['resname'], atom['name'], pdb[i][0], pdb[i][1], pdb[i][2], atom['segid'], atom['resid'], pdb[i][3]))
        prev_resid = atom['resid']
"""
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('psf', metavar='psffile', help='PSF file')
    parser.add_argument('pdb', metavar='pdbfile', help='PDB file')
    parser.add_argument('-o' , dest='output', help='output CRD file', default=None)
    args = parser.parse_args()

    title, psf = parse_psf(args.psf)
    pdb = parse_pdb(args.pdb)

    if not args.output: crd = sys.stdout
    else: crd = open(args.output, 'w')
    write_crd(title, psf, pdb, crd)


"""
import parmed as pmd
parm = pmd.load_file('pilb_ionized.psf')
parm.coordinates = pmd.load_file('pilb_ionized.pdb').coordinates

parm.save('charmm.crd', format='charmmcrd')