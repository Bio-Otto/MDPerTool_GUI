import os.path as op

global readme, changelog, __citation__, __long_description__, __url__, __author_email__


__version_info__ = ('0', '0', '1', '')

__author__ = 'Halil I. Ozdemir'

__credits__ = 'Marmara University - Ozbek Lab'

__description__ = 'A Software Tool for Investigation of Allosteric Communication within Protein Structures via Energy ' \
                  'Dissipation in Molecular Dynamics Simulations '

__url__ = 'https://github.com/bio-otto/MDPerTool_GUI'

__author_email__ = 'halil.ibrahim.oozdemir@gmail.com'

__version__ = '.'.join(__version_info__[:3])
if len(__version_info__) == 4:
    __version__ += __version_info__[-1]


with open(op.join(op.dirname(op.realpath(__file__)), 'README.md')) as readme_file:
    readme = readme_file.read()

# with open(op.join(op.dirname(op.realpath(__file__)), 'CHANGELOG.md')) as changelog_file:
#     changelog = changelog_file.read()

#with open(op.join(op.dirname(op.realpath(__file__)), 'CITATION.md')) as citation_file:
#    __citation__ = citation_file.read()

__long_description__ = """The Program applying Energy Dissipation Concept using powerfull OpenMM Molecular Dynamic 
Toolkit, which also supports the Cuda platform. Each residual energy calculation required for the concept can be 
calculated using OpenMMs flexible and useful infrastructure.In addition, you can use the package only for energy 
decomposition. For this, it will be sufficient to specify a XTC or a DCD file in the script. """
# desc = readme + '\n\n' + changelog + '\n\n' + __citation__
# try:
#     import pypandoc
#
#     __long_description__ = pypandoc.convert_text(desc, 'rst', format='md')
#     with open(op.join(op.dirname(op.realpath(__file__)), 'README.rst'), 'w') as rst_readme:
#         rst_readme.write(__long_description__)
# except (ImportError, OSError, IOError):
#     __long_description__ = desc