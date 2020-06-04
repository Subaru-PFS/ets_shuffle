import pkg_resources
__version__ = pkg_resources.get_distribution('hetdex-shuffle').version

__svn_revision__ = '$Revision: 211 $'

__full_version__ = '''Shuffle version: {}.
Svn revision: {}'''.format(__version__, __svn_revision__)
