__version_info__ = ('0', '0', '1', '')
__version__ = '.'.join(__version_info__[:3])
if len(__version_info__) == 4:
    __version__ += __version_info__[-1]
