from os.path import join as pjoin
from pkg_resources import resource_filename


def get_resource(fname=None):
    """ Function for getting `abagen` test resources """
    path = resource_filename('abagen', 'data')
    return pjoin(path, fname) if fname is not None else path
