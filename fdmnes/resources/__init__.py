import os 

def resource_filename(resource):
    """
    Args:
        resource(str): resource path relative to resource directory

    Returns:
        str: absolute resource path in the file system
    """
    return os.path.join(os.path.abspath(os.path.dirname(__file__)),
                            *resource.split('/'))

