
import os
from shutil import copy

def copy_or_write(input_item, destination_path, **kwargs):
    """
    Copy or write an input object (e.g., photosphere, transitions) to the given
    destination path.

    :param input_item:
        The relevant input object. This could be a photosphere, a set of transitions,
        or a path to one of those items.
    
    :param destination_path:
        The path where the input object should be written.
    """

    os.makedirs(os.path.dirname(destination_path), exist_ok=True)

    if isinstance(input_item, (str, bytes)):    
        copy(input_item, destination_path)

    else:
        input_item.write(destination_path, **kwargs)
        

