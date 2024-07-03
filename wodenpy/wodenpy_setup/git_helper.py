from subprocess import check_output
import importlib_resources
import os

def get_commandline_output(command_list):
    """
    Takes a command line entry separated into list entries, and returns the
    output from the command line as a string

    Parameters
    ----------
    command_list : list of strings
        list of strings that when combined form a coherent command to input into
        the command line

    Returns
    -------
    output : string
        the output result of running the command

    """
    output = check_output(command_list,universal_newlines=True).strip()
    return output

def make_gitdict():
    """
    Makes a dictionary containing key git information about the repo by running
    specific commands on the command line

    Returns
    -------
    git_dict : dictionary
        A dictionary containing git information with keywords: describe, date,
        branch

    """

    ##Try and get a git version. If this is a release version, it might not
    ##be in a git repo so try, otherwise return False
    try:
        git_dict = {
            'describe': get_commandline_output(["git", "describe", "--always"]),
            'date': get_commandline_output(["git", "log", "-1", "--format=%cd"]),
            'branch': get_commandline_output(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
        }
        
    except:
        git_dict = False

    return git_dict

def retrieve_gitdict():
    """
    Attempt to recover the git information that was created during pip install. If found, return the git dictionary. If not, return False

    Returns
    -------
    git_dict : dictionary
        A dictionary containing git information with keywords: describe, date,
        branch

    """

    file_path = importlib_resources.files('wodenpy').joinpath('wodenpy_gitinfo.txt')

    ##If things have been pip installed correctly
    if os.path.isfile(file_path):
        git_dict = {}
        
        with open(file_path, "r") as infile:
            for line in infile.read().split('\n'):
                
                if line == '':
                    pass
                else:
                    key,entry = line.split(',')
                    git_dict[key] = entry
                
    else:
        git_dict = False
        
    return git_dict
