"""
args.py

Simple commandline parser utilising argparse.

"""

import sys
from collections import OrderedDict

from sats.ui.options import format_defaults



def parse_args():
    """
    Parse arguments from the commandline.

    Returns
    -------
    args : dict of str
        All the options from the command line as a dict[section][option].
        values are kept as strings.
    """

    args = OrderedDict()

    description = "sats structure generation tool. Generate some structures"

    for arg in sys.argv[1:]:
        if arg in ['-h', '--help', 'help', 'ui.help']:
            print("usage: {0} [section.option=value ...]".format(sys.argv[0]))
            print("")
            print(description)
            print("")
            print("\n".join(format_defaults()))
            raise SystemExit
        elif '=' in arg:
            option, value = arg.split('=')
            section, option = option.split('.')
        else:
            value = "True"
            section, option = arg.split('.')

        section = section.strip()
        option = option.strip()
        value = value.strip()

        if section not in args:
            args[section] = OrderedDict()

        args[section][option] = value

_parsed_args = parse_args()


def get(section, option):
    """
    Return the value of the section.option if is has been specified on the
    command line. Otherwise return None.

    Parameters
    ----------
    section : str
        The section where the option is located.
    option : str
        The name of the option to retrieve.

    Returns
    -------
    value : str or None
        The value of the section.option option. The value will be a string,

    """
    if section in _parsed_args and option in _parsed_args[section][option]:
        return _parsed_args[section][option]
    else:
        return None