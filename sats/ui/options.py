"""
options.py

Read and interpret input parameters. The types are
inferred from the defaults.

"""

import re
from ConfigParser import SafeConfigParser
from collections import OrderedDict, namedtuple

from sats.ui import args

Option = namedtuple('Option', ['parser_function', 'default_value', 'doc',
                               'multiple', 'valid_values'])

sats_ini = SafeConfigParser()
sats_ini.read('sats.ini')
DEFAULTS = OrderedDict()


def add_option(section, option, default_value, parser_function, doc,
               multiple=False, valid_values=None):
    """
    Add a option to the global namespace. Use as a decorator for any functions
    that take need options, these will be made available in the global
    options module.

    Parameters
    ----------
    section : str
        Name of section in which to store the argument.
    option : str
        Name of the option.
    default_value : any
        Default value to assign to the option. May be any type that can
        be interpreted by the options module.
    parser_function : function
        Function to parse individual values of the data, e.g. bool, int,
        float, str...
    doc : str
        Description of the option.
    multiple : bool
        If True, multiple values will be returned as a tuple.
    valid_values : list, optional
        If valid_values are give, these will be the only options that can
        be used with this option and will be added to any other valid_values
        from other calls to this function.

    Returns
    -------
    wrapper : function
        Decorating function that returns the unmodified function
    """

    # No sections exist in the beginning
    if section not in DEFAULTS:
        DEFAULTS[section] = OrderedDict()

    # merge any existing valid values
    if valid_values and option in DEFAULTS[section]:
        valid_values.extend(DEFAULTS[section][option].valid_values)

    DEFAULTS[section][option] = Option(parser_function=parser_function,
                                       default_value=default_value, doc=doc,
                                       multiple=multiple,
                                       valid_values=valid_values)

    # Return a function that does nothing so we can use this as a decorator
    def wrapper(func):
        """Do nothing but return the original function."""
        return func

    return wrapper


def get(section, option=None):
    """
    Return the value for the specified option. The value type will be the
    type expected by the default value.

    Parameters
    ----------
    section : str
        The section where the option is located.
    option : str, optional
        The name of the option to retrieve. If not given, the section is
        split as section.option.

    Returns
    -------
    value : Any
        The value of the section.option option. The type is the type specified
        by the default option.
    """

    if option is None:
        section, option = section.split('.')

    try:
        signature = DEFAULTS[section][option]
    except KeyError:
        raise ValueError("Unknown option %s.%s" % (section, option))

    parser_function = DEFAULTS[section][option].parser_function

    if args.get(section, option) is not None:
        value = args.get(section, option)
    elif sats_ini.has_section(option) and sats_ini.has_option(section, option):
        value = sats_ini.get(section, option)
    else:
        value = DEFAULTS[section][option].default_value

    # Parse through options
    if signature.multiple:
        try:
            # split on spaces and commas, throw away brackets
            value = [x for x in re.split('[\s,\(\)\[\]]*', value) if x]
        except TypeError:
            # Already a list or single value
            pass

        try:
            value = [parser_function(x) for x in value]
        except TypeError:
            # Make a list of a single value
            value = [parser_function(value)]
    else:
        value = parser_function(value)

    return value


def format_defaults(width=80):
    """
    Return a list of strings that has all of the possible options
    available.

    :return:
    """
    defaults = []
    for section in sorted(DEFAULTS):
        for option in sorted(DEFAULTS[section]):
            defaults.append("{0} {1} {2}".format(section, option,
                                                 DEFAULTS[section][option]))

    return defaults
