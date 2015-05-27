"""
ui

Helpers for user interfacing with sats.
"""

import sys

from sats.ui.args import get
from options import format_defaults


description = "sats structure generation tool. Generate some structures"

if get('ui', 'help'):
    print("usage: {0} [section.option=value ...]".format(sys.argv[0]))
    print("")
    print(description)
    print("")
    print("\n".join(format_defaults()))
#    raise SystemExit
