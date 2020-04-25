"""Extentions for PDB2PQR Suite

This module provides various utilities for the PDB2PQR suite to be imported
into other Python scripts.

The "template" file is a template for interfacing with PDB2PQR.
"""

import logging
import pkgutil
from optparse import OptionGroup, OptionConflictError, Option


_LOGGER = logging.getLogger(__name__)
_extList = [name for _, name, _ in pkgutil.iter_modules(__path__)]
extDict = {}


for extName in _extList:
    extDict[extName] = __import__(extName,globals(),locals(),[],1)


def setupExtensionsOptions(parser):
    """
    Takes an instance of an OptionParser
    and adds the options for all extensions

    If an extension adds it's own options, those
    options are put in their own group.
    """

    if len(extDict) == 0:
        return None

    firstGroup = OptionGroup(parser,"Extension options")
    groups = [firstGroup]

    for extName, extModule in list(extDict.items()):
        helpArg = {}
        if hasattr(extModule, 'usage'):
            helpArg['help'] = extModule.usage()

        extOption = Option('--' + extName, action='append_const',
                             const = extName, dest = 'active_extensions', **helpArg )

        try:
            if hasattr(extModule, 'addExtensionOptions'):
                group = OptionGroup(parser, extName.capitalize() + " extension options")
                group.add_option(extOption)

                extModule.addExtensionOptions(group)

                if len(group.option_list) == 1:
                    opt = group.option_list[0]
                    group.remove_option(opt.get_opt_string())
                    firstGroup.add_option(opt)
                else:
                    groups.append(group)

            else:
                firstGroup.add_option(extOption)

        except OptionConflictError as value:
            _LOGGER.error('Error adding command line options for extension ' + extName + ' ' + '(' + str(value) + ')')

    for group in groups:
        parser.add_option_group(group)

    return groups

