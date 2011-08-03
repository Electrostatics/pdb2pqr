import pkgutil

from optparse import OptionGroup, OptionConflictError

_extList = [name for _, name, _ in pkgutil.iter_modules(__path__)]

extDict = {}

for extName in _extList:
    extDict[extName] = __import__(extName,globals(),locals(),[], -1)
    
def setupExtensionsOptions(parser):
    '''Takes an instance of an OptionParser 
    and adds the options for all extensions'''
    
    if len(extDict) == 0:
        return None
    
    group = OptionGroup(parser,"Extension options", "Options for output extensions.")
    for extName, extModule in extDict.items():
        helpArg = {}
        if hasattr(extModule, 'usage'):
            helpArg['help'] = extModule.usage()
        
        group.add_option('--' + extName, action='append_const', const = extName, dest = 'active_extentions', **helpArg )
        
        if hasattr(extModule, 'addExtensionOptions'):
            try:
                extModule.addExtensionOptions(group)
            except OptionConflictError, value:
                print 'Error adding command line option for extension ' + extName + ' ' + '(' + str(value) + ')'
    
    parser.add_option_group(group)
    
    return group