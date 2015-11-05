###############################################################################
#
# build_utils.py
#
# Builders and utilities to support scons building.
#
###############################################################################

from SCons.Script import *
import SCons.Errors
import re
from test_tools import (ComparePQRAction,
                        ComparePROPKAAction,
                        CompareTitCurvesAction,
                        CompareStringFunc,
                        CompareDirectoryFunc)

# Taken from the SubstInFileBuilder on the SCons Wiki.
# See http://www.scons.org/wiki/AutoconfRecipes

# Also changed to use string.replace as using re.sub was a dumb idea if
# we are going to be replicating the behavior of autotools. Now using
# re.sub is optional -KEM

def TOOL_SUBST(env):
    """Adds ACGenerateFile builder, which substitutes the keys->values of
    SUBST_DICT from the source to the target.

    The values of SUBST_DICT first have any construction variables expanded
    (its keys are not expanded).

    If a value of SUBST_DICT is a python callable function, it is called and
    the result is expanded as the value.

    If there's more than one source and more than one target, each target gets
    substituted from the corresponding source.
    """
    env.Append(TOOLS = 'SUBST')

    def do_subst_in_file(targetfile, sourcefile, dict, useRegex):
        """Replace all instances of the keys of dict with their values.
        For example, if dict is {'VERSION': '1.2345', 'BASE': 'MyProg'},
        then all instances of VERSION in the file will be replaced with
        1.2345 etc.
        """
        try:
            f = open(sourcefile, 'rb')
            contents = f.read()
            f.close()
        except:
            raise SCons.Errors.UserError, "Can't read source file %s"%sourcefile
        if useRegex:
            for (k,v) in dict.items():
                #In 2.6 re.sub has no "flags" argument.
                mo = re.compile(k, flags=re.DOTALL)
                contents = mo.sub(v, contents)
        else:
            for (k,v) in dict.items():
                contents = contents.replace(k, v)
        try:
            f = open(targetfile, 'wb')
            f.write(contents)
            f.close()
        except:
            raise SCons.Errors.UserError, "Can't write target file %s"%targetfile
        return 0 # success

    def subst_in_file(target, source, env):
        useRegex = env['SUBST_USE_REGEX'] if 'SUBST_USE_REGEX' in env else False
        if 'SUBST_DICT' not in env:
            raise SCons.Errors.UserError, "SubstInFile requires SUBST_DICT to be set."
        d = dict(env['SUBST_DICT']) # copy it
        for (k,v) in d.items():
            if callable(v):
                d[k] = env.subst(v())
            elif SCons.Util.is_String(v):
                d[k]=env.subst(v)
            else:
                raise SCons.Errors.UserError, "SubstInFile: key %s: %s must be a string or callable"%(k, repr(v))
        for (t,s) in zip(target, source):
            return do_subst_in_file(str(t), str(s), d, useRegex)

    def subst_in_file_string(target, source, env):
        """This is what gets printed on the console."""
        return '\n'.join(['Substituting vars from %s into %s'%(s, t)
                          for t,s in zip(target, source)])

    def subst_emitter(target, source, env):
        """Add dependency from substituted SUBST_DICT to target.
        Returns original target, source tuple unchanged.
        """
        d = env['SUBST_DICT'].copy() # copy it
        for (k,v) in d.items():
            if callable(v):
                d[k] = env.subst(v())
            elif SCons.Util.is_String(v):
                d[k]=env.subst(v)
        Depends(target, SCons.Node.Python.Value(d))
        return target, source

    subst_action=SCons.Action.Action(subst_in_file, subst_in_file_string)
    env['BUILDERS']['SubstInFile'] = Builder(action=subst_action, emitter=subst_emitter)

def CopySubAction(targetfile, sourcefile, dict, useRegex=False):
    """Replace all instances of the keys of dict with their values.
    For example, if dict is {'VERSION': '1.2345', 'BASE': 'MyProg'},
    then all instances of VERSION in the file will be replaced with
    1.2345 etc.
    """
    try:
        f = open(sourcefile, 'rb')
        contents = f.read()
        f.close()
    except:
        raise SCons.Errors.UserError, "Can't read source file %s"%sourcefile
    if useRegex:
        for (k,v) in dict.items():
            #In 2.6 re.sub has no "flags" argument.
            mo = re.compile(k, flags=re.DOTALL)
            contents = mo.sub(v, contents)
    else:
        for (k,v) in dict.items():
            contents = contents.replace(k, v)
    try:
        f = open(targetfile, 'wb')
        f.write(contents)
        f.close()
    except:
        raise SCons.Errors.UserError, "Can't write target file %s"%targetfile
    return 0 # success

def CopySubActionStringFunc(targetfile, sourcefile, dict, useRegex):
    return 'CopySubAction("%s", "%s")' % (targetfile, sourcefile)

def CompilePythonAction(targetfile, sourcefile):
    """Compile python into byte code.
    """
    try:
        import py_compile
    except ImportError:
        raise SCons.Errors.InternalError, "Couldn't import py_compile module"

    try:
        py_compile.compile(sourcefile, targetfile, doraise=True)
    except py_compile.PyCompileError:
        raise SCons.Errors.InternalError, "Couldn't compile {0}".format(sourcefile)



def CompilePythonActionStringFunc(targetfile, sourcefile):
    return 'Compiling python to bytecode ("%s", "%s")' % (targetfile, sourcefile)

CopySub = SCons.Action.ActionFactory( CopySubAction, CopySubActionStringFunc )

ComparePQR = SCons.Action.ActionFactory( ComparePQRAction, CompareStringFunc )
ComparePROPKA = SCons.Action.ActionFactory( ComparePROPKAAction, CompareStringFunc )
CompareTitCurves = SCons.Action.ActionFactory( CompareTitCurvesAction, CompareDirectoryFunc )

CompilePython = SCons.Action.ActionFactory( CompilePythonAction, CompilePythonActionStringFunc )




