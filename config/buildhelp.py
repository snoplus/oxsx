from __future__ import print_function
from dependency   import Dependency, VALID_FIELDS
try:
    from configparser import ConfigParser  # Python 3 
except ImportError:
    from ConfigParser import SafeConfigParser as ConfigParser  # Python 2
from SCons.Script import Split, Glob, Copy, File, Dir, Exit
from platform import system
from subprocess import check_call
import os

def read_dependencies(filename):
    cparse = ConfigParser()
    cparse.read(filename)
    
    dependencies = {}
    for dep_name in cparse.sections():
        try:
            libs   = cparse.get(dep_name, "libs")
            cheads = cparse.get(dep_name, "check_headers")
        except KeyError:
            print("Incomplete dependency spec for {0}, (needs libs & check_headers)".format(dep_name))
        dependencies[dep_name] = Dependency(dep_name, libs, cheads)
    return dependencies

def check_dependency(conf, dependency):
    print("\nChecking {0} dependencies..".format(dependency.name))
    for header in Split(dependency.check_headers):
        if not conf.CheckCXXHeader(header):
            print('!! Cannot locate header {0} ...'.format(header))
            Exit(0)
                
    # check libs
    for lib in Split(dependency.libs):
        if not conf.CheckLib(lib):
            print('!! Cannot locate library {0}'.format(lib))
	    print("Checking library", lib, "with env:")
            print("  LIBPATH:", conf.env['LIBPATH'])
            print("  LIBS:", conf.env['LIBS'])
	    print("  LINKFLAGS:", conf.env['LINKFLAGS'])
	    Exit(0)
        

def sanity_checks(conf):
    if not conf.CheckCXX():
        print('!! C++ compiler is not available...')
        Exit(0)


def parse_user_config(filename, dependencies):
    '''
    Read user_config and update the relevant fields
    '''
    cparse = ConfigParser()
    cparse.read(filename)
    for dep_name in cparse.sections():
        if not dep_name in dependencies.keys():
            print('!! user_config specifies unknown dependency "{0}"'.format(dep_name))
            Exit(0)
        
        # loop over options for that dependencies
        for opt in cparse.options(dep_name):
            if opt in VALID_FIELDS:
                dependencies[dep_name].__setattr__(opt, cparse.get(dep_name, opt))
            else:
                print("Unknown build option: {0} for dependency {1}".format(opt, dep_name))
    
def update_and_check_env(conf, dependencies):
    '''
    Update the build environment with the dependencies
    '''
    for dep in dependencies.values():
        if dep.header_path:
            conf.env.Append(CPPPATH = [dep.header_path])
        if dep.flags:
            conf.env.Append(CFLAGS = [dep.flags])
        if dep.libs:
            conf.env.Append(LIBS = Split(dep.libs))
        if dep.lib_path:
            conf.env.Append(LIBPATH = [dep.lib_path])
            conf.env.Append(RPATH   = [dep.lib_path])
            if conf.env["SYSTEM"] == "Darwin":
                conf.env.Append(LINKFLAGS = "-rpath {0}".format(dep.lib_path))
        check_dependency(conf, dep)
    
def create_gsl_cpy_commands(conf, dependencies, copy_folder):
    '''
    Create os dependent commands. On darwin: copy all gsl libs, fix
    the install names for dylibs using install_name_tool, and 
    replace lib path with the patched version. On linux: do nothing
    '''
    if conf.env["SYSTEM"] == "Darwin" and dependencies["gsl"].lib_path:
        lib_path = dependencies["gsl"].lib_path
        commands = []
        
        for lib in Glob(os.path.join(lib_path, "*")):
            new_path = os.path.join(copy_folder, 
                                    os.path.basename(lib.rstr()))
            action = [Copy("$TARGET", "$SOURCE")]
            if ("0.dylib" in lib.rstr()):
                action += [fix_dylib_for_darwin]

            kw = {
                'target' : '{0}'.format(new_path),
                'source' : '{0}'.format(lib),
                'action' : action
                }

            commands.append(kw)
        
        dependencies["gsl"].lib_path = Dir(copy_folder).abspath
        return commands

    else:
        return []

def fix_dylib_for_darwin(target, source, env):
    '''
    Fix the the install_names for darwin for all dylibs in target
    (here they are set to the abspath of the library)
    '''
    for t in target:
        abspath = File(t).abspath
        check_call("install_name_tool -id {0} {0}".format(abspath),
                   shell = True)
