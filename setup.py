import distutils.cmd
import distutils.log
import setuptools
import subprocess
import os


class InstallCommand(distutils.cmd.Command):
    """A custom command to install requested package to run Survey Strategy Support Pipeline"""

    description = 'Script to install requested package to run Survey Strategy Support Pipeline'
    user_options = [
        # The format is (long option, short option, description).
        ('package=', None, 'package to install'),
        ('branch=', None, 'git branch to consider')
    ]

    def initialize_options(self):
        """Set default values for options."""
        # Each user option must be listed here with their default value.
        self.package = ''
        self.branch = 'dev'

    def finalize_options(self):
        """Post-process options."""
        if self.package:
            if self.package not in ['metric', 'simulation']:
                print('{} impossible to install'.format(self.package))

    def run(self):
        """Run command."""
        """
        command = ['source']
        if self.release:
            command.append('%s' % self.release)
            # command.append(os.getcwd())
            # command.append('setup lsst_sims')
            self.announce(
                'Running command: %s' % str(command),
                level=distutils.log.INFO)
            subprocess.check_call(command)
        """
        self.install()

    def clone(self, pack):
        if not os.path.isdir(pack):
            cmd = 'git clone -b {} https://github.com/lsstdesc/{}.git'.format(
                self.branch, pack)
            os.system(cmd)

    def install(self):
        pwd = os.getcwd()
        """
        lib = 'h5py'
        thedir = '{}/lib/python3.6/site-packages/'.format(pwd, lib)
        if not os.path.isdir(thedir):
            print('The directory {} does not exist -> installing '.format(thedir))
            cmd = 'pip install --prefix={} {}==2.7.1'.format(pwd, lib)
            print(cmd)
            os.system(cmd)
        """
        git_packs = ['sn_utils', 'sn_catalog_simulations']
        npack = 1
        if self.package == 'simulation':
            npack = len(git_packs)

        for pack in git_packs[:npack]:
            self.clone(pack)


setuptools.setup(

    cmdclass={
        'install': InstallCommand,
    },
    # Usual setup() args.
    # ...
)
