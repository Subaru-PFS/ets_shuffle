#!/usr/bin/env python
"""
Created on Mon April 11 03:52:02 2016

This script is used to call do_shuffle.py for a single target.
It allows you to specify the target RA, Dec, the desired
IFU slot for it to land on, a small offset from the center of that
IFU if needed, and the track of the night.

Using the default shuffle.cfg, a visualization tool pops up for
the new shuffled coordinates and allows an inspection of the field.

Also produced are ACAM images along with the (RA, Dec) of the GP(s),
(RA, Dec) of the IHMP center, and the (x, y) pixel of the star(s) in the
ACAM image for the Telescope operator.

@author: gregz
"""
import builtins

import argparse as ap
import inspect
import logging
import os
# from pprint import pformat
import sys
import textwrap as tw
import warnings
import datetime
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time

import numpy
import configparser

from . import shuffle
from . import findStars
from . import __full_version__


class IFUException(ValueError):
    '''Exception raised when the IFU does not exist'''
    pass


def load_config(configfile):
    """Read configuration file

    Parameters
    ----------
    configfile : string
        name of the configuration file

    Returns
    -------
    config : :class:`~ConfigParser.ConfigParser` instance
    """
    config = configparser.ConfigParser(defaults={'xpa_method': 'local'})
    if not config.read(configfile):
        msg = ("File '{}' does not exist. Using the default one."
               " You can get the full set of configuration files with the"
               " command: ``shuffle_config``".format(configfile))
        warnings.warn(msg)
        config_def = os.path.join(os.path.dirname(__file__), 'configs',
                                  'shuffle.cfg')
        if not config.read(config_def):
            msg = ("ERROR: Failed to load the default configuration file.")
            raise ap.ArgumentTypeError(msg.format(configfile))

    return config


class RaToFloat(ap.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        '''Convert ra to floats.

        Parameters
        ----------
        parser : current parser
        namespace : name space
        value : string
            value from the command line
        option_string : string
        '''
        ra = values.split(':')
        try:
            if len(ra) == 1:
                ra_deg = float(ra[0])
            else:
                rah, ram, ras = ra
                ra_deg = float(rah) + float(ram)/60. + float(ras)/3600.
        except ValueError:
            raise parser.error('Ra can be provided as number or as string of'
                               ' colon separated numbers, like "hh:mm:ss".'
                               ' "{}" is not allowed'.format(values))
        setattr(namespace, self.dest, ra_deg)


class DecToFloat(ap.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        '''Convert the dec to float

        Parameters
        ----------
        parser : current parser
        namespace : name space
        value : string
            value from the command line
        option_string : string
        '''
        dec = values.strip().split(':')
        try:
            if len(dec) == 1:
                dec_deg = float(dec[0])
            else:
                decd, decm, decs = dec
                if decd.startswith('-'):
                    dec_deg = float(decd) - float(decm)/60. - float(decs)/3600.
                else:
                    dec_deg = float(decd) + float(decm)/60. + float(decs)/3600.
        except ValueError:
            raise parser.error('Dec can be provided as number or as string of'
                               ' colon separated numbers, like "dd:mm:ss".'
                               ' "{}" is not allowed'.format(values))
        setattr(namespace, self.dest, dec_deg)


class PositionalParser(ap.ArgumentParser):
    def error(self, message):
        msg = tw.dedent('''
                        Are you by any chance trying to provide a negative
                        ``dec`` in the form ``dd:mm:ss`` to ``%s``? If the
                        answer is yes, use the float representation of dec or
                        add ``--`` before the positional arguments. E.g.:

                            do_shuffle [options] 23:40:05.43 -1.32195
                                       radius {0,1,2} ifuslot [ra_offset]
                                       [dec_offset] [target]

                        or

                            do_shuffle [options] -- 23:40:05.43 -1:19:19.02
                                       radius {0,1,2} ifuslot [ra_offset]
                                       [dec_offset] [target]
                        ''' % (self.prog))
        print_(msg)

        super(PositionalParser, self).error(message)


def parse_args(argv=None):
    """Parse the command line arguments

    Parameters
    ----------
    argv : list of string
        arguments to parse; if ``None``, ``sys.argv`` is used

    Returns
    -------
    Namespace
        parsed arguments
    """
    description = "Shuffle the HETDEX shots"

    parser = PositionalParser(description=description,
                              formatter_class=ap.ArgumentDefaultsHelpFormatter,
                              )

    parser.add_argument("--verbose", '-v', action="count", default=0,
                        help="""Increase verbosity, can be called multiple
                        times""")
    parser.add_argument('--version', '-V', action='version',
                        version=__full_version__)

    parser.add_argument("ra", action=RaToFloat, help='''ra of the FP center in
                        hours.''')

    parser.add_argument("dec", action=DecToFloat,
                        help='dec of the FP center in degrees')

    parser.add_argument("pa", type=float,
                        help='PA of the telescope in degrees')

    parser.add_argument("radius", type=float,
                        help='''Search radius (in degrees) for proper guide
                        stars.''')

    parser.add_argument("--epoch", type=float,
                        help='''Epoch of target coordinates''',
                        default=2000.0)

    parser.add_argument("--pmRA", type=float,
                        help='''Proper motion in right ascension (mas/yr)''',
                        default=None)

    parser.add_argument("--pmDEC", type=float,
                        help='''Proper motion in delination (mas/yr)''',
                        default=None)

    parser.add_argument("--gaia2radius", type=float,
                        help='''Search radius for target match to gaia2 in "''',
                        default=30.)

    parser.add_argument("-c", "--config", help="""Name of the configuration
                        file. When parsing the command line, the file is loaded
                        into a configuration object""",
                        default="./shuffle.cfg", type=load_config)

    # ARGUMENTS IN CONFIG FILE

    parser.add_argument("--localcat", type=str, help="""Use a local catalog
                        for star searches""", default=None)

    parser.add_argument("--interactive", type=str, help="""Use interactive
                        mode for guide star search.""", default=None)

    parser.add_argument("--hpf", type=str,
                        help="""Use hpf additional selection""", default=None)

    parser.add_argument("--use_brightness", type=str,
                        help="""Use the brightness
                        of the probe stars as a selection criteria.  If False,
                        then the distance from the target angle is
                        minimized.""", default=None)

    parser.add_argument("--fplane_file", type=str, help="""Fplane file to be
                        used.""", default=None)

    parser.add_argument("--catalog", type=str, help="""Star catalog to use,
                        for example: GAIA""", default=None)

    for st in ['gp', 'wfs']:
        parser.add_argument("--" + st + "pickcol", type=int,
                            help="""Probe filter to be used.
                            For SDSS DR9 1 = u, 2 = g, 3 = r, 4 = i , 5 = z
                            For USNO A2  2 = B, 3 = R""", default=None)

    parser.add_argument("--visualize", type=str, help="""Make visualization
                        for the shot.""", default=None)

    parser.add_argument("--visualize_ACAM", type=str, help="""Make visualization
                        for the ACAM.""", default=None)

    parser.add_argument("--visualize_probestars", type=str,
                        help="""Make visualization for the probestars.""",
                        default=None)

    parser.add_argument("--acam", type=float, help="""Offset angle for ACAM""",
                        default=None)

    parser.add_argument("--fplane", type=float,
                        help="""Offset angle for fplane""",
                        default=None)

    parser.add_argument("--probes", type=float,
                        help="""Offset angle for probes""",
                        default=None)

    for o in ['guide1', 'guide2', 'wfs1', 'wfs2', 'ifu', 'acam', 'fplane',
              'allprobes', 'allcams']:
        for p in ['magadd', 'nstars', 'minsep', 'magmin', 'magmax']:
            parser.add_argument('--' + o + '_' + p,
                                help="""{} {}""".format(o, p), default=None)

    args = parser.parse_args(args=argv)

    if args.allprobes_magadd is not None:
        for o in ['guide1', 'guide2', 'wfs1', 'wfs2']:
            if getattr(args, o+'_'+'magadd') is None:
                setattr(args, o+'_'+'magadd', args.allprobes_magadd)

    if args.allcams_magadd is not None:
        for o in ['guide1', 'guide2', 'wfs1', 'wfs2', 'acam']:
            if getattr(args, o+'_'+'magadd') is None:
                setattr(args, o+'_'+'magadd', args.allcams_magadd)

    if args.fplane is not None:
        args.probes = args.fplane
    return args


def cl_to_config(args):
    """Copy the values of some command line options into the configuration.

    Values copied: force_* into the corresponding entries in the ``General``
    section if the value is not ``None``

    Parameters
    ----------
    args : Namespace
        parsed command line arguments

    Returns
    -------
    args : Namespace
        parsed command line arguments
    """
    # move not None values into the General, visualization, offsets section
    # of the configuration
    general_args, viz_args, offset_args, mag_args = [], [], [], []

    list_type = ['General', 'visualisation', 'offsets', 'MagLimits']
    for o in ['gp1', 'gp2', 'wfs1', 'wfs2']:
        general_args.append('force_'+o)
    general_args.extend(['localcat', 'interactive', 'use_brightness',
                         'fplane_file', 'gppickcol', 'wfspickcol', 'visualize',
                         'visualize_ACAM', 'catalog', 'hpf'])
    for o in ['g1', 'g2', 'w1', 'w2']:
        for p in ['min', 'max', 'targ']:
            general_args.append('dpatrol_' + o + p)

    viz_args.append('visualize_probestars')

    offset_args.extend(['acam', 'fplane', 'probes'])

    for o in ['guide1', 'guide2', 'wfs1', 'wfs2', 'ifu', 'acam', 'fplane']:
        for p in ['magadd', 'nstars', 'minsep', 'magmin', 'magmax']:
            mag_args.append(o+'_'+p)

    big_list = [general_args, viz_args, offset_args, mag_args]
    for list_name, p in zip(big_list, list_type):
        for o in list_name:
            value = getattr(args, o)
            if value is not None:
                args.config.set(p, o, value)

    return args


def check_fplane_file(args):
    '''If the fplane file does not exist, in the config file point to the one
    shipped with shuffle.

    Parameters
    ----------
    args : Namespace
        parsed configuration options

    Returns
    -------
    args : Namespace
        parsed configuration options
    '''
    fplane_file = args.config.get("General", "fplane_file")
    if not os.path.exists(fplane_file):
        msg = ("File '{}' does not exist. Using the default one."
               " You can get the full set of configuration files with the"
               " command: ``shuffle_config``".format(fplane_file))
        warnings.warn(msg)
        fplane_def = os.path.join(os.path.dirname(__file__), 'configs',
                                  'fplane.txt')
        args.config.set('General', 'fplane_file', fplane_def)
    return args


def make_dirs(args):
    '''To run shuffle, some directories are needed.
    The name of the directories is given in the configuration file.
    This function tries to make make them, and ignore any error

    Parameters
    ----------
    args : Namespace
        parsed configuration options
    '''
    for i in ['cache', 'images', 'acam_images', 'ds9_regions', 'gp_images']:
        dir_name = args.config.get("directories", i)
        try:
            os.mkdir(dir_name)
        except OSError:
            pass


def setup_logging(args):
    '''Set up a logger for shuffle with a name ``shuffle``.

    Use a StreamHandler to write to stdout and set the level to DEBUG if
    verbose is set from the command line
    '''
    fmt = '[%(levelname)s - %(asctime)s] %(message)s'
    if args.verbose == 0:
        level = logging.WARNING
    elif args.verbose == 1:
        level = logging.INFO
    else:
        level = logging.DEBUG
        fmt = '[%(levelname)s - %(filename)s - %(asctime)s] %(message)s'
    fmt = logging.Formatter(fmt)

    handler = logging.StreamHandler()
    handler.setFormatter(fmt)
    handler.setLevel(level)

    log = logging.getLogger('shuffle')
    log.setLevel(logging.DEBUG)
    log.addHandler(handler)


def set_xpa_method(config):
    '''Push the xpa method to the environment before loading ds9. By default is
    set to 'local'. It can be customised with the [ds9][xpa_method] option in
    the configuration file.

    Parameters
    ----------
    config : :class:`configparser.ConfigParser` instance
        configuration
    '''

    # set the xpa method; defaults to local to avoid locks when loading pyds9
    os.environ['XPA_METHOD'] = config.get('ds9', 'xpa_method')


def update_coords_for_proper_motion(ra, dec, epoch, log, args):
    current_epoch = Time(datetime.datetime.now()).decimalyear
    if (args.pmRA is not None) and (args.pmDEC is not None):
        deltaRA = ((current_epoch - epoch) * args.pmRA /
                   1e3 / 3600. / numpy.cos(dec * numpy.pi / 180.))
        deltaDE = ((current_epoch - epoch) * args.pmDEC /
                   1e3 / 3600.)
        ra += deltaRA
        dec += deltaDE
        log.info('Change in RA from proper motion: %0.2f"' %
                 (deltaRA * 3600.))
        log.info('Change in Dec from proper motion: %0.2f"' %
                 (deltaDE * 3600.))
        return ra, dec
    try:
        pmcat = findStars.queryGAIA2(ra, dec, args.gaia2radius / 3600.)
        c1 = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
        deltaRA = ((epoch - 2000.) * pmcat[:, -2] /
                   1e3 / 3600. / numpy.cos(dec * numpy.pi / 180.))
        deltaDE = ((epoch - 2000.) * pmcat[:, -1] /
                   1e3 / 3600.)
        c2 = SkyCoord((pmcat[:, 2]+deltaRA)*u.deg, (pmcat[:, 3]+deltaDE)*u.deg,
                      frame='icrs')
        idx, d2d, d3d = c1.match_to_catalog_sky(c2)
        if d2d.arcsec < 2.:
            deltaRA = ((current_epoch - epoch) * pmcat[idx, -2] /
                       1e3 / 3600. / numpy.cos(dec * numpy.pi / 180.))
            deltaDE = ((current_epoch - epoch) * pmcat[idx, -1] /
                       1e3 / 3600.)

            idx, d2d, d3d = c1.match_to_catalog_sky(c2)
            ra += deltaRA
            dec += deltaDE
            log.info('Change in RA from proper motion: %0.2f"' %
                     (deltaRA * 3600.))
            log.info('Change in Dec from proper motion: %0.2f"' %
                     (deltaDE * 3600.))
            return ra, dec
        else:
            return ra, dec
    except Exception:
        log.warning("Could NOT update astrometry with proper motion")
        return ra, dec


def main():

    args = parse_args()
    args = cl_to_config(args)

    args = check_fplane_file(args)
    setup_logging(args)
    set_xpa_method(args.config)
    make_dirs(args)

    log = logging.getLogger('shuffle')

    # read the parameters from command line
    # assume input ra in hms, and convert to degrees by *15
    ra = args.ra * 15
    dec = args.dec
    ra, dec = update_coords_for_proper_motion(ra, dec, args.epoch, log, args)
    radius = args.radius
    pa = args.pa
    obstime = args.obstime
    # create the data array
    data = numpy.array([1, ra, dec, pa, obstime, 1, 1, 1, None])
    # reshape to 2-d array (1,data.size)
    data = data.reshape(1, data.size)
    shuffle.do(args, data, radius=radius, start=1, orig_loc=[ra, dec],
               targetID=args.target)


if __name__ == '__main__':
    main()
