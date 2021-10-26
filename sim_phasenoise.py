import os
import sys
import numpy as np
import matplotlib.pyplot as plt
rng = np.random.default_rng(2021) #Setting the seed for reproducibility

def coords_to_string(coords):
    """Function to convert a (RA, DEC) set of coordinates to a string.

    Args:
      coords:
        List or tuple with RA and DEC.

    Returns:
      string with coordinates in 'J2000 XhXmXs +/-XdXmXs' format.
    """
    pra = coords[0]
    pdec = coords[1]
    Hpra = np.int(pra)
    Mpra = np.int((pra - Hpra) * 60.0)
    Spra = (pra - Hpra - Mpra / 60.0) * 3600.0
    pdsign = '+'
    if pdec < 0:
        pdsign = '-'
        pdec = -pdec
    Hpdec = np.int(pdec)
    Mpdec = np.int((pdec - Hpdec) * 60.0)
    Spdec = (pdec - Hpdec - Mpdec / 60.0) * 3600.0
    coord_string = 'J2000 {}h{}m{}s {}{}d{}m{}s'.format(
        str(Hpra), str(Mpra), str(Spra), pdsign, str(Hpdec), str(Mpdec),
        str(Spdec))
    return coord_string


def gaussian_comp_model(phase_center, comps, nu, name = 'Gauss_pointE.cl'):
    """Function that creates a component list to be used with simobserve.

    Args:
      phase_center:
        List or tuple with RA (in hours) and DEC (in degrees) of the phase
        center. Example:
          phase_center = (15.8, -35.0)
      comps:
        List with Gaussian components. Each item of the list should be another
        list with 6 items: flux, beam major axis, beam minor axis, beam PA,
        offset in RA, and offset in DEC. Example:
          comps = [[0.1,  0.001,  0.001,    0.0,  0.00,   0.00],
                   [1.0,  0.1,    0.01,     0.5,  1.0,    -1.0]]
      nu:
        Frequency (in Hz) of the emission of the model.
      name:
        String with name of the file where the components will be saved. Default
        is Gauss_pointE.cl
    """
    ngauss = len(comps)
    dcos = np.cos(phase_center[1] * np.pi / 180.0)

    # Make Gauss-pointE.cl file'
    os.system('rm -rf {}'.format(name))
    cl.done()
    for i in range(0, ngauss):
        FLUX = comps[i][0]
        MJS = str(comps[i][1]) + 'arcsec'
        MNS = str(comps[i][2]) + 'arcsec'
        PAS = str(comps[i][3]) + 'deg'
        ddelta = phase_center[1] + comps[i][5] / 3600.0
        dra = phase_center[0] + comps[i][4] / 3600.0 / 15.0 * dcos
        dirc = coords_to_string([dra, ddelta])
        cl.addcomponent(
            dir = dirc, flux = FLUX, fluxunit = 'Jy', freq = nu,
            shape = "Gaussian", majoraxis = MJS, minoraxis = MNS,
            positionangle = PAS)
    cl.rename(name)
    cl.close()


def phase_rms(antbl, scale_prms = 1.0, nu_scale = None):
    """Function that computes phase rms for a given baseline.

    The structure function here gives 30 deg phase rms at 10000m = 10km, which
    is the limit "go" conditions when performing a go-nogo at ALMA.

    Args:
      antbl:
        Antenna baseline in m.
      scale_prms:
        (Optional) Factor to scale the phase rms, in order to simulate better or
        worse observing conditions. Default is 1, producing a phase rms of
        30 deg at 10 km.
      nu_scale:
        (Optional) Frequency of the observation, in case one wants to simulate
        different frequencies given the same atmospheric conditions. If given,
        it will assume that the phase rms at 10 km is 30 deg*scale_prms at 100
        GHz.
    """
    prms = 1.0/52.83 * antbl**0.8 # phase rms ~0.8 power to 10 km
    prms[antbl >= 1e4] = 3.0 * antbl[antbl >= 1e4]**0.25 # phase rms `0.25 power beyond 10 km
    if nu_scale != None:
        scale_prms *= (nu_scale/100.0e9)

    return prms * scale_prms


def baseline_arrayconf(fconf, refant = None):
    """Reads an array configuration file and computes antenna distances to a
    reference one.

    Args:
      fconf:
        Path to array configuration file (CASA format)
      refant:
        (optional) Index of antenna reference, with respect to the antennas
        in configuration file. If not set, it will use the antenna that is
        the farthest from the array center.

    Returns:
      Array with distance to a reference antenna for each antenna, and index of
      reference antenna used.

    """
    f = open(fconf)
    skip_header = 0
    while f.readline()[0] == '#':
        skip_header += 1
    f.close()
    arrayconf = np.genfromtxt(fconf, skip_header = skip_header, dtype = str)
    x = arrayconf[:,0].astype(float)
    y = arrayconf[:,1].astype(float)
    z = arrayconf[:,2].astype(float)
    # pad = arrayconf[:,4]

    if refant == None:
      refant = np.argmax(x**2. + y**2. + z**2.)

    baseline2ref = np.sqrt((x - x[refant])**2. + (y - y[refant])**2.
    + (z - z[refant])**2.)

    return baseline2ref, refant


def get_antpos_ms(visms):
    """Gets the X and Y positions of the antennas in a measurement set.

    Mostly copied from task_plotants.py.

    Args:
      visms:
        Path to measurement set.

    Returns:
      Arrays with X and Y positions of the antennas, in m, and with respect to
      the array center.
    """
    from casatools import table, quanta, measures, msmetadata
    tb = table( )
    me = measures( )
    qa = quanta( )

    # Getting array position in Wgs84
    metadata = msmetadata()
    metadata.open(visms)
    telescope = metadata.observatorynames()[0]
    arrayPos = metadata.observatoryposition()
    metadata.close()
    arrayWgs84 = me.measure(arrayPos, 'WGS84')
    arrayLon, arrayLat, arrayAlt = [arrayWgs84[i]['value']
            for i in ['m0','m1','m2']]

    # Getting antenna positions
    tb.open(visms + '/ANTENNA')
    antPositions = np.array([me.position('ITRF', qa.quantity(x, 'm'),
            qa.quantity(y, 'm'), qa.quantity(z, 'm'))
            for (x, y, z) in tb.getcol('POSITION').transpose()])
    tb.close()
    # Filtering for the ones that are used
    tb.open(visms)
    ants1 = tb.getcol('ANTENNA1')
    ants2 = tb.getcol('ANTENNA2')
    tb.close()
    antIdsUsed = list(set(np.append(ants1, ants2)))
    antPositions = [antPositions[i] for i in antIdsUsed]
    # Getting antenna positions in Wgs84
    antWgs84s = np.array([me.measure(pos, 'WGS84') for pos in antPositions])

    # Getting antenna positions from array center in m
    antLons, antLats = [np.array( [pos[i]['value']
        for pos in antWgs84s]) for i in ['m0','m1']]
    radE = 6370000.
    antXs = (antLons - arrayLon) * radE * np.cos(arrayLat)
    antYs = (antLats - arrayLat) * radE

    return antXs, antYs


def baseline_ms(visms, refant = None):
    """Reads antenna positions and returns distances to a reference one.

    Note that it only uses x and y coordinates (i.e., it ignores distance
    in z axis).

    Args:
      visms:
        Path to measurement set.
      refant:
        (optional) Index of antenna reference, with respect to the antennas
        in configuration file. If not set, it will use the antenna that is
        the farthest from the array center.

    Returns:
      Array with distance to a reference antenna for each antenna, and index of
      reference antenna used.
    """
    antXs, antYs = get_antpos_ms(visms)

    if refant == None:
        refant = np.argmax(antXs**2. + antYs**2.)

    refX, refY = antXs[refant], antYs[refant]

    baseline2ref = np.sqrt((antXs - refX)**2. + (antYs - refY)**2.)

    return baseline2ref, refant


def apply_phasenoise(visms, fconf = None, refant = None, prms_0 = 1.0,
                     nu_scale = None, prms_timescale = 1, apply = True,
                     verbose = False, plots = False,
                     rng = np.random.default_rng(2021)):
    """Computes and applies atmospheric phase noise to a measurement set.

    Creates a copy of input MS with corrupted phases.

    Args:
      visms:
        Path to measurement set
      fconf:
        (optional) Path to array configuration file (CASA format) in case user
        wants to use it to get antenna positions. If not given, it will read the
        antenna positions directly from the measurement set.
      refant:
        (optional) Index of antenna reference, with respect to the antennas
        in configuration file. If not set, it will use the antenna that is
        the farthest from the array center.
      prms_0:
        (Optional) Factor to determine the phase rms at 10 km. It can be given
        as a float, in which case it will work as a factor scaling the phase
        rms, with 1 being 30 degrees at 10 km. It can also be given as a string
        of the form 'XXdeg'. If set as a list, it will be interpreted as the
        scales for the different timescales at which phase noise will be
        introduced. If so, prms_timescale must also be set as a list. Default is
        1.0.
      nu_scale:
        (Optional) Frequency of the observation, in case one wants to simulate
        different frequencies given the same atmospheric conditions. If given,
        it will assume that the phase rms set by prms_0 is at 100 GHz.
      prms_timescale:
        (Optional) Timescale(s) at which phase rms will be introduced, in units
        of integration times. Default is 1 (i.e., phase noise will be introduced
        everye integration point). If set as a list, it will introduce phase
        noise multiple times at different timescales.
      apply:
        If True, it will apply the phase corruption. Default is True.
      verbose:
        If True, it will print to the screen some details. Default is False.
      plots:
        If True, it will plot the phase noise. Default is False.
      rng:
        Random number generator to be used (for reproducibility).
    """
    if fconf == None:
        baseline2ref, refant = baseline_ms(visms)
    else:
        baseline2ref, refant = baseline_arrayconf(fconf, refant)

    nants = len(baseline2ref)

    if type(prms_0) is list:
        if type(prms_0[0]) is str:
            prms_0 = list(map(lambda x: float(x.split('deg')[0]) /30.0, prms_0))
        # Scaling prms_0 to make phase rms consistent at the different timescales
        for i in range(len(prms_0) - 1, 0, -1):
            prms_0[i] = np.sqrt(prms_0[i]**2. - prms_0[i-1]**2.)
    elif type(prms_0) is str:
        prms_0 = float(prms_0.split('deg')[0]) / 30.0

    tb.close()
    os.system('rm -rf atmos.cal')
    gaincal(vis = visms,
            caltable = 'atmos.cal',
            refant = str(refant),
            minsnr = 0.001,
            solint = 'int',
            calmode = 'p')

    # Obtain certain 'atmos.cal' columns
    #   and fill in antenna-based antmospheric phases with
    tb.open('atmos.cal', nomodify=False)
    ants_cal = tb.getcol('ANTENNA1')
    time = tb.getcol('TIME')
    cparam = tb.getcol('CPARAM')
    npoints = len(ants_cal)
    nints = int(len(time) / nants)

    if type(prms_0) is not list:
        # If introducing phase noise in only one timescale
        # prms = phase_rms(baseline2ref[ants_cal], prms_0, nu_scale)
        # perror = rng.normal(0, prms)
        prms = phase_rms(baseline2ref, prms_0, nu_scale)
        if prms_timescale == 1:
            perror = rng.normal(0.0, prms[ants_cal])
        else:
            perror = np.zeros_like(ants_cal, dtype='float64')
            ntimesteps = round(nints / prms_timescale)
            for i in range(ntimesteps):
                perror_i = rng.normal(0.0, prms)
                start = i * prms_timescale * nants
                end = (i + 1) * prms_timescale * nants
                perror[start:end] = np.concatenate(prms_timescale * [perror_i])
    else:
        perror = np.zeros_like(ants_cal, dtype='float64')
        prms = []
        for tscale in range(len(prms_timescale)):
            prms.append(phase_rms(baseline2ref, prms_0[tscale], nu_scale))
            if prms_timescale[tscale] == 1:
                perror += rng.normal(0, prms[tscale][ants_cal])
            else:
                ntimesteps = round(nints / prms_timescale[tscale])
                for i in range(ntimesteps):
                    perror_i = rng.normal(0.0, prms[tscale])
                    start = i * prms_timescale[tscale] * nants
                    end = (i + 1) * prms_timescale[tscale] * nants
                    perror[start:end] += np.concatenate(prms_timescale[tscale]
                        * [perror_i])

    # PRINT OUT ANT PHASE ERRORS
    if verbose:
        print('     ant   time  baseline   rmsphase\n')
        for i in range(0,npoints):
            print('{} {} {:e} {:e} {:e} {:e}\n'.format(i, ants_cal[i],
                time[i]-time[0], baseline2ref[ants_cal[i]], prms[i], perror[i]))

    atmerror = ( np.cos(perror * np.pi / 180.0)
        + np.sin(perror * np.pi / 180.0)*1j )
    cparam[0][0] = atmerror  #X POL
    cparam[1][0] = atmerror  #Y POL  ASSUMED SAME

    tb.putcol('CPARAM', cparam)
    tb.flush()
    tb.close()

    # plot random atmospjajajheric phases
    if plots:
        # plotcal(caltable = 'atmos.cal',
        #     xaxis = 'time',
        #     yaxis = 'phase',
        #     iteration = 'antenna',
        #     plotrange = [0,0,-180,180],
        #     subplot = 331)
        fig = plt.figure(figsize=(10,6))
        fig.subplots_adjust(wspace=0.4, hspace=0.0)
        ax1 = fig.add_subplot(212)
        ax1.plot(baseline2ref[ants_cal], perror, '.', color='orange', alpha=0.2)
        ax1.set_xlabel('Baseline length [m]')
        ax1.set_ylabel('Phase noise [deg]')
        ax1.axhline(0.0, linestyle=':', color='black')
        ax1.set_ylim([-np.max(np.abs(perror)), np.max(np.abs(perror))])
        ax2 = fig.add_subplot(211)
        if type(prms_0) is not list:
            timescale = round(prms_timescale * (time[nants] - time[0]))
            ax2.plot(baseline2ref, prms[:len(baseline2ref)], '.',
                label = 'Timescale: {}s'.format(timescale))
        else:
            prms_plot = 0.0
            for tscale in range(len(prms_timescale)):
                timescale = round(prms_timescale[tscale] * (time[nants] -
                    time[0]))
                prms_plot = np.sqrt(prms_plot**2. +
                    prms[tscale][:len(baseline2ref)]**2.)
                ax2.plot(baseline2ref, prms_plot, '.',
                    label = 'Timescale: {}s'.format(timescale))
        ax2.legend()

        # ax2.set_xlabel('Baseline length [m]')
        ax2.set_xticklabels([])
        ax2.set_ylabel('Phase rms [deg]')
        ax2.axhline(0.0)
        plt.savefig('{}_phase_errors.pdf'.format(os.path.basename(visms)[:-3]),
                    dpi=350 , bbox_inches='tight')
        plt.close(fig)

    #apply the atmospheric phases to uvdata
    if apply:
        os.system('rm -rf {}'.format(visms[:-3]+'_phasenoisy.ms'))
        os.system('cp -r {} {}'.format(visms, visms[:-3]+'_noisy_temp.ms'))
        applycal(vis = visms[:-3]+'_noisy_temp.ms',
                 gaintable = 'atmos.cal')
        split(vis = visms[:-3] + '_noisy_temp.ms',
              outputvis = visms[:-3] + '_phasenoisy.ms')
        os.system('rm -rf {}'.format(visms[:-3]+'_noisy_temp.ms'))


def simobserve_phase(fconf, pwv, obstime, model_image = None,
                     phase_center = None, comps = None, nu_comp = None,
                     mapsize = "", bandwidth = 8, prms_0 = 1.0, nu_scale = None,
                     integration = '40s', prms_timescale = 1, verbose = False,
                     plots = False, seed = 2021, overwrite = False, **kwargs):
    """Function to simulate ALMA observations, including phase atmospheric
    noise.

    It can be called using a model image in FITS format (giving the path to
    the image in the model_image argument), OR using the arguments phase_center,
    comps, and nu to construct a model image from a list of Gaussian components.

    Args:
      model_image:
        Name of model image in FITS format.
      phase_center:
        List or tuple with RA (in hours) and DEC (in degrees) of the phase
        center. Example:
          phase_center = (15.8, -35.0)
      comps:
        List with Gaussian components. Each item of the list should be another
        list with 6 items: flux, beam major axis, beam minor axis, beam PA,
        offset in RA, and offset in DEC. Example:
          comps = [[0.1,  0.001,  0.001,    0.0,  0.00,   0.00],
                   [1.0,  0.1,    0.01,     0.5,  1.0,    -1.0]]
      nu_comp:
        Frequency (in Hz) of the emission of the model.
      mapsize:
        Angular size of model map, defined as 'XXarcsec'. Default is "", which
        simobserve will interpret as to span the full model.
      fconf:
        Path to array configuration file (CASA format).
      pwv:
        Precipitable water vapor to be used with simobserve in CASA.
      obstime:
        Total observing time to be used with simobserve in CASA.
      bandwidth:
        Bandwidth of model (either FITS file or component list), in GHz.
      prms_0:
        (Optional) Factor to determine the phase rms at 10 km. It can be given
        as a float, in which case it will work as a factor scaling the phase
        rms, with 1 being 30 degrees at 10 km. It can also be given as a string
        of the form 'XXdeg'. If set as a list, it will be interpreted as the
        scales for the different timescales at which phase noise will be
        introduced. If so, prms_timescale must also be set as a list. Default is
        1.0.
      nu_scale:
        (Optional) Frequency of the observation, in case one wants to simulate
        different frequencies given the same atmospheric conditions. If given,
        it will assume that the phase rms set by prms_0 is at 100 GHz.
      integration:
        (Optional) Integration time used during the simulation, in sec. Default
        is 40s.
      prms_timescale:
        (Optional) Timescale(s) at which phase rms will be introduced. It can be
        given in 'XXs' format, or in units of integration times. Default is 1
        (i.e., phase noise will be introduced everye integration point). If set
        as a list, it will introduce phase noise multiple times at different
        timescales.
      verbose:
        If True, it will print to the screen some details. Default is False.
      plots:
        If True, it will produce a few plots. To do so, it will need the
        vis_tools package (https://github.com/emaciasq/vis_tools). Default is
        False.
      overwrite:
        If True it will overwrite any output files already in the folder.
        Default is False.
      seed:
        Seed for random number generator (for reproducibility).
    """
    rng = np.random.default_rng(seed)
    if 'project' in kwargs:
        project = kwargs['project']
    else:
        project = model_image[:-5] # removing the .fits at the end

    if type(prms_0) is list:
        if type(prms_timescale) is not list:
            raise IOError('If multiple scale factors are provided, then '+
                          'multiple phase rms timescales must be provided too.')
        if type(prms_timescale[0]) is str:
            integration_int = int(integration[:-1])
            for i in range(len(prms_timescale)):
                t_i = int(prms_timescale[i][:-1])
                prms_timescale[i] = int(t_i / integration_int)
    else:
        if type(prms_timescale) is list:
            raise IOError('If multiple phase rms timescales are provided, then'+
                          ' multiple scale factors must also be provided.')
        if type(prms_timescale) is str:
            integration_int = int(integration[:-1])
            t_i = int(prms_timescale[:-1])
            prms_timescale = int(t_i / integration_int)


    if type(bandwidth) is not str:
        bandwidth = '{}GHz'.format(bandwidth)

    if model_image != None:
        if (phase_center != None) or (comps != None) or (nu_comp != None):
            print('WARNING: A model image was provided together with a list of'
                  ' components. The code will use the model image and disregard'
                  ' the list of components.')
        if (os.path.exists(project) == False) or (overwrite == True):
            simobserve(project = project,
                    antennalist = fconf,
                    skymodel = model_image,
                    inwidth = bandwidth,
                    user_pwv = pwv,
                    seed = seed,
                    integration = integration,
                    totaltime = obstime)
        elif overwrite == False:
            print('Simulated project exists. Will skip simobserve step')

    elif (phase_center != None) and (comps != None) and (nu_comp != None):
        if (os.path.exists(project) == False) or (overwrite == True):
            gaussian_comp_model(phase_center, comps, nu_comp,
                name = project + '.cl')
            simobserve(project = project,
                    antennalist = fconf,
                    complist = project + '.cl',
                    compwidth = bandwidth,
                    direction = coords_to_string(phase_center),
                    mapsize = mapsize,
                    user_pwv = pwv,
                    seed = seed,
                    integration = integration,
                    totaltime = obstime)
        elif overwrite == False:
            print('Simulated project exists. Will skip simobserve step')
    else:
        raise IOError('Input model not correctly set. Input model must be '
            'either a model image in FITS format (giving the path to the image '
            'in the model_image argument), OR a component list set using the '
            'arguments phase_center, comps, AND nu.')

    ms_thermalnoise = '{}/{}.{}.noisy.ms'.format(
        project, project, os.path.basename(fconf)[:-4])

    apply_phasenoise(ms_thermalnoise, prms_0 = prms_0, nu_scale = nu_scale,
                     prms_timescale = prms_timescale, apply = True,
                     verbose = verbose, plots = plots, rng = rng)

    ms_phasenoise = '{}/{}.{}.noisy_phasenoisy.ms'.format(
        project, project, os.path.basename(fconf)[:-4])

    if plots:
        try:
            import vis_tools
            execfile('{}/CASA_vis_tools.py'.format(
                os.path.dirname(vis_tools.__file__)))
        except ModuleNotFoundError:
            print('Cannot find vis_tools in your path. The module vis_tools ' +
            'and accompanying script CASA_vis_tools.py is needed to run ' +
            'these plots. You can find them at ' +
            'https://github.com/emaciasq/vis_tools')
        except FileNotFoundError:
            print('Cannot find CASA_vis_tools.py, which is needed to run ' +
            'these plots. Please put the script CASA_vis_tools.py in the same' +
            ' path where vis_tools.py is located. You can find both at ' +
            'https://github.com/emaciasq/vis_tools')

        visobs = get_vis_obs(ms_thermalnoise)
        visobs_pnoisy = get_vis_obs(ms_phasenoise)

        # Histogram of baselines
        fig = plt.figure(figsize=(8,4))
        ax = fig.add_subplot(111)
        ax.hist(visobs.uvwave/1000., bins=20)
        ax.set_xlabel(r'Baseline length [k$\lambda$]')
        ax.set_ylabel(r'N')
        plt.savefig(
            '{}_baselines.pdf'.format(os.path.basename(ms_thermalnoise)[:-3]),
            dpi=350, bbox_inches='tight')
        plt.close(fig)

        visobs.deproject(0.0, 0.0)
        visobs.bin_vis(nbins=150)
        visobs_pnoisy.deproject(0.0, 0.0)
        visobs_pnoisy.bin_vis(nbins=150)

        # Visibility plot
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111)
        # ax.plot(
        #   visobs_pnoisy.rho/1000., visobs_pnoisy.r, '.', markeredgecolor=None,
        #   color='orange', alpha=0.1, markersize=2, label='With phase noise')
        # ax.plot(
        #   visobs.rho/1000., visobs.r, '.', markeredgecolor=None, color='blue',
        #   alpha=0.1, markersize=2, label='Without phase noise')
        ax.errorbar(
            visobs_pnoisy.bin_centers/1000., visobs_pnoisy.r_binned,
            yerr = visobs_pnoisy.r_sigma, color='orange', markersize=2,
            label='With phase noise')
        ax.errorbar(
            visobs.bin_centers/1000., visobs.r_binned, yerr = visobs.r_sigma,
            color='blue', markersize=2, label='Without phase noise')
        ax.set_ylabel('Re [Jy]')
        ax.set_xlabel(r'$\rho$ [k$\lambda$]')
        ax.legend(fancybox=True)
        plt.savefig(
            '{}_visibilities.pdf'.format(
            os.path.basename(ms_thermalnoise)[:-3]), dpi=350 ,
            bbox_inches='tight')


