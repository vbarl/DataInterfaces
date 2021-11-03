#
# Alternative functions to what is inside assp.py
#
2021-11-03  Vasileios Barlakas <vasileios.barlakas@chalmers.se>
        * DataInterfaces/Python/assp.py:
          * ssdb2assp:
          Added functionality for retrieving the asymmetry parameter in case
          of ARO ice and TRO liquid and ice hydrometeors.
          * assp_import_ssdb:
          Added functionality for retrieving the asymmetry parameter in case
          of ARO ice and TRO liquid and ice hydrometeors.
        * DataInterfaces/Python/rttov.py:
          * get_assp:
          Added functionality for retrieving the asymmetry parameter. Yet,
          this is not returned.
        * demo_ssp4arts.py:
        Updated, due to extra output (asymmetry parameter)

        * DataInterfaces/Python/assp.py:
          * ssdb2asspAndG:
          Equivalent to ssdb2assp, but inluding the asymmetry parameter as an extra 
          output in case of oriented ice hydrometeors
        * DataInterfaces/Python/assp.py:
          * get_full_assp:     
          Equivalent to assp_import_ssdb, but inluding the asymmetry parameter as an 
          extra output  
          * ssdb2assp_full:
          Added functionality for retrieving the asymmetry parameter in case
          of ARO ice and TRO liquid and ice hydrometeors.
          
        * DataInterfaces/Python/assp.py:
          * calc_g4aro:
          Reads data from SSDB for one habit (only oriented ice) and derives the
          asymmetry parameter.
###################################################################################
###################################################################################
#
# Externally deriving the asymmetry parameter
# 
# FOR utils
#
def calc_g4aro(habit_id, orientation,
                     habit_folder=None,
                     size_range=None, size_type='dveq',
                     freq_range=None, temp_range=None,
                     allow_nodata=False):
  #2021-11-03 Vasileios Barlakas: Derives the asymmetry parameter 
  #                               for ARO ice hydrometeors following
  #                               M. Brath's draft.                          
  """Reads data from SSDB database for one habit (oriented ice) and derives
  the asymmetry parameter g.
  
  The reading can be restricted in terms of size, frequency
  and temperature. For frequency and temperature, data are selected as: 
    limit1 <= data <= limit2
  while for frequency data at the lower limit is excluded:
    limit1 < data <= limit2.

  Default is to issue an error as soon as data are missing for a frequency
  and temperature combination. Allow missing data by setting the optional
  argument *allow_nodata* to true. Note that sets of frequencies and
  temperatures are allowed to differ between sizes, independently of how
  *allow_nodata* is set. For example, for one size there can be a single
  temperature, while other sizes have a temperature grid with several
  elements.
  
  Parameters
  ----------
  habit_id: int
    Habit id number.
  orientation: str
    Descriptor of orientation to explore for the chosen habit.
  habit_folder: str
    Full path to a habit folder.
    Temporary add for testing az.random not-yet-in-DB-structure. If used,
    provide a dummy habit_id.
  allow_nodata: bool
    See above. Default is false.
  size_range: 2-element list
    Particle size limits [unit: m]. Ignore data outside these limits. If not
    given, no limits are applied, i.e. all available data is considered.
  size_type: str
    Quantity used for size cropping. Allowed options are 'dveq' (default),
    'dmax' and 'mass'.
  freq_range: 2-element list
    Frequency limits [unit: Hz]. Ignore data outside these limits. If not
    given, no limits are applied, i.e. all available data is considered.
  temp range: 2-element list
    Temperature limits [unit: K]. Ignore data outside these limits. If not
    given, no limits are applied, i.e. all available data is considered.

  Returns
  -------
  g: numpy array
    Array of asymmetry parameter values
  """
  try:
    from . import sph
  except ImportError:
    try: #relative import does not work from scripts in the package folder
      import sph
      # from Python import sph
    except ImportError:
      raise Exception( 'Module *sph* of database interface required but not found.' )
  except Exception as e:
      print('Module *assp* requires module *sph*, but failed to import it:\n%s' %str(e))

  if (habit_folder is None):
    data = utils.ssdb_import_habit( habit_id, orientation,
                                  allow_nodata=allow_nodata,
                                  size_range=size_range, size_type=size_type,
                                  freq_range=freq_range, temp_range=temp_range )
  else:
    data = utils.ssdb_import_habit( habit_id, orientation,
                                    habit_folder=habit_folder,
                                  allow_nodata=allow_nodata,
                                  size_range=size_range, size_type=size_type,
                                  freq_range=freq_range, temp_range=temp_range )


  grid_size_theta=90

  nf = data[0]['freq'].size
  nt = data[0]['temp'].size
  theta_inc   = data[0]['SSD'][0,0]['SingleScatteringData']['za_inc']['data']

  g    = np.empty((len(data),nf,nt,len(theta_inc)))*np.NAN
  g_i  = np.empty((len(theta_inc)))*np.NAN
  
  
  # Loop particle sizes and retrive data
  for isize in np.arange(len(data)):     # size loop
    
    try:
     ssd = data[isize]['SSD']
     for f in np.arange(nf):              # freq loop
      for t in np.arange(nt):             # temp loop
        phamat_real = ssd[f,t]['SingleScatteringData']['phaMat_data_real']['data'][0,0,:,:]
        phamat_imag = ssd[f,t]['SingleScatteringData']['phaMat_data_imag']['data'][0,0,:,:]
      
        ph_sph      = np.complex128(phamat_real+1j*phamat_imag)
      
        for i in range(len(theta_inc)):
          # get truncation level of spherical harmonics series
          max_ntrunc = sph.get_lmax(ph_sph[i,:])
          if max_ntrunc % 2 ==0:
             min_grid_size=2*max_ntrunc+2
          else:
             min_grid_size=2*(max_ntrunc+1) 
          if min_grid_size<grid_size_theta:
             min_grid_size=grid_size_theta 
          elif grid_size_theta<min_grid_size:
             grid_size_theta=min_grid_size       
          if min_grid_size % 2 == 1:
             min_grid_size=min_grid_size+1
          # set up spherical harmonics object
          sph_object = sph.Spharmt(min_grid_size*2,min_grid_size,max_ntrunc,1,gridtype='regular')
          
          # rotate spherical harmonics series, so that forward direction points toward northpole
          # and backward direction towards southpole.
          phase_sph_rot = sph_object.yrotate(ph_sph[i,:], -theta_inc[i]*np.pi/180)
          g_i[i]        = np.real(phase_sph_rot[1]/phase_sph_rot[0]/np.sqrt(4*np.pi)*2)
        g[isize,f,t,...]              = np.array(g_i) 
    except:
     print('For particle #%i, no data available -'
              ' skipping this particle (ie. no equivalent entry in S & M).\n'
              %(isize))
  return g
###################################################################################
###################################################################################
#
# ONLY for ICE ARO provides the asymmetry
#
def ssdb2asspAndG(SSD, freq, temp, nodata,interpm='linear'):
  #2017-03-22 Jana Mendrok
  #2017-11-28 Robin Ekelund: Fixed sph importing error
  #2021-02-19 Vasileios Barlakas: Extended towards liquid oriented hydrometeors
  #2021-11-03 Vasileios Barlakas: Added functionality for deriving the asymmetry
  #                               parameter for ARO ice hydrometeors following
  #                               M. Brath's draft.
  """Equivalent to ssdb2assp's functionality, but it additionally returns the
  asymmetry parameter for azimuthally oriented ice hydrometeors.

  The input SSD data should be imported through utils.ssdb_import_data using
  the grid_sort option.

  Parameters
  ----------
  SSD: 2-D numpy array of dictionary (dim: [freq,temp])
    Single scattering data over frequencies and temperatures.
  freq: 1-D numpy array
    Frequency grid of the SSD data.
  temp: 1-D numpy array
    Temperature grid of the SSD data.
  nodata: numpy array
    Boolean matrix, flagging empty SSD elements.
  interpm: str
    Method for zenith angle interpolation.

    Allowed are:
     TRO: 'linear', 'nearest', 'zero', 'slinear', 'quadratic','cubic'
     from scipy's interp1d (where the latter four use spline interpolation of
     zeroth to third order, respectively) as well as 'pchip' applying scipy's
     PchipInterpolator (supposed to provide similar, but not necessarily
     identical results to Matlab's pchip interpolation from interp1).
     ARO (liquid only): 'linear' and 'nearest' on the basis of RegularGridInterpolator
    Default is 'linear'.

  Returns
  -------
  S: object (typhon.arts.scattering.SingleScatteringData
    ARTS-type single scattering data of a specific scattering element (SingleScatteringData).
  M: object (typhon.arts.scattering.ScatteringMetaData)
    ARTS-type scattering meta data of a specific scattering element (ScatteringMetaData).
  G: List containing the corresponding asymmetry parameter values
  -------
  TODOS:
      * Add functionality for oriented liquid hydrometeors
      * Add functionality for totally randomly oriented hydrometeors (rttov.assp2g)
  -------
  """
  # Basic sanity check of input
  if not( SSD.shape[0]==freq.size ):
    raise Exception( 'Mismatch in size between *SSD* and *freq*.' )
  if not( SSD.shape[1]==temp.size ):
    raise Exception( 'Mismatch in size between *SSD* and *temp*.' )
  if not( all('SingleScatteringData' in d for d in SSD.reshape(-1)) and
          all('ShapeData' in d for d in SSD.reshape(-1)) ):
    raise Exception( \
      'At least one element of *SSD* does not seem to hold the expected data format.' )

  # Set M
  #  (here from the first element in SSD; later we check that the crucial meta
  #   data is consistent between all SSD)
  if use_typhon:
    M = tas.ScatteringMetaData()
    if (SSD.size>0):
      M.description         = 'Meta data for '+SSD[0,0]['ShapeData']['description']
      M.source              = SSD[0,0]['ShapeData']['source']
      M.refr_index          = SSD[0,0]['ShapeData']['refrIndex_model']
      assert( (len(SSD[0,0]['ShapeData']['mass']['dim'])==0) and
              (SSD[0,0]['ShapeData']['mass']['data'].size==1) ), \
        'Mass data has wrong dimensions.'
      assert( (len(SSD[0,0]['ShapeData']['diameter_max']['dim'])==0) and
              (SSD[0,0]['ShapeData']['diameter_max']['data'].size==1) ), \
        'Max. diameter data has wrong dimensions.'
      assert( (len(SSD[0,0]['ShapeData']['diameter_vol_eq']['dim'])==0) and
              (SSD[0,0]['ShapeData']['diameter_vol_eq']['data'].size==1) ), \
        'Vol. equ. diameter data has wrong dimensions.'
      assert( (len(SSD[0,0]['ShapeData']['diameter_area_eq_aerodynamical']['dim'])==0) and
              (SSD[0,0]['ShapeData']['diameter_area_eq_aerodynamical']['data'].size==1) ), \
        'Aerodyn. area equ. diameter data has wrong dimensions.'
      M.mass                = np.float(SSD[0,0]['ShapeData']['mass']['data'])
      M.diameter_max        = np.float(SSD[0,0]['ShapeData']['diameter_max']['data'])
      M.diameter_volume_equ = np.float(SSD[0,0]['ShapeData']['diameter_vol_eq']['data'])
      M.diameter_area_equ_aerodynamical \
                            = np.float(SSD[0,0]['ShapeData']['diameter_area_eq_aerodynamical']['data'])
  else:
    raise Exception( 'Non-typhon *ScatteringMetaData* not yet implemented.' )

  phase = SSD[0,0]['ShapeData']['phase']

  # Basic data of S
  if use_typhon:
    S = tas.SingleScatteringData()
    S.version     = 3
    if (SSD.size>0):
      try:
        ptypeID = db_ptype[SSD[0,0]['SingleScatteringData']['orient_type']]
      except KeyError:
        raise Exception( \
          "Database ptype '%s' is unknown." %SSD[0,0]['SingleScatteringData']['orient_type'] )
      try:
        S.ptype       = a_ptype[ptypeID]
      except KeyError:
        raise Exception( \
          "No ARTS ptype equivalent defined for database ptype '%s' (internal ID=%i)." \
          %(SSD[0,0]['SingleScatteringData']['orient_type'],ptypeID) )
      S.description = SSD[0,0]['ShapeData']['description']
      S.f_grid  = freq
      S.T_grid  = temp
  else:
    raise Exception( 'Non-typhon *SingleScatteringData* not yet implemented' )

  # Fill S according to ptype
  nf = freq.size
  nt = temp.size

  if (SSD.size>0):

    #####
    # totally random orientation
    #####
    if ( S.ptype==a_ptype[20] ):
      raise Exception( 'This hydrometeor is characterized by total random orientation' )
    #####
    # azimuthally random orientation
    #####
    elif ( S.ptype==a_ptype[30] ):
      #####
      # Separate according to phase: liquid (gridded data) vs ice (spherical harmonics)
      #####
      if phase=='ice':

          try:
            from . import sph
          except ImportError:
            try: #relative import does not work from scripts in the package folder
              import sph
              # from Python import sph
            except ImportError:
              raise Exception( 'Module *sph* of database interface required but not found.' )
          except Exception as e:
              print('Module *assp* requires module *sph*, but failed to import it:\n%s' %str(e))

          S.za_grid=SSD[0,0]['SingleScatteringData']['za_inc']['data']

          for f in np.arange(freq.size):
            for t in np.arange(temp.size):
              if (not nodata[f,t]):
                # A number of asserts to double-check the input
                assert( abs(M.mass-SSD[f,t]['ShapeData']['mass']['data'])<1e-12 ), \
                  'Mass info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
                  %(freq[f]*1e-9,temp[t])
                assert( abs(M.diameter_max-SSD[f,t]['ShapeData']['diameter_max']['data'])<1e-12 ), \
                  'Dmax info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
                  %(freq[f]*1e-9,temp[t])
                assert( abs(M.diameter_volume_equ-SSD[f,t]['ShapeData']['diameter_vol_eq']['data'])<1e-12 ), \
                  'Dveq info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
                  %(freq[f]*1e-9,temp[t])

                try:
                  ptypeID = db_ptype[SSD[f,t]['SingleScatteringData']['orient_type']]
                except KeyError:
                  raise Exception( \
                    "Database ptype for SSD at f=%.1fGHz and T=%.1fK ('%s') is unknown." \
                    %(freq[f]*1e-9,temp[t],SSD[f,t]['SingleScatteringData']['orient_type']) )
                try:
                  ptypeName = a_ptype[ptypeID]
                except KeyError:
                  raise Exception( \
                    "No ARTS ptype equivalent defined for database ptype '%s' (internal ID=%i) for SSD at f=%.1fGHz and T=%.1fK." \
                    %(SSD[f,t]['SingleScatteringData']['orient_type'],ptypeID,freq[f]*1e-9,temp[t]) )
                assert( S.ptype==ptypeName ), \
                  'ptype for SSD at f=%.1fGHz and T=%.1fK inconsistent with initial value.' \
                  %(freq[f]*1e-9,temp[t])

                if not( abs(freq[f]-SSD[f,t]['SingleScatteringData']['frequency']['data'])<1e3 ):
                  raise Exception( \
                    'Freq info in SSD at f=%.1fGHz and T=%.1fK inconsistent with freq grid value.' \
                    %(freq[f]*1e-9,temp[t]) )
                if not( abs(temp[t]-SSD[f,t]['SingleScatteringData']['temperature']['data'])<0.001 ):
                  raise Exception( \
                    'Temp info in SSD at f=%.1fGHz and T=%.1fK inconsistent with temp grid value.' \
                    %(freq[f]*1e-9,temp[t]) )
                assert( (SSD[f,t]['SingleScatteringData']['za_inc']['data'][0]==0.) and
                        (SSD[f,t]['SingleScatteringData']['za_inc']['data'][-1]==180.) ), \
                  'za_inc at f=%.1fGHz and T=%.1fK does not fulfill basic requirement\n' + \
                  'of first value being at 0deg and/or last value at 180deg.' %(freq[f]*1e-9,temp[t])

                # for now we require all data to be on identical za_inc grids
                assert( (S.za_grid==SSD[f,t]['SingleScatteringData']['za_inc']['data']).all() ), \
                  'za_inc grid f=%.1fGHz and T=%.1fK inconsistent with global za_grid.' \
                  %(freq[f]*1e-9,temp[t])

          pha_inds = SSD[0,0]['SingleScatteringData']['phaMat_index']['data']
          ext_inds = SSD[0,0]['SingleScatteringData']['extMat_index']['data']
          abs_inds = SSD[0,0]['SingleScatteringData']['absVec_index']['data']

          # Set scattering grid sizes
          #   ARTS requires za_inc==za_sca, sph requires (as implemented)
          #    aa_sca[0-180]=za_sca. Hence here we fix everything to the za_inc, which
          #    is the coarser grid in our ADDA calcs.
          #   Basically, that could be replaced by two (independent) input parameters,
          #    one for za, one for aa. However, that requires more interpolations
          #    (and rotations?) in sph.
          nza = len(SSD[0,0]['SingleScatteringData']['za_inc']['data'])
          g    = np.empty((nf,nt,len(S.za_grid)))*np.NAN
          g_i  = np.empty((len(S.za_grid)))*np.NAN
          got_grid=False
          for f in np.arange(freq.size):
            for t in np.arange(temp.size):
              if (not nodata[f,t]):

                # extract optical properties
                ext_mat = SSD[f,t]['SingleScatteringData']['extMat_data']['data']
                abs_vec = SSD[f,t]['SingleScatteringData']['absVec_data']['data']

                phamat_real = SSD[f,t]['SingleScatteringData']['phaMat_data_real']['data']
                phamat_imag = SSD[f,t]['SingleScatteringData']['phaMat_data_imag']['data']
                phamat_sph = phamat_real+1j*phamat_imag
               
                # get truncation level
                lmax = sph.get_lmax(SSD[f,t]['SingleScatteringData']['sph_coeffs_scat']['data'])
                # transform spherical harmonics to grid
                [phase_matrix_reg,theta_s,phi_s] = \
                  sph.create_regular_representation(phamat_sph, lmax,
                                                    grid_size_theta=nza, gridtype='regular')
                  
                ###############################################################                                    
                # Settings for asymetry parameter
                ph_real     = phamat_real[0,0,:,:]
                ph_imag     = phamat_imag[0,0,:,:]
                ph_sph      = np.complex128(ph_real+1j*ph_imag)
                
                theta_inc = S.za_grid
                grid_size_theta=90
                
                for i in range(len(theta_inc)):
                  # get truncation level of spherical harmonics series
                  max_ntrunc = sph.get_lmax(ph_sph[i,:])
                  if max_ntrunc % 2 ==0:
                      min_grid_size=2*lmax+2
                  else:
                      min_grid_size=2*(lmax+1)  
                  if min_grid_size<grid_size_theta:
                      min_grid_size=grid_size_theta
                  elif grid_size_theta<min_grid_size:
                      grid_size_theta=min_grid_size
                  if min_grid_size % 2 == 1:
                      min_grid_size=min_grid_size+1
                  # set up spherical harmonics object
                  sph_object = sph.Spharmt(min_grid_size*2,min_grid_size,max_ntrunc,1,gridtype='regular')
                  # rotate spherical harmonics series, so that forward direction points toward northpole
                  # and backward direction towards southpole.
                  phase_sph_rot = sph_object.yrotate(ph_sph[i,:], -theta_inc[i]*np.pi/180)
                  g_i[i]        = np.real(phase_sph_rot[1]/phase_sph_rot[0]/np.sqrt(4*np.pi)*2)
                ###############################################################
                assert( (S.za_grid==theta_s).all() ), \
                  'Spherical harmonics output polar grid at f=%.1fGHz and T=%.1fK' + \
                  ' inconsistent with global za_grid.' %(freq[f]*1e-9,temp[t])

                # allocate and set angle grids
                if got_grid==False:
                  S.aa_grid=phi_s
                  S.ext_mat_data = np.empty((nf,nt,len(S.za_grid),1,3))*np.NAN
                  S.abs_vec_data = np.empty((nf,nt,len(S.za_grid),1,2))*np.NAN
                  S.pha_mat_data = np.empty((nf,nt,len(S.za_grid),len(S.aa_grid),
                                            len(S.za_grid),1,16))*np.NAN
                  got_grid=True

                # Fill data fields
                S.ext_mat_data[f,t,...] = np.transpose(ext_mat,(2,1,0))
                S.abs_vec_data[f,t,...] = np.transpose(abs_vec,(2,1,0))
                S.pha_mat_data[f,t,...] = np.transpose(phase_matrix_reg,(3,4,2,1,0))
                g[f,t,...]              = np.array(g_i) # VB
      elif phase=='liquid':
          raise Exception( 'This is a liquid hydrometeor with azimuth random orientation' )
    #####
    # unknown particle type
    #####
    else:
      raise Exception( "Ptype '%s' not (yet?) implemented." %S.ptype )

  else: # if (SSD.size>0)
    print('No SSD data available.\n'
          'Returning default-filled SingleScatteringData and ScatteringMetaData'
          ' (grids partly empty, data fields all empty).')

  return S, M, g
  
###################################################################################
#
# Alternative to assp_import_ssdb, with both options for importing the data
#
def get_full_assp(habit_id, orientation,
                     habit_folder=None,
                     size_range=None, size_type='dveq',
                     freq_range=None, temp_range=None,
                     allow_nodata=False):
  # 2021-11-03 Vasileios Barlakas: Added functionality for azimuthally random 
  #                                orientation following J. Mendrok's draft.
  """
  Equivalent to assp_import_ssdb's functionality, but additionally returns the 
  asymmetry parameter for azimuthally oriented ice hydrometeors.
  ----------
  Parameters
  ----------
  habit_id: int
    Habit id number.
  orientation: str
    Descriptor of orientation to explore for the chosen habit.
  habit_folder: str
    Full path to a habit folder.
    Temporary add for testing az.random not-yet-in-DB-structure. If used,
    provide a dummy habit_id.
  allow_nodata: bool
    See above. Default is false.
  size_range: 2-element list
    Particle size limits [unit: m]. Ignore data outside these limits. If not
    given, no limits are applied, i.e. all available data is considered.
  size_type: str
    Quantity used for size cropping. Allowed options are 'dveq' (default),
    'dmax' and 'mass'.
  freq_range: 2-element list
    Frequency limits [unit: Hz]. Ignore data outside these limits. If not
    given, no limits are applied, i.e. all available data is considered.
  temp range: 2-element list
    Temperature limits [unit: K]. Ignore data outside these limits. If not
    given, no limits are applied, i.e. all available data is considered.

  Returns
  -------
  S: list of objects (of SingleScatteringData)
    Array of ARTS-type single scattering data, one entry per particle size.
  M: list of objects (of ScatteringMetaData)
    Array of ARTS-type scattering meta data, one entry per particle size.
  C: List of asymmetry parameter
  ----------
  """
  # import data
  if (habit_folder is None):
    data = utils.ssdb_import_habit( habit_id, orientation,
                                  allow_nodata=allow_nodata,
                                  size_range=size_range, size_type=size_type,
                                  freq_range=freq_range, temp_range=temp_range )
  else:
    data = utils.ssdb_import_habit( habit_id, orientation,
                                    habit_folder=habit_folder,
                                  allow_nodata=allow_nodata,
                                  size_range=size_range, size_type=size_type,
                                  freq_range=freq_range, temp_range=temp_range )

  # Loop particle sizes and convert data
  S = []
  M = []
  G = []
  print( 'Converting extracted SSDB data to ARTS format data.' )
  print( '%i scattering elements to process' %len(data) )
  for i in np.arange(len(data)):
    if ((i+1)%5==0):
      print( '  processing element %i' %(i+1) )
    try:
      if (data[i]['SSD'].size>0):
      
        # g for only ARO ice
        #S1,M1,G1 = ssdb2asspAndG( data[i]['SSD'], data[i]['freq'],
        #                   data[i]['temp'], data[i]['nodata'] )
                           
        # g for both TRO and ARO (not ARO liquid)
        S1,M1,G1 = ssdb2assp_full( data[i]['SSD'], data[i]['freq'],
                           data[i]['temp'], data[i]['nodata'] )
        S.append( S1 )
        M.append( M1 )
        G.append( G1 )
      else:
        try:
          print('For this particle (%s=%.0f um), no SSD available -'
                ' skipping this particle (ie. no equivalent entry in S & M).\n'
                %(size_type,data[i][size_type]*1e6))
        except:
          print('For particle #%i, no data available -'
                ' skipping this particle (ie. no equivalent entry in S & M).\n'
                %(i))
    except KeyError:
        print('For particle #%i, no data available -'
              ' skipping this particle (ie. no equivalent entry in S & M).\n'
              %(i))

  return S, M, G
###################################################################################
def ssdb2assp_full(SSD, freq, temp, nodata,interpm='linear'):
  #2017-03-22 Jana Mendrok
  #2017-11-28 Robin Ekelund: Fixed sph importing error
  #2021-02-19 Vasileios Barlakas: Extended towards liquid oriented hydrometeors
  #2021-11-03 Vasileios Barlakas: Added functionality for deriving the asymmetry
  #                               parameter following M. Brath's draft.
  """Converts internal SSD to ARTS single scattering properties. It additionally
     provides the asymmetry parameter as an extra output.

  The input SSD data should be imported through utils.ssdb_import_data using
  the grid_sort option.

  Parameters
  ----------
  SSD: 2-D numpy array of dictionary (dim: [freq,temp])
    Single scattering data over frequencies and temperatures.
  freq: 1-D numpy array
    Frequency grid of the SSD data.
  temp: 1-D numpy array
    Temperature grid of the SSD data.
  nodata: numpy array
    Boolean matrix, flagging empty SSD elements.
  interpm: str
    Method for zenith angle interpolation.

    Allowed are:
      TRO: 'linear', 'nearest', 'zero', 'slinear', 'quadratic','cubic'
      from scipy's interp1d (where the latter four use spline interpolation of
      zeroth to third order, respectively) as well as 'pchip' applying scipy's
      PchipInterpolator (supposed to provide similar, but not necessarily
      identical results to Matlab's pchip interpolation from interp1).
      ARO (liquid only): 'linear' and 'nearest' on the basis of RegularGridInterpolator
    Default is 'linear'.

  Returns
  -------
  S: list of objects (of SingleScatteringData)
    Array of ARTS-type single scattering data, one entry per particle size.
  M: list of objects (of ScatteringMetaData)
    Array of ARTS-type scattering meta data, one entry per particle size.
  g: List containing the corresponding asymmetry parameter values.
  ----------
  TODOS:
      * Add functionality for oriented liquid hydrometeors
  """
  #####
  #FIXME: implement non-typhon version!
  #####
  
  import rttov
  
  # Basic sanity check of input
  if not( SSD.shape[0]==freq.size ):
    raise Exception( 'Mismatch in size between *SSD* and *freq*.' )
  if not( SSD.shape[1]==temp.size ):
    raise Exception( 'Mismatch in size between *SSD* and *temp*.' )
  if not( all('SingleScatteringData' in d for d in SSD.reshape(-1)) and
          all('ShapeData' in d for d in SSD.reshape(-1)) ):
    raise Exception( \
      'At least one element of *SSD* does not seem to hold the expected data format.' )

  # Set M
  #  (here from the first element in SSD; later we check that the crucial meta
  #   data is consistent between all SSD)
  if use_typhon:
    M = tas.ScatteringMetaData()
    if (SSD.size>0):
      M.description         = 'Meta data for '+SSD[0,0]['ShapeData']['description']
      M.source              = SSD[0,0]['ShapeData']['source']
      M.refr_index          = SSD[0,0]['ShapeData']['refrIndex_model']
      assert( (len(SSD[0,0]['ShapeData']['mass']['dim'])==0) and
              (SSD[0,0]['ShapeData']['mass']['data'].size==1) ), \
        'Mass data has wrong dimensions.'
      assert( (len(SSD[0,0]['ShapeData']['diameter_max']['dim'])==0) and
              (SSD[0,0]['ShapeData']['diameter_max']['data'].size==1) ), \
        'Max. diameter data has wrong dimensions.'
      assert( (len(SSD[0,0]['ShapeData']['diameter_vol_eq']['dim'])==0) and
              (SSD[0,0]['ShapeData']['diameter_vol_eq']['data'].size==1) ), \
        'Vol. equ. diameter data has wrong dimensions.'
      assert( (len(SSD[0,0]['ShapeData']['diameter_area_eq_aerodynamical']['dim'])==0) and
              (SSD[0,0]['ShapeData']['diameter_area_eq_aerodynamical']['data'].size==1) ), \
        'Aerodyn. area equ. diameter data has wrong dimensions.'
      M.mass                = np.float(SSD[0,0]['ShapeData']['mass']['data'])
      M.diameter_max        = np.float(SSD[0,0]['ShapeData']['diameter_max']['data'])
      M.diameter_volume_equ = np.float(SSD[0,0]['ShapeData']['diameter_vol_eq']['data'])
      M.diameter_area_equ_aerodynamical \
                            = np.float(SSD[0,0]['ShapeData']['diameter_area_eq_aerodynamical']['data'])
  else:
    raise Exception( 'Non-typhon *ScatteringMetaData* not yet implemented.' )

  phase = SSD[0,0]['ShapeData']['phase']

  # Basic data of S
  if use_typhon:
    S = tas.SingleScatteringData()
    S.version     = 3
    if (SSD.size>0):
      try:
        ptypeID = db_ptype[SSD[0,0]['SingleScatteringData']['orient_type']]
      except KeyError:
        raise Exception( \
          "Database ptype '%s' is unknown." %SSD[0,0]['SingleScatteringData']['orient_type'] )
      try:
        S.ptype       = a_ptype[ptypeID]
      except KeyError:
        raise Exception( \
          "No ARTS ptype equivalent defined for database ptype '%s' (internal ID=%i)." \
          %(SSD[0,0]['SingleScatteringData']['orient_type'],ptypeID) )
      S.description = SSD[0,0]['ShapeData']['description']
      S.f_grid  = freq
      S.T_grid  = temp
  else:
    raise Exception( 'Non-typhon *SingleScatteringData* not yet implemented' )

  # Fill S according to ptype
  nf = freq.size
  nt = temp.size

  if (SSD.size>0):

    #####
    # totally random orientation
    #####
    if ( S.ptype==a_ptype[20] ):

      pha_inds = [ [1,2,0,0], [2,3,0,0], [0,0,4,-5], [0,0,5,6] ]

      # Empty aa_grid
      S.aa_grid = np.array([])

      # Create a common za_grid
      # (instead of using union1d over and over, we concatenate all za_scat data
      #  and derive the unique value from complete array. to be tested whether
      #  this is indeed faster than union1d use.)
      za_grid = np.array([])
      for f in np.arange(freq.size):
        for t in np.arange(temp.size):
          if (not nodata[f,t]):
            # A number of asserts to double-check the input
            assert( abs(M.mass-SSD[f,t]['ShapeData']['mass']['data'])<1e-12 ), \
              'Mass info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
              %(freq[f]*1e-9,temp[t])
            assert( abs(M.diameter_max-SSD[f,t]['ShapeData']['diameter_max']['data'])<1e-12 ), \
              'Dmax info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
              %(freq[f]*1e-9,temp[t])
            assert( abs(M.diameter_volume_equ-SSD[f,t]['ShapeData']['diameter_vol_eq']['data'])<1e-12 ), \
              'Dveq info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
              %(freq[f]*1e-9,temp[t])

            try:
              ptypeID = db_ptype[SSD[f,t]['SingleScatteringData']['orient_type']]
            except KeyError:
              raise Exception( \
                "Database ptype for SSD at f=%.1fGHz and T=%.1fK ('%s') is unknown." \
                %(freq[f]*1e-9,temp[t],SSD[f,t]['SingleScatteringData']['orient_type']) )
            try:
              ptypeName = a_ptype[ptypeID]
            except KeyError:
              raise Exception( \
                "No ARTS ptype equivalent defined for database ptype '%s' (internal ID=%i) for SSD at f=%.1fGHz and T=%.1fK." \
                %(SSD[f,t]['SingleScatteringData']['orient_type'],ptypeID,freq[f]*1e-9,temp[t]) )
            assert( S.ptype==ptypeName ), \
              'ptype for SSD at f=%.1fGHz and T=%.1fK inconsistent with initial value.' \
              %(freq[f]*1e-9,temp[t])

            if not( abs(freq[f]-SSD[f,t]['SingleScatteringData']['frequency']['data'])<1e3 ):
              raise Exception( \
                'Freq info in SSD at f=%.1fGHz and T=%.1fK inconsistent with freq grid value.' \
                %(freq[f]*1e-9,temp[t]) )
            if not( abs(temp[t]-SSD[f,t]['SingleScatteringData']['temperature']['data'])<0.001 ):
              raise Exception( \
                'Temp info in SSD at f=%.1fGHz and T=%.1fK inconsistent with temp grid value.' \
                %(freq[f]*1e-9,temp[t]) )
            assert( pha_inds==(SSD[f,t]['SingleScatteringData']['phaMat_index']['data']).tolist() ), \
              'phamat indexing info in SSD at f=%.1fGHz and T=%.1fK inconsistent ptype.' \
              %(freq[f]*1e-9,temp[t])

            assert( (SSD[f,t]['SingleScatteringData']['za_scat']['data'][0]==0.) and
                    (SSD[f,t]['SingleScatteringData']['za_scat']['data'][-1]==180.) ), \
              'za_scat at f=%.1fGHz and T=%.1fK does not fulfill basic requirement\n' + \
              'of first value being at 0deg and/or last value at 180deg.' %(freq[f]*1e-9,temp[t])
            za_grid = np.append(za_grid,SSD[f,t]['SingleScatteringData']['za_scat']['data'])
      S.za_grid = np.unique(za_grid)
      nza = S.za_grid.size
      S.pha_mat_data = np.empty((nf,nt,nza,1,1,1,6))*np.NAN
      S.ext_mat_data = np.empty((nf,nt,1,1,1))*np.NAN
      S.abs_vec_data = np.empty((nf,nt,1,1,1))*np.NAN
      g              = np.empty((nf,nt))*np.NAN
      for f in np.arange(freq.size):
        for t in np.arange(temp.size):
          if (not nodata[f,t]):
            # Fill data fields

            # Old, linear za-interpolation only version:
            #for i in np.arange(S.pha_mat_data.shape[-1]):
            #  S.pha_mat_data[f,t,:,0,0,0,i] = \
            #    np.interp( S.za_grid,
            #              SSD[f,t]['SingleScatteringData']['za_scat']['data'],
            #              SSD[f,t]['SingleScatteringData']['phaMat_data']['data'][i,0,0,0,:] )

            # selectable za-interpolation method:
            if (interpm=='pchip'):
              if not( SSD[f,t]['SingleScatteringData']['za_scat'].size>2 ):
                raise Exception( \
                  "Interpolation method '%s' requires at least 3 sample points" + \
                  " for evaluating the function,\n" + \
                  "but SSD at f=%.1fGHz and T=%.1fK provides only %i" + \
                  " zenith angle sample points." \
                  %(interpm,freq[f]*1e-9,temp[t],i,
                    SSD[f,t]['SingleScatteringData']['za_scat'].size) )
              fi = spi.PchipInterpolator(SSD[f,t]['SingleScatteringData']['za_scat']['data'],
                                          SSD[f,t]['SingleScatteringData']['phaMat_data']['data'],
                                          axis=4)
              S.pha_mat_data[f,t,...] = fi(S.za_grid).T
            else:
              fi = spi.interp1d(SSD[f,t]['SingleScatteringData']['za_scat']['data'],
                                SSD[f,t]['SingleScatteringData']['phaMat_data']['data'],
                                axis=4,kind=interpm,assume_sorted=True)
              S.pha_mat_data[f,t,...] = fi(S.za_grid).T

            S.ext_mat_data[f,t,...] = SSD[f,t]['SingleScatteringData']['extMat_data']['data']
            S.abs_vec_data[f,t,...] = SSD[f,t]['SingleScatteringData']['absVec_data']['data']
            g[f,t]                  = rttov.assp2g(S)
    #####
    # azimuthally random orientation
    #####
    elif ( S.ptype==a_ptype[30] ):
      #####
      # Separate according to phase: liquid (gridded data) vs ice (spherical harmonics)
      #####
      if phase=='ice':

          try:
            from . import sph
          except ImportError:
            try: #relative import does not work from scripts in the package folder
              import sph
              # from Python import sph
            except ImportError:
              raise Exception( 'Module *sph* of database interface required but not found.' )
          except Exception as e:
              print('Module *assp* requires module *sph*, but failed to import it:\n%s' %str(e))

          S.za_grid=SSD[0,0]['SingleScatteringData']['za_inc']['data']

          for f in np.arange(freq.size):
            for t in np.arange(temp.size):
              if (not nodata[f,t]):
                # A number of asserts to double-check the input
                assert( abs(M.mass-SSD[f,t]['ShapeData']['mass']['data'])<1e-12 ), \
                  'Mass info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
                  %(freq[f]*1e-9,temp[t])
                assert( abs(M.diameter_max-SSD[f,t]['ShapeData']['diameter_max']['data'])<1e-12 ), \
                  'Dmax info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
                  %(freq[f]*1e-9,temp[t])
                assert( abs(M.diameter_volume_equ-SSD[f,t]['ShapeData']['diameter_vol_eq']['data'])<1e-12 ), \
                  'Dveq info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
                  %(freq[f]*1e-9,temp[t])

                try:
                  ptypeID = db_ptype[SSD[f,t]['SingleScatteringData']['orient_type']]
                except KeyError:
                  raise Exception( \
                    "Database ptype for SSD at f=%.1fGHz and T=%.1fK ('%s') is unknown." \
                    %(freq[f]*1e-9,temp[t],SSD[f,t]['SingleScatteringData']['orient_type']) )
                try:
                  ptypeName = a_ptype[ptypeID]
                except KeyError:
                  raise Exception( \
                    "No ARTS ptype equivalent defined for database ptype '%s' (internal ID=%i) for SSD at f=%.1fGHz and T=%.1fK." \
                    %(SSD[f,t]['SingleScatteringData']['orient_type'],ptypeID,freq[f]*1e-9,temp[t]) )
                assert( S.ptype==ptypeName ), \
                  'ptype for SSD at f=%.1fGHz and T=%.1fK inconsistent with initial value.' \
                  %(freq[f]*1e-9,temp[t])

                if not( abs(freq[f]-SSD[f,t]['SingleScatteringData']['frequency']['data'])<1e3 ):
                  raise Exception( \
                    'Freq info in SSD at f=%.1fGHz and T=%.1fK inconsistent with freq grid value.' \
                    %(freq[f]*1e-9,temp[t]) )
                if not( abs(temp[t]-SSD[f,t]['SingleScatteringData']['temperature']['data'])<0.001 ):
                  raise Exception( \
                    'Temp info in SSD at f=%.1fGHz and T=%.1fK inconsistent with temp grid value.' \
                    %(freq[f]*1e-9,temp[t]) )
                assert( (SSD[f,t]['SingleScatteringData']['za_inc']['data'][0]==0.) and
                        (SSD[f,t]['SingleScatteringData']['za_inc']['data'][-1]==180.) ), \
                  'za_inc at f=%.1fGHz and T=%.1fK does not fulfill basic requirement\n' + \
                  'of first value being at 0deg and/or last value at 180deg.' %(freq[f]*1e-9,temp[t])

                # for now we require all data to be on identical za_inc grids
                assert( (S.za_grid==SSD[f,t]['SingleScatteringData']['za_inc']['data']).all() ), \
                  'za_inc grid f=%.1fGHz and T=%.1fK inconsistent with global za_grid.' \
                  %(freq[f]*1e-9,temp[t])

          pha_inds = SSD[0,0]['SingleScatteringData']['phaMat_index']['data']
          ext_inds = SSD[0,0]['SingleScatteringData']['extMat_index']['data']
          abs_inds = SSD[0,0]['SingleScatteringData']['absVec_index']['data']

          # Set scattering grid sizes
          #   ARTS requires za_inc==za_sca, sph requires (as implemented)
          #    aa_sca[0-180]=za_sca. Hence here we fix everything to the za_inc, which
          #    is the coarser grid in our ADDA calcs.
          #   Basically, that could be replaced by two (independent) input parameters,
          #    one for za, one for aa. However, that requires more interpolations
          #    (and rotations?) in sph.
          nza = len(SSD[0,0]['SingleScatteringData']['za_inc']['data'])
          g    = np.empty((nf,nt,len(S.za_grid)))*np.NAN
          g_i  = np.empty((len(S.za_grid)))*np.NAN
          got_grid=False
          for f in np.arange(freq.size):
            for t in np.arange(temp.size):
              if (not nodata[f,t]):

                # extract optical properties
                ext_mat = SSD[f,t]['SingleScatteringData']['extMat_data']['data']
                abs_vec = SSD[f,t]['SingleScatteringData']['absVec_data']['data']

                phamat_real = SSD[f,t]['SingleScatteringData']['phaMat_data_real']['data']
                phamat_imag = SSD[f,t]['SingleScatteringData']['phaMat_data_imag']['data']
                phamat_sph = phamat_real+1j*phamat_imag
               
                # get truncation level
                lmax = sph.get_lmax(SSD[f,t]['SingleScatteringData']['sph_coeffs_scat']['data'])
                # transform spherical harmonics to grid
                [phase_matrix_reg,theta_s,phi_s] = \
                  sph.create_regular_representation(phamat_sph, lmax,
                                                    grid_size_theta=nza, gridtype='regular')
                  
                ###############################################################                                    
                # Settings for asymetry parameter
                  
                ph_real     = phamat_real[0,0,:,:]
                ph_imag     = phamat_imag[0,0,:,:]
                ph_sph      = np.complex128(ph_real+1j*ph_imag)
                
                theta_inc = S.za_grid
                grid_size_theta=90
                
                for i in range(len(theta_inc)):
                  # get truncation level of spherical harmonics series
                  max_ntrunc = sph.get_lmax(ph_sph[i,:])
                  if max_ntrunc % 2 ==0:
                      min_grid_size=2*lmax+2
                  else:
                      min_grid_size=2*(lmax+1)  
                  if min_grid_size<grid_size_theta:
                      min_grid_size=grid_size_theta
                  elif grid_size_theta<min_grid_size:
                      grid_size_theta=min_grid_size
                  if min_grid_size % 2 == 1:
                      min_grid_size=min_grid_size+1
                  # set up spherical harmonics object
                  sph_object = sph.Spharmt(min_grid_size*2,min_grid_size,max_ntrunc,1,gridtype='regular')
                  # rotate spherical harmonics series, so that forward direction points toward northpole
                  # and backward direction towards southpole.
                  phase_sph_rot = sph_object.yrotate(ph_sph[i,:], -theta_inc[i]*np.pi/180)
                  g_i[i]        = np.real(phase_sph_rot[1]/phase_sph_rot[0]/np.sqrt(4*np.pi)*2)
                ###############################################################
                assert( (S.za_grid==theta_s).all() ), \
                  'Spherical harmonics output polar grid at f=%.1fGHz and T=%.1fK' + \
                  ' inconsistent with global za_grid.' %(freq[f]*1e-9,temp[t])

                # allocate and set angle grids
                if got_grid==False:
                  S.aa_grid=phi_s
                  S.ext_mat_data = np.empty((nf,nt,len(S.za_grid),1,3))*np.NAN
                  S.abs_vec_data = np.empty((nf,nt,len(S.za_grid),1,2))*np.NAN
                  S.pha_mat_data = np.empty((nf,nt,len(S.za_grid),len(S.aa_grid),
                                            len(S.za_grid),1,16))*np.NAN
                  got_grid=True

                # Fill data fields
                S.ext_mat_data[f,t,...] = np.transpose(ext_mat,(2,1,0))
                S.abs_vec_data[f,t,...] = np.transpose(abs_vec,(2,1,0))
                S.pha_mat_data[f,t,...] = np.transpose(phase_matrix_reg,(3,4,2,1,0))
                g[f,t,...]              = np.array(g_i) # VB
      elif phase=='liquid':
          S.aa_grid = np.array([])
          S.za_grid = np.array([])

          pha_inds = [ [1,5,9,13], [2,6,10,14], [3,7,11,15], [4,8,12,16] ]

          for f in np.arange(freq.size):
            for t in np.arange(temp.size):
              if (not nodata[f,t]):
                # A number of asserts to double-check the input
                assert( abs(M.mass-SSD[f,t]['ShapeData']['mass']['data'])<1e-12 ), \
                  'Mass info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
                  %(freq[f]*1e-9,temp[t])
                assert( abs(M.diameter_max-SSD[f,t]['ShapeData']['diameter_max']['data'])<1e-12 ), \
                  'Dmax info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
                  %(freq[f]*1e-9,temp[t])
                assert( abs(M.diameter_volume_equ-SSD[f,t]['ShapeData']['diameter_vol_eq']['data'])<1e-12 ), \
                  'Dveq info for SSD at f=%.1fGHz and T=%.1fK inconsistent with meta data value.' \
                  %(freq[f]*1e-9,temp[t])
                try:
                  ptypeID = db_ptype[SSD[f,t]['SingleScatteringData']['orient_type']]
                except KeyError:
                  raise Exception( \
                    "Database ptype for SSD at f=%.1fGHz and T=%.1fK ('%s') is unknown." \
                    %(freq[f]*1e-9,temp[t],SSD[f,t]['SingleScatteringData']['orient_type']) )
                try:
                  ptypeName = a_ptype[ptypeID]
                except KeyError:
                  raise Exception( \
                    "No ARTS ptype equivalent defined for database ptype '%s' (internal ID=%i) for SSD at f=%.1fGHz and T=%.1fK." \
                    %(SSD[f,t]['SingleScatteringData']['orient_type'],ptypeID,freq[f]*1e-9,temp[t]) )
                assert( S.ptype==ptypeName ), \
                  'ptype for SSD at f=%.1fGHz and T=%.1fK inconsistent with initial value.' \
                  %(freq[f]*1e-9,temp[t])

                if not( abs(freq[f]-SSD[f,t]['SingleScatteringData']['frequency']['data'])<1e3 ):
                  raise Exception( \
                    'Freq info in SSD at f=%.1fGHz and T=%.1fK inconsistent with freq grid value.' \
                    %(freq[f]*1e-9,temp[t]) )
                if not( abs(temp[t]-SSD[f,t]['SingleScatteringData']['temperature']['data'])<0.001 ):
                  raise Exception( \
                    'Temp info in SSD at f=%.1fGHz and T=%.1fK inconsistent with temp grid value.' \
                    %(freq[f]*1e-9,temp[t]) )
                assert( pha_inds==(SSD[f,t]['SingleScatteringData']['phaMat_index']['data']).tolist() ), \
                  'phamat indexing info in SSD at f=%.1fGHz and T=%.1fK inconsistent ptype.' \
                  %(freq[f]*1e-9,temp[t])

                assert( (SSD[f,t]['SingleScatteringData']['za_scat']['data'][0]==0.) and
                        (SSD[f,t]['SingleScatteringData']['za_scat']['data'][-1]==180.) ), \
                  'za_scat at f=%.1fGHz and T=%.1fK does not fulfill basic requirement\n' + \
                  'of first value being at 0deg and/or last value at 180deg.' %(freq[f]*1e-9,temp[t])
                S.za_grid = np.append(S.za_grid,SSD[f,t]['SingleScatteringData']['za_scat']['data'])
                S.aa_grid = np.append(S.aa_grid,SSD[f,t]['SingleScatteringData']['aa_scat']['data'])

          # Unique grid
          S.za_grid = np.unique(S.za_grid)
          S.aa_grid = np.unique(S.aa_grid)

          nza = S.za_grid.size
          naa = S.aa_grid.size
          npha= 16
          S.ext_mat_data = np.empty((nf,nt,nza,1,3))*np.NAN
          S.abs_vec_data = np.empty((nf,nt,nza,1,2))*np.NAN
          S.pha_mat_data = np.empty((nf,nt,nza,naa,nza,1,npha))*np.NAN
          g              = np.empty((nf,nt,len(S.za_grid)))*np.NAN
          for f in np.arange(freq.size):
            for t in np.arange(temp.size):
              if (not nodata[f,t]):

                exti = spi.interp1d(SSD[f,t]['SingleScatteringData']['za_scat']['data'],
                                    SSD[f,t]['SingleScatteringData']['extMat_data']['data'],
                                    axis=2,kind=interpm,assume_sorted=True)
                S.ext_mat_data[f,t,...] = exti(S.za_grid).T

                absi = spi.interp1d(SSD[f,t]['SingleScatteringData']['za_scat']['data'],
                                    SSD[f,t]['SingleScatteringData']['absVec_data']['data'],
                                    axis=2,kind=interpm,assume_sorted=True)
                S.abs_vec_data[f,t,...] = absi(S.za_grid).T

                zasq,aasq,zaiq = np.meshgrid(S.za_grid,S.aa_grid,S.za_grid, indexing='ij')

                for ind in np.arange(npha):
                    f_npha = SSD[f,t]['SingleScatteringData']['phaMat_data']['data'][ind,0,:,:,:]
                    fi = spi.RegularGridInterpolator((SSD[f,t]['SingleScatteringData']['za_scat']['data'],
                                                      SSD[f,t]['SingleScatteringData']['aa_scat']['data'],
                                                      SSD[f,t]['SingleScatteringData']['za_inc']['data']),f_npha, method = interpm)

                    # Change the order of the points to match the input shape of the interpolation function.
                    extended_grid = np.rollaxis(np.array([zasq,aasq,zaiq]),0,4)
                    interpolated  = fi(extended_grid)
                    S.pha_mat_data[f,t,:,:,:,0,ind] = interpolated

    #####
    # unknown particle type
    #####
    else:
      raise Exception( "Ptype '%s' not (yet?) implemented." %S.ptype )

  else: # if (SSD.size>0)
    print('No SSD data available.\n'
          'Returning default-filled SingleScatteringData and ScatteringMetaData'
          ' (grids partly empty, data fields all empty).')

  return S, M, g
