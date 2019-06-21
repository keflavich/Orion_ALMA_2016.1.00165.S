from astropy import units as u
ms_basepath = '/export/home/rng90003/orion/2016.1.00165.S/measurement_sets/'
ms_basepath = '/lustre/aginsbur/orion/2016.1.00165.S/measurement_sets/'


mses = {'b3': ['member.uid___A001_X88e_X1d9_calibrated.ms'],
        'b6': ['uid___A001_X88e_X1d3_calibrated.ms',
               'uid___A002_Xb925ef_X4334_calibrated.ms'],
        'b7': ['band7.ms', 'band7_lb.ms'],
       }

imaging_parameters = {
    'b3': {'imsize': [3840, 3840],
           'cell': ['0.02arcsec'],
           'threshold': '1mJy',
           'niter': 1000,
          },
    'b6':{'imsize': [3600, 3600],
          'cell': ['0.008arcsec'],
          'threshold': '10mJy',
          'niter': 1000,
         },
    'b7':{'imsize': [3840, 4720],
          'cell': ['0.006arcsec'],
          'threshold': '15mJy',
          'niter': 1000,
         }
}


line_to_image_list = [
    {'name': 'SiOv=0_2-1', 'frequency': 86.84696*u.GHz,'band': 3, },
    {'band': 6, 'name': 'H2Ov2=1_5(5,0)-6(4,3)', 'frequency': 232.6867*u.GHz, 'velocity_range': (-30, 40)},
    {'name': 'SiOv=0_8-7', 'frequency': 347.330824*u.GHz,'band': 7, },
    {'name': 'SiOv=1_2-1', 'frequency': 86.24343*u.GHz,'band': 3, },
    {'name': 'SiS_12-11', 'frequency': 217.817644*u.GHz, 'band': 6, },
    {'name': '29SiOv=0_5-4', 'frequency': 214.3857577*u.GHz,'band': 6, },
    {'name': '29SiOv=0_2-1', 'frequency': 85.759199*u.GHz,'band': 3, },
    {'name': '29SiOv=2_2-1', 'frequency': 85.640452*u.GHz,'band': 3, },
    {'name': 'SiOv=0_5-4', 'frequency': 217.10498*u.GHz,'band': 6, },
    {'name': 'SiOv=1_5-4', 'frequency': 215.59595*u.GHz,'band': 6, },
    {'name': 'H30a', 'frequency': 231.900928*u.GHz,'band': 6, },
    {'name': 'SO2_2-1_1', 'frequency': 86.09395*u.GHz,'band': 3, },
    {'name': 'SO3_2-2_1', 'frequency': 99.299905*u.GHz,'band': 3, },
    {'name': 'SOv=1_3_2-2_1', 'frequency': 98.758192*u.GHz,'band': 3, },
    {'name': 'SO8_8-7_7', 'frequency': 344.310612*u.GHz,'band': 7, },
    {'name': 'SO9_8-8_7', 'frequency': 346.528481*u.GHz,'band': 7, },
    {'name': '29SiOv=0_8-7', 'frequency': 342.9808425*u.GHz,'band': 7, },
    {'name': '29SiOv=1_8-7', 'frequency': 340.6118623*u.GHz,'band': 7, },
    {'name': '29SiOv=2_8-7', 'frequency': 338.2451561*u.GHz,'band': 7, },
    {'name': 'SiOv=1_8-7', 'frequency': 344.916247*u.GHz,'band': 7, },
    {'name': 'SiOv=2_8-7', 'frequency': 342.504607*u.GHz,'band': 7, },
    {'name': '30SiOv=0_8-7', 'frequency': 338.9300437*u.GHz,'band': 7, },
    {'name': '30SiOv=1_8-7', 'frequency': 336.6029763*u.GHz,'band': 7, },
    {'name': '30SiOv=2_8-7', 'frequency': 334.2781299*u.GHz,'band': 7, },
    {'name': 'SiS19-18', 'frequency': 344.779481*u.GHz,'band': 7, },
    {'name': 'H13CN1-0', 'frequency': 86.3401764*u.GHz,'band': 3, },
    {'name': 'HC15N1-0', 'frequency': 86.0549664*u.GHz,'band': 3, },
    {'name': 'HCN1-0', 'frequency': 88.6318473*u.GHz,'band': 3, },
    {'name': 'HC3N11-10', 'frequency': 100.076386*u.GHz,'band': 3, },
    {'name': 'HC3Nv7=1_11-10', 'frequency': 100.3224109*u.GHz,'band': 3, },
    {'name': 'CS2-1', 'frequency': 97.980953*u.GHz,'band': 3, },
    {'name': 'SO4_5-4_4', 'frequency': 100.029565*u.GHz,'band': 3, },
    {'name': 'SOv=1_4_5-4_4', 'frequency': 98.317606*u.GHz,'band': 3, },
    ]
for entry in line_to_image_list:
    if 'velocity_range' not in entry:
        entry['velocity_range'] = (-30, 40)
