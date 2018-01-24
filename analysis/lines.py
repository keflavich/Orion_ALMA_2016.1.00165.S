from astropy import units as u
disk_lines = {'Si34S_13-12':229.5008677*u.GHz,
              #'CH3OCHO_19(3,17)-18(2,16)E':229.590456*u.GHz,
              'Unknown_1':230.321535*u.GHz,
              'Unknown_2':230.780241*u.GHz,
              'Unknown_3':229.2474253*u.GHz, # maybe acetone?
              'U229.682': 229.682*u.GHz,
              'U229.819': 229.819*u.GHz,
              'U230.728': 230.728*u.GHz,
              'U230.970': 230.970*u.GHz,
              'H2Ov2=1_5(5,0)-6(4,3)': 232.6867*u.GHz,
              'Unknown_4': 232.5105*u.GHz,
              'Unknown_5': 217.9795*u.GHz,
              'SiS_12-11': 217.817644*u.GHz,
              'Unknown_6': 217.6660212*u.GHz,
              'Unknown_7': 217.5473*u.GHz, # deuterated formic acid?
              'HC3N_24-23': 218.324788*u.GHz,
              'Unknown_8': 214.7417*u.GHz, # Si33S?
              'Unknown_9': 214.9396309*u.GHz,
              'Unknown_10': 215.0092408*u.GHz, # CNHCO?  KCl?
              '29SiOv=0_5-4': 214.3857577*u.GHz,
              'Unknown_11': 214.54879*u.GHz,
              'Unknown_12': 214.637167*u.GHz, # SiSv=3 12-11?
              #'Unknown_13': 215.88712*u.GHz,
              '13CH3OH_4(2,2)-3(1,2)': 215.886979*u.GHz,
              'Unknown_14': 233.170797*u.GHz,
              'Unknown_15': 233.6083343*u.GHz,
              'Unknown_16': 232.16300227817646*u.GHz,
              '29SiOv=0_2-1': 85.759199*u.GHz,
              '29SiOv=2_2-1': 85.640452*u.GHz,
              'SiOv=1_2-1': 86.24343*u.GHz,
              'SiOv=0_2-1': 86.84696*u.GHz,
              'Unknown_B3_1': 87.2662137*u.GHz,
              'Unknown_B3_2': 87.1707681*u.GHz,
              'Unknown_B3_3': 89.2210022*u.GHz,
              'Unknown_B3_4': 89.1509141*u.GHz,
              'Unknown_B3_5': 88.5912644*u.GHz,
              'SiOv=0_5-4': 217.10498*u.GHz,
              'SiOv=1_5-4': 215.59595*u.GHz,
              'H30a': 231.900928*u.GHz, # not actually detected
              'U346.3': 346.3708*u.GHz,
              'U344.4': 344.4537*u.GHz,
              'U335.0': 335.0523*u.GHz,
              'U335.5': 335.5057*u.GHz,
              'U335.1': 335.1307*u.GHz,
              'U333.0': 333.0126*u.GHz,
              '29SiOv=0_8-7': 342.9808425*u.GHz,
              '29SiOv=1_8-7': 340.6118623*u.GHz,
              '29SiOv=2_8-7': 338.2451561*u.GHz,
              'SiOv=0_8-7': 347.330824*u.GHz,
              'SiOv=1_8-7': 344.916247*u.GHz,
              'SiOv=2_8-7': 342.504607*u.GHz,
              '30SiOv=0_8-7': 338.9300437*u.GHz,
              '30SiOv=1_8-7': 336.6029763*u.GHz,
              '30SiOv=2_8-7': 334.2781299*u.GHz,
             }

outflow_lines = {'Si34S_13-12': 229.499086*u.GHz,
                 'H2Ov2=1_5(5,0)-6(4,3)': 232.6867*u.GHz,
                 'SiOv=1_5-4': 215.59595*u.GHz,
                 'SiOv=0_5-4': 217.10498*u.GHz,
                 '29SiOv=0_5-4': 214.3857577*u.GHz,
                 '29SiOv=0_8-7': 342.9808425*u.GHz,
                 '29SiOv=1_8-7': 340.6118623*u.GHz,
                 '29SiOv=2_8-7': 338.2451561*u.GHz,
                 'SiOv=0_8-7': 347.330824*u.GHz,
                 'SiOv=1_8-7': 344.916247*u.GHz,
                 'SiOv=2_8-7': 342.504607*u.GHz,
                 '30SiOv=0_8-7': 338.9300437*u.GHz,
                 '30SiOv=1_8-7': 336.6029763*u.GHz,
                 '30SiOv=2_8-7': 334.2781299*u.GHz,
                }

absorbers = {'12CO2-1': 230.538*u.GHz,
            }

texnames = {'Si34S_13-12': 'Si$^{34}$S 13-12',
            #'CH3OCHO_19(3,17)-18(2,16)E':229.590456*u.GHz,
            'Unknown_1':'U230.322', #230.321535*u.GHz,
            'Unknown_2':'U230.780', #230.780241*u.GHz,
            'Unknown_3':'U229.247', #229.2474253*u.GHz, # maybe acetone?
            'U229.682': 'U229.682',
            'U229.819': 'U229.819',
            'U230.728': 'U230.728',
            'U230.970': 'U230.970',
            'H2Ov2=1_5(5,0)-6(4,3)': 'H$_2$O v$_2$=1 $5_{5,0}-6_{4,3}$',
            'Unknown_4': 'U232.511', #232.5105*u.GHz,
            'Unknown_5': 'U217.780', #217.9795*u.GHz,
            'SiS_12-11': 'SiS 12-11',
            'Unknown_6': 'U217.666', #217.6660212*u.GHz,
            'Unknown_7': 'U217.547', #217.5473*u.GHz, # deuterated formic acid?
            'HC3N_24-23': 'HC$_3$N 24-23',
            'Unknown_8': 'U214.742', #214.7417*u.GHz, # Si33S?
            'Unknown_9': 'U214.940', #214.9396309*u.GHz,
            'Unknown_10': 'U215.009', #215.0092408*u.GHz, # CNHCO?  KCl?
            '29SiOv=0_5-4': '$^{29}$SiO v=0 J=5-4', #214.3857577*u.GHz,
            'Unknown_11': 'U214.549', #214.54879*u.GHz,
            'Unknown_12': 'U214.637', #214.637167*u.GHz, # SiSv=3 12-11?
            #'Unknown_13': 215.88712*u.GHz,
            '13CH3OH_4(2,2)-3(1,2)': '$^{13}$CH$_3$OH $4_{2,2}-3_{1,2}$',
            'Unknown_14': 'U233.171', #233.170797*u.GHz,
            'Unknown_15': 'U233.608', #233.6083343*u.GHz,
            'Unknown_16': 'U232.163', #232.16300227817646*u.GHz,
            '29SiOv=0_2-1': '$^{29}$SiO v=0 J=2-1',
            '29SiOv=2_2-1': '$^{29}$SiO v=2 J=2-1',
            'SiOv=1_2-1': 'SiO v=1 J=2-1',
            'SiOv=0_2-1': 'SiO v=0 J=2-1',
            'Unknown_B3_1': 'U87.266', #87.2662137*u.GHz,
            'Unknown_B3_2': 'U87.171', #87.1707681*u.GHz,
            'Unknown_B3_3': 'U89.221', #89.2210022*u.GHz,
            'Unknown_B3_4': 'U89.151', #89.1509141*u.GHz,
            'Unknown_B3_5': 'U88.591', #88.5912644*u.GHz,
            'SiOv=0_5-4': 'SiO v=0 J=5-4',
            'SiOv=1_5-4': 'SiO v=1 J=5-4',
            'H30a': 'H30$\\alpha$', #231.900928*u.GHz, # not actually detected
            'U346.3': 'U346.371', #346.3708*u.GHz,
            'U344.4': 'U344.454', #344.4537*u.GHz,
            'U335.0': 'U335.052', #335.0523*u.GHz,
            'U335.5': 'U335.506', #335.5057*u.GHz,
            'U335.1': 'U335.131', #335.1307*u.GHz,
            'U333.0': 'U333.013', #333.0126*u.GHz,
            '29SiOv=0_8-7': '$^{29}$SiO v=0 J=8-7',
            '29SiOv=1_8-7': '$^{29}$SiO v=1 J=8-7',
            '29SiOv=2_8-7': '$^{29}$SiO v=2 J=8-7',
            'SiOv=0_8-7': 'SiO v=0 J=8-7',
            'SiOv=1_8-7': 'SiO v=1 J=8-7',
            'SiOv=2_8-7': 'SiO v=2 J=8-7',
            '30SiOv=0_8-7': '$^{30}$SiO v=0 J=8-7',
            '30SiOv=1_8-7': '$^{30}$SiO v=1 J=8-7',
            '30SiOv=2_8-7': '$^{30}$SiO v=2 J=8-7',
            '12CO2-1': '$^{12}$CO 2-1',
            }
