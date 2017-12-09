from astropy import units as u
disk_lines = {'Si34S_13-12':229.5008677*u.GHz,
              'CH3OCHO_19(3,17)-18(2,16)E':229.590456*u.GHz,
              'Unknown_1':230.321535*u.GHz,
              'Unknown_2':230.780241*u.GHz,
              'Unknown_3':229.2474253*u.GHz, # maybe acetone?
              'H2Ov2=1_5(5,0)-6(4,3)': 232.6867*u.GHz,
              'Unknown_4': 232.511*u.GHz,
              'SiS_12-11': 217.817644*u.GHz,
              'Unknown_5': 217.9802311*u.GHz,
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
             }

outflow_lines = {'Si34S_13-12': 229.499086*u.GHz,
                 'H2Ov2=1_5(5,0)-6(4,3)': 232.6867*u.GHz,
                 'SiOv=1_5-4': 215.59595*u.GHz,
                 'SiOv=0_5-4': 217.10498*u.GHz,
                 '29SiOv=0_5-4': 214.3857577*u.GHz,
                }
