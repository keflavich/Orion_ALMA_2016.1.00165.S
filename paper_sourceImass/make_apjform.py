with open('sourceImass.tex','r') as fhr:
    with open('sourceImass_apjform.tex', 'w') as fhw:
        for line in fhr:
            if 'figures/' in line:
                fhw.write(line.replace("figures/",""))
            else:
                fhw.write(line)
