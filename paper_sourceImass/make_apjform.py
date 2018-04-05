import os
import tarfile
import shutil

if os.path.exists('apj_tar'):
    shutil.rmtree('apj_tar')
os.mkdir('apj_tar')

with open('sourceImass.tex','r') as fhr:
    with open('sourceImass_apjform.tex', 'w') as fhw:
        for line in fhr:
            if 'figures/' in line:
                fhw.write(line.replace("figures/",""))
                fn = line.split("figures/")[-1].replace("{","").replace("}","").strip()
                print(fn)
                if not os.path.exists('apj_tar/{0}'.format(fn)):
                    os.link('figures/{0}'.format(fn), 'apj_tar/{0}'.format(fn))
            else:
                fhw.write(line)


critical_files = [
    'apjmacros.tex',
    'authors.tex',
    'continuum_beams.tex',
    'continuum_fit_parameters.tex',
    'cube_metadata.tex',
    'gitstuff.tex',
    'image_metadata.tex',
    'macros.tex',
    'obs_metadata.tex',
    'preface.tex',
    'solobib.tex',
    'sourceImass_apjform.tex',
    'unknown_line_freqs.tex',
    'extracted.bib',
    'sourceImass.bbl',
]

for fn in critical_files:
    if not os.path.exists('apj_tar/{0}'.format(fn)):
        os.link(fn, 'apj_tar/{0}'.format(fn))

with tarfile.open('apj.tgz', mode='w:gz') as tf:
    tf.add('apj_tar')
