Readme for the salt analysis
----------------------------

These files need to be run in this order:

 * line_stack_overlay_salts.py
 * find_salts_in_band.py
 * line_stack_analysis_3bands.py
 * kcl_rotation_diagram.py
 * salt_detection_stats.py

These are support files:

 * lines.py
 * salt_tables.py

The 'lines.py' table is a file created by hand containing line rest
frequencies.  Since it's done by hand, it isn't fully self-consistent and
probably contains some errors.  However, it's _only_ used for labeling figures
in ``line_stack_overlay_salts.py`` and for obtaining latex-formatted names of
species.  The detection and fitting analysis is done automatically.
