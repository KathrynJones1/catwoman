.. :changelog:
1.1.0 (26-01-26)
~~~~~~~~~~~~~~~~~~~
-[NE] upgraded catwoman to account for Mans Holmberg's upgrades and bug fixes (see PR #19).
-[NE] changed version to prep for major upgrade.

1.0.14 (19-12-22)
~~~~~~~~~~~~~~~~~~~
-updated version number to match pypi distributions.

1.0.13 (28-10-22)
~~~~~~~~~~~~~~~~~~~
-fixed installation bug (issue #15).

1.0.12 (24-08-20)
~~~~~~~~~~~~~~~~~~~
-fixed a bug with the phi angle being corrected for orbital motion
-fixed a bug when, due to orbital motion, phi increases above 90 degrees or below -90 degrees

1.0.10 + 1.0.11 (12-08-20)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
-made sure numpy>=1.16.2 is installed when installing catwoman

1.0.9 (10-08-20)
~~~~~~~~~~~~~~~~~
-actually removed unused files (see below)

1.0.8 (10-08-20)
~~~~~~~~~~~~~~~~~
-removed unused files from directory so they aren't include in installation

1.0.7 (29-06-20)
~~~~~~~~~~~~~~~~~
-fixed a bug with installing catwoman straight from source archive, involved removing unused python files from import

1.0.6 (27-06-20)
~~~~~~~~~~~~~~~~~~
- fixed setup.py file so catwoman specifies and downloads dependencies automatically when pip installed in a new environment

