*******************************************
Tools for the package setup and development
*******************************************

All scripts, that are intended to be use in the development/deployment process
are located here. Scripts that are used to produce images for the documentation are
located in the "images" sub-folder. The scripts within the "images" sub-folder usually
have the same name as the .png file that they produce and therefore are not documented
in this README.

Note: Makefile is not located here, but it heavily relies on the scripts in this
directory.

All scripts are intended to be run within the root directory of the project.

Scripts
=======

*   check-release-metadata.py
    Performs several check of the source code metadata before each release.
