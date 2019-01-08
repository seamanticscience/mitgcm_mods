# mitgcm_mods

Modifications to the core MITgcm code written by Jonathan Lauderdale for various projects.

This code is for checkpoint67f of MITgcm and it supercedes the the previous versions for checkpoint63m and checkpoint66j, which are older versions of MITgcm:
<a href="https://zenodo.org/badge/latestdoi/87344333"><img src="https://zenodo.org/badge/87344333.svg" alt="DOI"></a>

Some code in the bgc folder adds preformed tracer functionality to MITgcm (GCHEM_OPTIONS.h, gchem_forcing_sep.F and gchem_preformed_tracers.F and increase PTRACERS_num in PTRACERS_SIZE.h)

One might use this code (partially or in its entirety) by:

1.) Downloading the correct model version:

`git clone https://github.com/MITgcm/MITgcm.git`

`cd MITgcm`

`git checkout -b my_model checkpoint67f`

OR download from http://mitgcm.org/download/other_checkpoints/MITgcm_c67f.tar.gz, extract, and `cd MITgcm_c67f` into the model directory.

2.) Cloning the mitgcm_mods repository.

`git clone https://github.com/seamanticscience/mitgcm_mods.git`

`cd mitgcm_mods`

`git checkout -b my_mods master`

3.) Running the genmake2 script like so (for the 3 degree global model with biogeochem and kpp run on a single processor):

 /path/to/MITgcm/tools/genmake2 -rootdir /path/to/MITgcm/ -mods '/path/to/MITgcm/mitgcm_mods/with_kpp /path/to/MITgcm/mitgcm_mods/bgc /path/to/MITgcm/mitgcm_mods/code /path/to/MITgcm/mitgcm_mods/sp_code'
 
 4.) make depend (The individual "mod" directories will be symlinked into the build directory in reverse order.)
 
 5.) make
 
 Please let me know if you find any errors, or have trouble compiling!
 
