# mitgcm_mods
Modifications to the core MITgcm code that I have written

This code is currently for an older version of MITgcm (63m, http://mitgcm.org/download/other_checkpoints/MITgcm_c63m.tar.gz) 

Some code in the bgc folder adds preformed tracer functionality to MITgcm (GCHEM_OPTIONS.h, gchem_forcing_sep.F and gchem_preformed_tracers.F and increase PTRACERS_num in PTRACERS_SIZE.h)

One might use this code (partially or in its entirety) by:

1.) Downloading the correct model version.

2.) Cloning the mitgcm_mods repository.

3.) Running the genmake2 script like so (for the 3 degree global model with biogeochem and kpp run on a single processor):

 /path/to/MITgcm/tools/genmake2 -rootdir /path/to/MITgcm/ -mods '/path/to/MITgcm/mitgcm_mods/with_kpp /path/to/MITgcm/mitgcm_mods/bgc /path/to/MITgcm/mitgcm_mods/code /path/to/MITgcm/mitgcm_mods/sp_code'
 
 4.) make depend (The individual "mod" directories will be symlinked into the build directory in reverse order.)
 
 5.) make
 
 Please let me know if you find any errors, or have trouble compiling!
