#!/bin/bash

# Set current working dir to the directory of this script
cd "$(dirname "$0")"

for f in 'IFC_QPE' 'IFC_QPF' 'IFC_WHATIF_NORAIN' 'IFC_WHATIF_2IN24H'
do

echo "Cleaning " $f.

# Clean empty std out/err files
find out/$f.* -maxdepth 1 -size 0 | xargs -r -d '\n' rm

# Clean non empty std out/err files
ls -1tr out/$f.o* | head -n -5 | xargs -r -d '\n' rm
ls -1tr out/$f.e* | head -n -5 | xargs -r -d '\n' rm

done

# Clean .str files
ls -1tr out/forcing_rain_qpe_*.str | head -n -5 | xargs -r -d '\n' rm
ls -1tr out/forcing_rain_qpf_*.str | head -n -5 | xargs -r -d '\n' rm

# Clean .rec files
ls -1tr out/state_*.h5 | head -n -5 | xargs -r -d '\n' rm
ls -1tr out/forecast_norain_*.h5 | head -n -5 | xargs -r -d '\n' rm
ls -1tr out/forecast_2in24h_*.h5 | head -n -5 | xargs -r -d '\n' rm
ls -1tr out/forecast_qpf_*.h5 | head -n -5 | xargs -r -d '\n' rm
