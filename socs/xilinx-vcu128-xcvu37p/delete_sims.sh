#!/bin/bash

# Absolute path to the directory
SIM_DIR="/home/manas/NEW/esp_manas/accelerators/catapult_hls/kalman_filter_sysc_catapult/hw/sim"
# SIM_DIR="/home/manas/NEW/esp_manas/accelerators/catapult_hls/mac_sysc_catapult/hw/sim"

# List of files to delete
FILES=(
    "accelerator_output.txt"
    "kalman_filter.o"
    "log.txt"
    "sc_main.d"
    "testbench.d"
    "kalman_filter.d"
    "kalman_filter_sysc_catapult"
    "sc_main.o"
    "testbench.o"
)

# FILES=(
#     "accelerator_output.txt"
#     "mac.o"
#     "log.txt"
#     "sc_main.d"
#     "testbench.d"
#     "mac.d"
#     "mac_sysc_catapult"
#     "sc_main.o"
#     "testbench.o"
# )

# Check if the directory exists
if [ -d "$SIM_DIR" ]; then
    echo "Deleting files from $SIM_DIR..."
    for FILE in "${FILES[@]}"; do
        FILE_PATH="$SIM_DIR/$FILE"
        if [ -e "$FILE_PATH" ]; then
            rm -f "$FILE_PATH"
            echo "Deleted: $FILE_PATH"
        else
            echo "File not found: $FILE_PATH"
        fi
    done
    echo "File deletion completed."
else
    echo "Directory $SIM_DIR does not exist."
    exit 1
fi
