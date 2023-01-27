#!/bin/bash

# Start from XToYggHbb_looper/.
rm -rf tmp_test_package/
mkdir -p tmp_test_package/
cd tmp_test_package/

tar xf $1 # Use absolute directory
cd XToYggHbb_looper/

# Try running a command
cd cpp/
bash runOutput_XToYggHbb.sh temp_data_test 2018 0 1 0 HHbbgg
