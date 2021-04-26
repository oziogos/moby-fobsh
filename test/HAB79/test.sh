#!/bin/bash

# create pot and psf files
/deploy/utils/gaff_batch.py

#
/deploy/utils/deploy_fobsh_test_suite.py

#
/deploy/utils/run_fobsh_tests.py

#
/deploy/utils/correlate.py
