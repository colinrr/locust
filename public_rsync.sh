#!/bin/bash

# Basic safe-ish rsync for script in public repo.
# Anything not a .m file needs extra treatment
# Add -n option for dry run

rsync -avrh *.m ../locust-public/