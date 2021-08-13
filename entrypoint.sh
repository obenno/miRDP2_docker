#!/bin/bash --login
# The --login ensures the bash configuration is loaded,
# enabling Conda.
set -euo pipefail
source $(conda info --base)/etc/profile.d/conda.sh
conda activate miRDP2
exec /app/miRDP2-v1.1.4_pipeline.bash "$@"
