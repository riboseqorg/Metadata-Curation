
# entrypoint.sh
#!/bin/bash
set -e

# Ensure we're in the right directory
cd /app

# Create necessary directories if they don't exist
mkdir -p /app/temp_files
mkdir -p /app/SraRunInfo
mkdir -p /app/output

# Ensure proper permissions
chmod -R 777 /app/temp_files /app/SraRunInfo /app/output

# Run the Python script first if it's needed for SraRunInfo
if [ -f "pull_bioprojects.py" ]; then
    python3 pull_bioprojects.py
fi

# Run the main R script with error handling
Rscript metadata_main_script.R
exit_code=$?

if [ $exit_code -ne 0 ]; then
    echo "Error: R script failed with exit code $exit_code"
    # You might want to add additional error handling here
fi

exit $exit_code