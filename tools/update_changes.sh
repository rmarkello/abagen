#!/bin/bash
#
# Collects the pull-requests since the latest release and
# aranges them in the CHANGELOG txt file.
#
# This is a script to be run before releasing a new version.
#
# Usage:
#
#    $ /bin/bash update_changes.sh v0.0.1
#
# Originally authored by @oesteban (github.com/oesteban) for fmriprep, licensed
# under BSD-3 (github.com/poldracklab/fmriprep/blob/master/LICENSE)

# Setting      # $ help set
set -u         # Treat unset variables as an error when substituting.
set -x         # Print command traces before executing command.

if [ $# -ne 1 ]; then
    echo "Need to provide a version tag"
    exit 1
fi

# Allow running this from whatever directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $(dirname "${DIR}")

# Check whether the Upcoming release header is present
head -1 CHANGELOG | grep -q Upcoming
UPCOMING=$?
if [[ "$UPCOMING" == "0" ]]; then
    head -n3 CHANGELOG >> newchanges
fi

# Elaborate today's release header
HEADER="$1 ($(date '+%B %d, %Y'))"
echo $HEADER >> newchanges
echo $( printf "%${#HEADER}s" | tr " " "=" ) >> newchanges
echo "" >> newchanges

# Search for PRs since previous release
git log --grep="Merge pull request" `git describe --tags --abbrev=0`..HEAD --pretty='format:  * %b %s' | sed 's/Merge pull request \#\([^\d]*\)\ from \([^/]*\)\/.*/(\#\1), @\2/' >> newchanges
echo "" >> newchanges
echo "" >> newchanges

# Add back the Upcoming header if it was present
if [[ "$UPCOMING" == "0" ]]; then
    tail -n+4 CHANGELOG >> newchanges
else
    cat CHANGELOG >> newchanges
fi

# Replace old CHANGELOG with new file
mv newchanges CHANGELOG
