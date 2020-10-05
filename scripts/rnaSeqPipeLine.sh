#!/bin/bash
#!/bin/bash

# rreggiar@ucsc.edu


scriptName=$(basename $0)
if [ $# -lt 6 ]; then
    echo "error: usage $scriptName trim=T salmon=T star=T inputDir salmonIndex starGenome adapterChoice "
    echo "example $scriptName "
    exit 1
fi

for module in trimmomaticRun.sh salmonRun.sh starRun.sh; do
	if [[ -f "$module" ]]; then
		echo "$module" "is present"
	else
		echo "$module" "is absent"
		exit 1
	fi
done

dateStamp="$(bash dateStamp.sh)"
echo "script: $scriptName"
echo "time: $dateStamp"
