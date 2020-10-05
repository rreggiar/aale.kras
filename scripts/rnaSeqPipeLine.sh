#!/bin/bash
#!/bin/bash

# rreggiar@ucsc.edu


scriptName=$(basename $0)
if [ $# -lt 3 ]; then
    echo "error: usage $scriptName trim=T salmon=T star=T inputDir salmonIndex starGenome adapterChoice "
    echo "example $scriptName "
    exit 1
fi

for module in trimmomaticRun.sh salmonRun.sh starRun.sh dateStamp.sh; do
	if [[ -f "$module" ]]; then
		echo "$module" "is present"
	else
		echo "$module" "is absent"
		exit 1
	fi
done

dateStamp="$(bash dateStamp.sh)"

trimArg="$1"
salmonArg="$2"
starArg="$3"
inputDir="$4"

if [[ "$trimArg"=='trim=T' ]]; then
	if [[ -f "$5" ]]; then
		adapterChoice="$5"
	else
		echo 'you must provide an adapter choice to trim [$5]'
		exit 1
	fi
fi

if [[ "$salmonArg"=='salmon=T' ]]; then
	if [[ -f "$6" ]]; then
		salmonIndex="$6"
	else
		echo 'you must provide a salmon index [$6]'
		exit 1
	fi
fi

if [[ "$starArg"=='star=T' ]]; then
	if [[ -f "$7" ]]; then
		starGenome="$7"
	else
		echo 'you must provide a STAR indexed genome [$7]'
		exit 1
	fi
fi

echo "script: $scriptName"
echo "time: $dateStamp"
echo "adapterChoice: $adapterChoice"
echo "salmonIndex: $salmonIndex"
echo "starGenome: $starGenome"


