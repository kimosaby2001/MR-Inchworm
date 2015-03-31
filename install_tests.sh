#!/bin/bash

echo "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""
echo 'Performing Unit Tests of Build'
echo ' '
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"

if [ -e "MR_Inchworm/mr_inchworm" ]
then
	echo "MR-Inchworm:      has been Installed Properly"
else
	echo "MR-Inchworm Installation appears to have FAILED"
fi
if [ -e "Fasta_Splitter/Fasta_Splitter" ]
then
	echo "Fasta Splitter:   has been Installed Properly"
else
	echo "C Fasta Splitter Installation appears to have FAILED"
fi

if [ -e "fastool/fastool" ]
then
	echo "fastool:          has been Installed Properly"
else
	echo "fastool Installation appears to have FAILED"
fi

if [ -e "Collectl/bin/collectl" ]
then
        echo "collectl:         has been Installed Properly"
else
        echo "collectl Installation appears to have FAILED"
fi

