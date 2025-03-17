#!/bin/bash

if [ -d "digested_peptide_library" ]; then
    rm -r digested_peptide_library
else
    echo "digested_peptide_library: No such file or directory. It has already been deleted."
fi

if [ -d "digested_glycopeptide_library" ]; then
    rm -r digested_glycopeptide_library
else
    echo "digested_glycopeptide_library: No such file or directory. It has already been deleted."
fi

if [ -d "glycosite_library" ]; then
    rm -r glycosite_library
else
    echo "glycosite_library: No such file or directory. It has already been deleted."
fi

if [ -d "logs" ]; then
    rm -r logs
else
    echo "logs: No such file or directory. It has already been deleted."
fi

echo "Digest libraries have been deleted. Have a Fantastic Day!"
