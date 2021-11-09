#!/bin/sh
export CONFFILE=$1
source ${CONFFILE}

RESULT=0
for vcf in ${DBDATADIR}/chr*.vcf.gz
do
   if [ ! -e ${vcf} ]; then
      # bash returns the regex when no files are found so need to
      # check that the returned file exists
      continue
   fi
   if [[ "$RESULT" == 1 ]]; then
      # break out of for loop once tabix has been run
      break
   fi
   # if all's good, load samples by parsing tabix output
   tabix -h ${vcf} 0:0-0 | python ${PYLDIR}/load_samples.py --assaytype=${ASSAYTYPE}
   RESULT=1
done

if [[ "$RESULT" == 0 ]]; then
  echo "Error - no VCF files found in ${DBDATADIR}."
fi
