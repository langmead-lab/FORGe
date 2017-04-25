#!/bin/bash

dr=$(dirname "$0")
[ -z "$1" ] && echo "Specify destination directory as argument"
# can be "wet" (submit immediately) or
# "dry" (just print the qsub commands but don't run them)
type="$2"
stop_after="$3"
[ -z "$stop_after" ] && stop_after=2000

j=1
for i in `find . -name '*.sam' | grep -v 'fastq\.sam'` ; do
    dn=$(dirname "$i")
    bn=$(basename "$i")
    bbn=`echo ${bn} | sed 's/\.sam$//'`
    dest="$1/${dn}"
    fulldr=`pwd`
    if [ ! -f "${dest}/${bbn}.csv" ] ; then
        mkdir -p ${dest}
cat <<EOF >${dest}/.${bbn}.sh
#!/bin/sh
#PBS -l walltime=02:00:00
#PBS -e ${dest}/.${bn}.sh.e
#PBS -o ${dest}/.${bn}.sh.o

cd ${fulldr}/${dn}
sh ${dr}/annotate_sam.sh ${bn} ${dest}
EOF
        cmd="qsub ${dest}/.${bbn}.sh"
        if [ "$2" = "wet" ] ; then
            eval ${cmd}
            sleep 1
        else
            echo ${cmd}
        fi
        j=$(expr $j + 1)
        if [ $j -gt $stop_after ] ; then
            break
        fi
    else
        echo "Skipping \"$i\" because destination csv exists"
    fi
done
