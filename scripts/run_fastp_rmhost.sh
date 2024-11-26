#!/usr/bin/bash
shopt -s expand_aliases

if [ $# -ne 3 ];then
    echo -e "Usage: $0 [fq|fq1,fq2] [host] [out_prefix]\e[0m"
    echo -e " \033[31mOptional host:\033[0m:"
    printf "%16s:    %-20s %-20s\n" human Homo_sapiens GCF_000001405.40_GRCh38
    printf "%16s:    %-20s %-20s\n" pig Sus_scrofa GCA_000003025.6
    printf "%16s:    %-20s %-20s\n" chicken Gallus_gallus GCF_016699485.2
    exit 2
fi

fq=$1; host=$2; out=$3; trds=12
alias fastp="/share/data1/software/binary/fastp"
alias bowtie2="/share/data1/software/bowtie2-2.5.4/bowtie2"

#index
declare -A host_index
host_index["human"]="/share/data1/database/genome_host/Homo_sapiens_human/Homo_sapiens_human"
host_index["pig"]="/share/data1/database/genome_host/Sus_scrofa_pig/Sus_scrofa_pig"
host_index["chicken"]="/share/data1/database/genome_host/Gallus_gallus_chicken/Gallus_gallus_chicken"
index=${host_index["${host}"]}

( [ -f ${out}_map2host.log ] && grep -q 'overall' ${out}_map2host.log ) &&\
    echo -e "[$(date +%Y-%m-%d\ %H:%M:%S)] Skip sample: $out." && exit 0

if [[ ${fq} =~ "," ]];then
    fq1=$(echo $fq | cut -d "," -f1)
    fq2=$(echo $fq | cut -d "," -f2)
    fastp -w $trds -q 20 -u 30 -n 5 -y -Y 30 -l 80 --trim_poly_g \
        -i $fq1 -I $fq2 -o ${out}_clean_1.fq.gz -O ${out}_clean_2.fq.gz -h /dev/null -j /dev/null 1>${out}_fastp.log 2>&1
    bowtie2 --end-to-end --mm --fast -p $trds -x $index --no-head -1 ${out}_clean_1.fq.gz -2 ${out}_clean_2.fq.gz 2> ${out}_map2host.log |\
        perl -ne 'chomp;@s=split /\s+/;if($s[1]==77){print "\@$s[0]/1\n$s[9]\n+\n$s[10]\n";}elsif($s[1]==141){print STDERR "\@$s[0]/2\n$s[9]\n+\n$s[10]\n";}' \
        > >(pigz> ${out}_rmhost_1.fq.gz) 2> >(pigz > ${out}_rmhost_2.fq.gz)
    chmod 444 ${out}_fastp.log ${out}_map2host.log ${out}_rmhost_1.fq.gz ${out}_rmhost_2.fq.gz
    rm ${out}_clean_1.fq.gz ${out}_clean_2.fq.gz
else
    fastp -w $trds -q 20 -u 30 -n 5 -y -Y 30 -l 80 --trim_poly_g \
        -i $fq -o ${out}_clean.fq.gz -h /dev/null -j /dev/null 1>${out}_fastp.log 2>&1
    bowtie2 --end-to-end --mm --fast -p $trds -x $index --no-head -U ${out}_clean.fq.gz 2> ${out}_map2host.log |\
        perl -ne 'chomp;@s=split /\s+/;if($s[1]==4){print "\@$s[0]\n$s[9]\n+\n$s[10]\n";}' \
        > >(pigz > ${out}_rmhost.fq.gz)
    chmod 444 ${out}_fastp.log ${out}_map2host.log ${out}_rmhost.fq.gz
    rm ${out}_clean.fq.gz
fi