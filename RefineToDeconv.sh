#! /bin/bash


echo '===================================================================='


for dir in */
do
	cd $dir
	for file in *
	do
		if [[ $file == *.sam ]]
			then
			samp=$(echo $file | cut -d "." -f 1 )
			# echo $samp
			if [[ ! -f ${samp}_covar_deconv.tsv ]] && [[ ! -f ${samp}_covars.tsv ]]
				then
				python /mnt/g/MU_WW/SAM_Refiner/SAM_Refiner.py -r /mnt/g/MU_WW/SARS2/GP.fasta --collect 0 --indel 0 --nt_call 0 --read 0 -S $file --chim_rm 0 --deconv 0
			fi
		fi
	done
	python /mnt/g/MU_WW/SAM_Refiner/SAM_Refiner.py --collect 0 --chim_rm 0
	rm *covars.tsv
	cd ..
done