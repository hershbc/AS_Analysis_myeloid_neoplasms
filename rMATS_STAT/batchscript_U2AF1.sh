#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hershbc@ccf.org
#SBATCH --job-name=rmats_stat_U2AF1
#SBATCH -N1
#SBATCH --exclusive
module load python/2.7.13
srun  mkdir output/AML_U2AF1_Q157_Hotspot_versus_healthy_bm_pvalue
srun  python2.7 rMATS_unpaired.py AML_U2AF1_Q157_Hotspot_versus_healthy_bm.txt output/AML_U2AF1_Q157_Hotspot_versus_healthy_bm_pvalue/ 20 0.05 > output/AML_U2AF1_Q157_Hotspot_versus_healthy_bm_pvalue/log.txt
srun  mkdir output/AML_U2AF1_S34_Hotspot_versus_healthy_bm_pvalue
srun  python2.7 rMATS_unpaired.py AML_U2AF1_S34_Hotspot_versus_healthy_bm.txt output/AML_U2AF1_S34_Hotspot_versus_healthy_bm_pvalue/ 20 0.05 > output/AML_U2AF1_S34_Hotspot_versus_healthy_bm_pvalue/log.txt
srun  mkdir output/CMML_U2AF1_Q157_Hotspot_versus_healthy_bm_pvalue
srun  python2.7 rMATS_unpaired.py CMML_U2AF1_Q157_Hotspot_versus_healthy_bm.txt output/CMML_U2AF1_Q157_Hotspot_versus_healthy_bm_pvalue/ 20 0.05 > output/CMML_U2AF1_Q157_Hotspot_versus_healthy_bm_pvalue/log.txt
srun  mkdir output/CMML_U2AF1_S34_Hotspot_versus_healthy_bm_pvalue
srun  python2.7 rMATS_unpaired.py CMML_U2AF1_S34_Hotspot_versus_healthy_bm.txt output/CMML_U2AF1_S34_Hotspot_versus_healthy_bm_pvalue/ 20 0.05 > output/CMML_U2AF1_S34_Hotspot_versus_healthy_bm_pvalue/log.txt
srun  mkdir output/MDS_U2AF1_Q157_Hotspot_versus_healthy_bm_pvalue
srun  python2.7 rMATS_unpaired.py MDS_U2AF1_Q157_Hotspot_versus_healthy_bm.txt output/MDS_U2AF1_Q157_Hotspot_versus_healthy_bm_pvalue/ 20 0.05 > output/MDS_U2AF1_Q157_Hotspot_versus_healthy_bm_pvalue/log.txt
srun  mkdir output/MDS_U2AF1_S34_Hotspot_versus_healthy_bm_pvalue
srun  python2.7 rMATS_unpaired.py MDS_U2AF1_S34_Hotspot_versus_healthy_bm.txt output/MDS_U2AF1_S34_Hotspot_versus_healthy_bm_pvalue/ 20 0.05 > output/MDS_U2AF1_S34_Hotspot_versus_healthy_bm_pvalue/log.txt
srun  mkdir output/MDS-MPN-U_U2AF1_Q157_Hotspot_versus_healthy_bm_pvalue
srun  python2.7 rMATS_unpaired.py MDS-MPN-U_U2AF1_Q157_Hotspot_versus_healthy_bm.txt output/MDS-MPN-U_U2AF1_Q157_Hotspot_versus_healthy_bm_pvalue/ 20 0.05 > output/MDS-MPN-U_U2AF1_Q157_Hotspot_versus_healthy_bm_pvalue/log.txt
