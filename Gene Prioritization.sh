# Run magma score

#!/bin/bash
#SBATCH --job-name=magma_sel3
#SBATCH --array=1-22
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --partition=normal
#SBATCH --output=magma_sel3_chr%a.out
#SBATCH --error=magma_sel3_chr%a.err

# 路径设置
magma_bin=/data/circ_ukb/AHua/project/GWAS/regenie/MAGMA/magma
bfile_dir=/data/circ_ukb/AHua/project/GWAS/regenie/MAGMA
sumstats_dir=/data/circ_ukb/AHua/project/GWAS/regenie/height/metaanalysis
gene_annot=${bfile_dir}/MGBB_hg38.genes.annot

chr=${SLURM_ARRAY_TASK_ID}

features=( {add in all features}

)

for feature in "${features[@]}"; do
    sumstats=${sumstats_dir}/height_${feature}_all_metal1.magma.txt
    echo "Processing $feature chr${chr}"

    $magma_bin \
      --bfile ${bfile_dir}/GSA53K.chr${chr}.maf01 \
      --pval $sumstats ncol=N \
      --gene-annot $gene_annot \
      --out ${bfile_dir}/height_${feature}.chr${chr}
done



#Run POPs

#!/bin/bash
#SBATCH --job-name=pops_all
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=pops_all.out
#SBATCH --error=pops_all.err

for feature in \
height_X45_Elongation \
height_X45_Flatness \
height_X45_LeastAxisLength \
height_X45_Maximum2DDiameterColumn \
height_X45_Sphericity \
height_X46_Maximum2DDiameterSlice \
height_X46_Maximum3DDiameter \
height_X46_MeshVolume \
height_X46_MinorAxisLength \
height_X46_Sphericity \
height_X47_MajorAxisLength \
height_X47_Maximum2DDiameterColumn \
height_X47_Sphericity \
height_X47_SurfaceArea \
height_X47_SurfaceVolumeRatio \
height_X48_Flatness \
height_X48_LeastAxisLength \
height_X48_Maximum3DDiameter \
height_X48_MeshVolume \
height_X48_Sphericity \
height_X48_SurfaceVolumeRatio \
height_X48_VoxelVolume \
height_X7_Elongation \
height_X7_LeastAxisLength \
height_X7_Maximum2DDiameterSlice \
height_X7_MinorAxisLength \
height_X7_SurfaceVolumeRatio
do
  echo " Running $feature"
  python /medpop/esp2/AHua/pops/pops.py \
    --gene_annot_path /medpop/esp2/AHua/pops/NCBI38.gene.entrez.tss.reordered.space.loc \
    --feature_mat_prefix /medpop/esp2/AHua/pops/features/pops_features \
    --num_feature_chunks 1 \
    --magma_prefix /medpop/esp2/AHua/pops/gene/${feature}.all \
    --out_prefix /medpop/esp2/AHua/pops/pops_output/${feature}
done
