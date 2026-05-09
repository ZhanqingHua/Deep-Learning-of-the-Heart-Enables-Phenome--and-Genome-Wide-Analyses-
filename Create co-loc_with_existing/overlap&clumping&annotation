#find out overlap SNPs

#!/bin/bash
cd /data/circ_ukb/AHua/project/GWAS/regenie/height/metaanalysis

# reference file
ref_file="/data/circ_ukb/AHua/project/GWAS/regenie/comparegwas/categorical-20002-1075_hg38_withID.filtered_ids.sorted.txt"
sort "$ref_file" -o "$ref_file"


features=(
X7_LeastAxisLength
X7_Maximum2DDiameterSlice
X7_SurfaceVolumeRatio
X45_LeastAxisLength
X45_Sphericity
X46_MeshVolume
X47_Flatness
X48_Flatness
X48_MeshVolume
X46_MinorAxisLength
X46_Sphericity
X47_MajorAxisLength
X47_Sphericity
X48_LeastAxisLength
X48_SurfaceVolumeRatio
X7_Elongation
X7_MinorAxisLength
X45_Elongation
X45_Flatness
X45_Maximum2DDiameterColumn
X46_Maximum2DDiameterSlice
X46_Maximum3DDiameter
X47_Maximum2DDiameterColumn
X47_SurfaceArea
X47_SurfaceVolumeRatio
X48_Maximum3DDiameter
X48_Sphericity
X48_VoxelVolume
)


for feature in "${features[@]}"; do
  input="height_${feature}_p_lt5e7.txt"
  snplist="height_${feature}_p_lt5e7.snps.txt"
  sortedlist="height_${feature}_p_lt5e7.sorted.snps.txt"
  overlap="height_${feature}_overlap_with_lifted.txt"

  if [[ ! -f "$input" ]]; then
    echo "skip"
    continue
  fi


  tail -n +2 "$input" | cut -f1 > "$snplist"
  sort "$snplist" -o "$sortedlist"
  comm -12 "$sortedlist" "$ref_file" > "$overlap"

  echo "$feature: $(wc -l < "$overlap") overlapping SNPs"
done

# 2. clumping to lead SNPs
awk 'BEGIN{OFS="\t"} {$2=$1":"$4; print}' g1000_all.hg38.bim > g1000_all.hg38.chrpos.bim
cp g1000_all.hg38.bed g1000_all.hg38.chrpos.bed
cp g1000_all.hg38.fam g1000_all.hg38.chrpos.fam
awk 'NR>1{print $3":"$11}' with_POS_hg38.txt > extract_chrpos.txt
awk 'BEGIN{OFS="\t"; print "SNP","P"} NR>1{print $3":"$11,$7}' with_POS_hg38.txt > sumstats.txt

 /medpop/esp2/btruong/Tools/plink \
>   --bfile my_snps_nodup \
>   --clump sumstats.txt \
>   --clump-snp-field SNP \
>   --clump-field P \
>   --clump-p1 5e-8 \
>   --clump-p2 0.05 \
>   --clump-r2 0.1 \
>   --clump-kb 250 \
>   --out my_clump_r2_01

# 3. Find the nearest gene
awk '
BEGIN{FS=OFS="\t"}
NR==FNR{
  if(FNR>1){
    chr=$3
    gene_chr[chr,++n[chr]]=$0
  }
  next
}
{
  snp_chr=$1
  snp_pos=$2
  snp_id=$3

  best=""
  min_dist=""

  for(i=1;i<=n[snp_chr];i++){
    split(gene_chr[snp_chr,i],a,OFS)
    start=a[4]
    end=a[5]
    tss=a[6]

    if(end >= snp_pos-500000 && start <= snp_pos+500000){
      dist=tss-snp_pos
      if(dist<0){dist=-dist}

      if(min_dist=="" || dist<min_dist){
        min_dist=dist
        best=gene_chr[snp_chr,i]
      }
    }
  }

  if(best==""){
    print snp_chr, snp_pos, snp_id, "NA", "NA", "NA", "NA", "NA", "NA"
  } else {
    split(best,a,OFS)
    print snp_chr, snp_pos, snp_id, a[2], a[1], a[4], a[5], a[6], min_dist
  }
}
' /medpop/esp2/btruong/Tools/gene_annot_hg38.txt lead30_chrpos.txt > lead30_nearest_gene_500kb_hg38.txt
