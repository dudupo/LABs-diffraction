




ffmpeg -i $1 -filter_complex "[0:v]eq=contrast=1:brightness=0:saturation=0:gamma=1000:
gamma_r=1:gamma_g=1:gamma_b=1:gamma_weight=1[outv]" -map [outv] $2 


# bio/tif/mgyi_1-20210605T083057Z-004/mgyi_1/Pos0/img_000000000_Phasefast_000.tif
# bio/tif/mgyi_1-20210605T083057Z-004/mgyi_1/Pos0/img_000000000_YFPFast_000.tif
# bio/tif/mgyi_1-20210605T083057Z-004/mgyi_1/Pos0/img_000000166_YFPFast_000.tif

# ffmpeg -i tif/img_000000047_YFPFast_000.tif -filter_complex "extractplanes=y+u+v[y][u][v];   \
# [y]histeq=strength=0.3:intensity=1[lumaeq];   \
# [lumaeq][u][v]mergeplanes=0x001020:yuv420p[out]" -map "[out]" out2.tif 

# blur()
# {
#     ifn=$1
#     #
#     sigma=${sigma:-13}
#     steps=${steps:-4}
#     #
#     ffmpeg  -y -i "${ifn}" -filter_complex "
#     gblur=sigma=${sigma}:steps=${steps}
#     "  -an "$ifn.blur.png"
# }
#     #-map '[v]'
#     # ,drawbox=c=blue[b];

#     # [0v2][b]overlay=30:580[v]

# max=10
# for i in `seq 1 $max`
# do
#     echo y | ffmpeg -i p$i.jpeg -f lavfi -i color=0xBBBBBB:s=1600x1200 -f lavfi -i color=black:s=1600x1200 -f lavfi -i color=white:s=1600x1200 -filter_complex threshold out.p$i.png
#     echo y | ffmpeg -i p$i.jpeg -f lavfi -i color=0x111111:s=1600x1200 -f lavfi -i color=black:s=1600x1200 -f lavfi -i color=white:s=1600x1200 -filter_complex threshold cen.p$i.png
#     echo y | blur out.p$i.png
#     echo y | ffmpeg -i out.p$i.png.blur.png -f lavfi -i color=0x999999:s=1600x1200 -f lavfi -i color=black:s=1600x1200 -f lavfi -i color=white:s=1600x1200 -filter_complex threshold final.out.p$i.png

    
# done 
