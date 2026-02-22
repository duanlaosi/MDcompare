#!/bin/bash
STRUCTURE_FILE=$1
TRAJECTORY_FILE=$2
SAVE_TIME=$3
ColorID=$4
MovieDir=$5

mkdir -p "$MovieDir"
vmd -dispdev text -eofexit -args "$STRUCTURE_FILE" "$TRAJECTORY_FILE" "$SAVE_TIME" "$ColorID" "$MovieDir" << EOF

# get argv 0 1 2 3 4
set structure_file [lindex \$argv 0]
set trajectory_file [lindex \$argv 1]
set save_time [lindex \$argv 2]
set colorID [lindex \$argv 3]
set movie_dir [lindex \$argv 4]


mol new \$structure_file type gro waitfor all

mol addfile \$trajectory_file type xtc waitfor all



mol modstyle 0 0 NewCartoon
mol modcolor 0 0 ColorID \$colorID
display projection orthographic
color Display Background white
display backgroundmode transparent
display depthcue on
display shadows on
display ambientocclusion on
axes location off

set step 10
set framerate 25
set nf [molinfo top get numframes]
set frame_idx 0
for {set i 0} {\$i < \$nf} {incr i \$step} {
    animate goto \$i
    set fname [format "%s/frame_%04d.tga" \$movie_dir \$frame_idx]
    render TachyonInternal \$fname
    incr frame_idx
}

EOF

if [ ! -f "$MovieDir/frame_0000.tga" ]; then
  echo "Error: No frames generated in $MovieDir. Please check the VMD script execution."
  exit 1
fi

ffmpeg -y -framerate 25 -i "$MovieDir/frame_%04d.tga" -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -preset medium -crf 23 -pix_fmt yuv420p "$MovieDir/traj_movie_${SAVE_TIME}.mp4"

rm $MovieDir/frame_*.tga

echo "Movie saved to $MovieDir/traj_movie_${SAVE_TIME}.mp4"
