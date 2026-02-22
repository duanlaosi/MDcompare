#!/bin/bash
STRUCTURE_FILE=$1
StructureDir=$2
Color_ID=$3

mkdir -p "$StructureDir"
BASE_FILENAME=$(basename "$STRUCTURE_FILE" | cut -d '.' -f 1)
vmd -dispdev text -eofexit -args "$STRUCTURE_FILE" "$StructureDir" "$Color_ID"<< EOF
# get argv 0 1 2
set structure_file [lindex \$argv 0]
set structure_dir [lindex \$argv 1]
set color_ID [lindex \$argv 2]


# get structure file
mol new \$structure_file type gro waitfor all

# calculatee secondary structure
mol ssrecalc top

mol modstyle 0 0 NewCartoon
mol modcolor 0 0 ColorID \$color_ID
mol modselect 0 0 protein
display projection orthographic
color Display Background white
display backgroundmode transparent
# enhance visualization effects
display depthcue on
display shadows on
display ambientocclusion on
axes location off

# get the base filename
set base_filename [file rootname [file tail \$structure_file]]
set output_image_path [format "%s/%s.tga" \$structure_dir \$base_filename]
render TachyonInternal \$output_image_path

EOF

if [ ! -f "$StructureDir/$BASE_FILENAME.tga" ]; then
    echo "Error: No image generated: ${StructureDir}/${BASE_FILENAME}.tga. Please check the VMD script execution."
    exit 1
fi

ffmpeg -i "$StructureDir/$BASE_FILENAME.tga" "$StructureDir/$BASE_FILENAME.png"

echo "Image saved to $StructureDir/$BASE_FILENAME.png"
