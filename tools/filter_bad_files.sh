ls -larth | grep "M Sep" | awk '{print $9}' | xargs -I {} mv {} good_files/
