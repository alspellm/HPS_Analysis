for file in *; do mv $file "$(cut -d'-' -f1 <<<${file})_$(cut -d'-' -f2 <<<${file})"; done
