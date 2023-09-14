for file in *; do mv $file ${file:0:1}0${file:1:3}0${file:4}; done
