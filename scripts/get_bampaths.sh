for dir in ERR*/; do
        filename=$(find ${dir} -maxdepth 1 -name "*.bam" -print -quit | xargs basename);
        path=${dir}${filename};
        echo $path >> bampaths.txt;
done
